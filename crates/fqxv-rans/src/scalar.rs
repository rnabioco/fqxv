//! Portable scalar rANS coder — the correctness reference.
//!
//! Order-0, 32 interleaved states, byte-wise renormalization. The 32 states are
//! independent, which is exactly what lets the SSE4.2/AVX2 backends map them
//! onto vector lanes later; this scalar path defines the byte-exact format they
//! must match.
//!
//! Stream layout (little-endian):
//! ```text
//! [u8 order][u64 n]                       // n = number of input bytes
//! if n > 0:
//!   [256 x u16 freq]                      // normalized order-0 table
//!   [32  x u32 state]                     // encoder's final states
//!   [renorm bytes ...]                    // consumed front-to-back on decode
//! ```

use crate::model::{Model, N_STATES, RANS_L, SCALE_BITS, TOTFREQ};
use crate::{Error, Result};

const ORDER0: u8 = 0;

/// Encode `src` with an order-0 model using the scalar coder.
pub(crate) fn encode_order0(src: &[u8]) -> Vec<u8> {
    let mut counts = [0u32; 256];
    for &b in src {
        counts[b as usize] += 1;
    }
    let model = Model::from_counts(&counts);

    // Encode symbols in reverse; state `i % N_STATES` handles symbol `i`.
    let mut states = [RANS_L; N_STATES];
    let mut renorm: Vec<u8> = Vec::new();
    for i in (0..src.len()).rev() {
        let s = src[i] as usize;
        let f = u32::from(model.freq[s]);
        let c = u32::from(model.cum[s]);
        let x = &mut states[i % N_STATES];
        // Renormalize down into [x_max/256, x_max) so the C-step lands back in
        // the normalized interval.
        let x_max = ((RANS_L >> SCALE_BITS) << 8) * f;
        let mut v = *x;
        while v >= x_max {
            renorm.push((v & 0xff) as u8);
            v >>= 8;
        }
        *x = ((v / f) << SCALE_BITS) + (v % f) + c;
    }

    let mut out = Vec::with_capacity(1 + 8 + 512 + 128 + renorm.len());
    out.push(ORDER0);
    out.extend_from_slice(&(src.len() as u64).to_le_bytes());
    if !src.is_empty() {
        for s in 0..256 {
            out.extend_from_slice(&model.freq[s].to_le_bytes());
        }
        for x in states {
            out.extend_from_slice(&x.to_le_bytes());
        }
        // Encoder emitted bytes in reverse symbol order; reversing the whole
        // sequence yields the order the forward decoder consumes.
        renorm.reverse();
        out.extend_from_slice(&renorm);
    }
    out
}

/// Decode a stream produced by [`encode_order0`].
pub(crate) fn decode(src: &[u8]) -> Result<Vec<u8>> {
    let mut r = Reader::new(src);
    let order = r.u8()?;
    let n = r.u64()? as usize;
    if n == 0 {
        return Ok(Vec::new());
    }
    if order != ORDER0 {
        return Err(Error::NotImplemented("rANS order-1 decode (M1 follow-up)"));
    }

    let mut freq = [0u16; 256];
    for f in &mut freq {
        *f = r.u16()?;
    }
    let model = Model::from_freqs(freq);

    let mut states = [0u32; N_STATES];
    for x in &mut states {
        *x = r.u32()?;
    }

    let mut out = vec![0u8; n];
    for (i, o) in out.iter_mut().enumerate() {
        let x = &mut states[i % N_STATES];
        let slot = (*x & (TOTFREQ - 1)) as usize;
        let s = model.slot2sym[slot];
        *o = s;
        let f = u32::from(model.freq[s as usize]);
        let c = u32::from(model.cum[s as usize]);
        let mut v = f * (*x >> SCALE_BITS) + (*x & (TOTFREQ - 1)) - c;
        while v < RANS_L {
            v = (v << 8) | u32::from(r.u8()?);
        }
        *x = v;
    }
    Ok(out)
}

/// Minimal forward reader over the encoded stream.
struct Reader<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn new(buf: &'a [u8]) -> Self {
        Reader { buf, pos: 0 }
    }
    fn take<const N: usize>(&mut self) -> Result<[u8; N]> {
        let end = self.pos + N;
        let slice = self
            .buf
            .get(self.pos..end)
            .ok_or(Error::Malformed("truncated stream"))?;
        self.pos = end;
        Ok(slice.try_into().expect("slice length checked"))
    }
    fn u8(&mut self) -> Result<u8> {
        Ok(self.take::<1>()?[0])
    }
    fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.take::<2>()?))
    }
    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take::<4>()?))
    }
    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.take::<8>()?))
    }
}
