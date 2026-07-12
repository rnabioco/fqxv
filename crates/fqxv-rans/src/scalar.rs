//! Portable scalar rANS coder — the correctness reference.
//!
//! Order-0 and order-1, 32 interleaved states, 16-bit word renormalization
//! (the "Nx16" design). The 32 states are independent, which is exactly what
//! lets the SSE4.2/AVX2 backends map them onto vector lanes; this scalar path
//! defines the byte-exact format they must match.
//!
//! Stream layout (little-endian):
//! ```text
//! [u8 order][u64 n]                       // n = number of input bytes
//! if n > 0 and order == 0:
//!   [256 x u16 freq]                      // normalized order-0 table
//!   [32  x u32 state]                     // encoder's final states
//!   [u16 renorm words ...]                // consumed front-to-back on decode
//! if n > 0 and order == 1:
//!   [32 byte context-presence bitmap]
//!   [256 x u16 freq] per present context  // ascending context order
//!   [32  x u32 state]
//!   [u16 renorm words ...]
//! ```

use crate::model::{EncSym, Model, N_STATES, RANS_L, SCALE_BITS, TOTFREQ};
use crate::{Error, Result};

const ORDER0: u8 = 0;
const ORDER1: u8 = 1;

/// Encode `src` with an order-0 model using the scalar coder.
pub(crate) fn encode_order0(src: &[u8]) -> Vec<u8> {
    let mut counts = [0u32; 256];
    for &b in src {
        counts[b as usize] += 1;
    }
    let model = Model::from_counts(&counts);
    let enc = model.enc_table();

    let mut states = [RANS_L; N_STATES];
    let mut renorm: Vec<u16> = Vec::with_capacity(src.len());
    for i in (0..src.len()).rev() {
        let s = src[i] as usize;
        let x = &mut states[i % N_STATES];
        encode_symbol(x, &enc[s], &mut renorm);
    }

    let mut out = Vec::new();
    out.push(ORDER0);
    out.extend_from_slice(&(src.len() as u64).to_le_bytes());
    if !src.is_empty() {
        for s in 0..256 {
            out.extend_from_slice(&model.freq[s].to_le_bytes());
        }
        finish(&mut out, &states, &mut renorm);
    }
    out
}

/// Encode `src` with an order-1 model (context = previous byte).
///
/// One normalized model per occurring previous-byte context; a 32-byte presence
/// bitmap lets us store tables only for contexts that occur.
pub(crate) fn encode_order1(src: &[u8]) -> Vec<u8> {
    let n = src.len();
    let mut counts = vec![[0u32; 256]; 256];
    let mut prev = 0u8;
    for &b in src {
        counts[prev as usize][b as usize] += 1;
        prev = b;
    }

    let mut models: Vec<Option<Model>> = (0..256).map(|_| None).collect();
    // Division-free encode constants per present context (boxed so the 8 KiB
    // table per context doesn't bloat the `Option` layout).
    let mut enc: Vec<Option<Box<[EncSym; 256]>>> = (0..256).map(|_| None).collect();
    for ctx in 0..256 {
        if counts[ctx].iter().any(|&c| c != 0) {
            let m = Model::from_counts(&counts[ctx]);
            enc[ctx] = Some(Box::new(m.enc_table()));
            models[ctx] = Some(m);
        }
    }

    let mut states = [RANS_L; N_STATES];
    let mut renorm: Vec<u16> = Vec::with_capacity(n);
    for i in (0..n).rev() {
        let ctx = if i == 0 { 0 } else { src[i - 1] } as usize;
        let e = enc[ctx].as_ref().expect("context counted during encode");
        let s = src[i] as usize;
        let x = &mut states[i % N_STATES];
        encode_symbol(x, &e[s], &mut renorm);
    }

    let mut out = Vec::new();
    out.push(ORDER1);
    out.extend_from_slice(&(n as u64).to_le_bytes());
    if n != 0 {
        let mut bitmap = [0u8; 32];
        for ctx in 0..256 {
            if models[ctx].is_some() {
                bitmap[ctx / 8] |= 1 << (ctx % 8);
            }
        }
        out.extend_from_slice(&bitmap);
        for m in models.iter().flatten() {
            for s in 0..256 {
                out.extend_from_slice(&m.freq[s].to_le_bytes());
            }
        }
        finish(&mut out, &states, &mut renorm);
    }
    out
}

/// Encode one symbol into state `x`: renormalize down, then apply the rANS map.
///
/// The map ([`EncSym::apply`]) is reciprocal-based where the state is provably
/// `< 2^31`, exact division otherwise; the precomputed `x_max` bounds
/// renormalization into `[x_max >> 16, x_max)` so the result lands back in
/// `[RANS_L, RANS_L << 16)`. `x_max` can reach 2^32, hence the u64 compare.
#[inline]
fn encode_symbol(x: &mut u32, sym: &EncSym, renorm: &mut Vec<u16>) {
    let mut v = *x;
    while u64::from(v) >= sym.x_max {
        renorm.push((v & 0xffff) as u16);
        v >>= 16;
    }
    *x = sym.apply(v);
}

/// Write the final states and the (reversed) renorm words.
fn finish(out: &mut Vec<u8>, states: &[u32; N_STATES], renorm: &mut [u16]) {
    for x in states {
        out.extend_from_slice(&x.to_le_bytes());
    }
    // Encoder emitted words in reverse symbol order; reversing yields the order
    // the forward decoder consumes.
    renorm.reverse();
    for w in renorm.iter() {
        out.extend_from_slice(&w.to_le_bytes());
    }
}

/// Decode a stream produced by [`encode_order0`] or [`encode_order1`].
pub(crate) fn decode(src: &[u8]) -> Result<Vec<u8>> {
    let mut r = Reader::new(src);
    let order = r.u8()?;
    let n = r.u64()? as usize;
    if n == 0 {
        return Ok(Vec::new());
    }
    match order {
        ORDER0 => decode_order0(&mut r, n),
        ORDER1 => decode_order1(&mut r, n),
        _ => Err(Error::Malformed("unknown order tag")),
    }
}

fn decode_order0(r: &mut Reader<'_>, n: usize) -> Result<Vec<u8>> {
    let mut freq = [0u16; 256];
    for f in &mut freq {
        *f = r.u16()?;
    }
    let model = Model::from_freqs(freq);

    let mut states = read_states(r)?;
    let mut out = vec![0u8; n];
    for (i, o) in out.iter_mut().enumerate() {
        let x = &mut states[i % N_STATES];
        let s = model.slot2sym[(*x & (TOTFREQ - 1)) as usize];
        *o = s;
        step_state(x, &model, s, r)?;
    }
    Ok(out)
}

fn decode_order1(r: &mut Reader<'_>, n: usize) -> Result<Vec<u8>> {
    let mut bitmap = [0u8; 32];
    for b in &mut bitmap {
        *b = r.u8()?;
    }
    let mut models: Vec<Option<Model>> = (0..256).map(|_| None).collect();
    for ctx in 0..256 {
        if (bitmap[ctx / 8] >> (ctx % 8)) & 1 == 1 {
            let mut freq = [0u16; 256];
            for f in &mut freq {
                *f = r.u16()?;
            }
            models[ctx] = Some(Model::from_freqs(freq));
        }
    }

    let mut states = read_states(r)?;
    let mut out = vec![0u8; n];
    let mut prev = 0u8;
    for (i, o) in out.iter_mut().enumerate() {
        let ctx = if i == 0 { 0 } else { prev } as usize;
        let model = models[ctx]
            .as_ref()
            .ok_or(Error::Malformed("symbol references absent context"))?;
        let x = &mut states[i % N_STATES];
        let s = model.slot2sym[(*x & (TOTFREQ - 1)) as usize];
        *o = s;
        step_state(x, model, s, r)?;
        prev = s;
    }
    Ok(out)
}

fn read_states(r: &mut Reader<'_>) -> Result<[u32; N_STATES]> {
    let mut states = [0u32; N_STATES];
    for x in &mut states {
        *x = r.u32()?;
    }
    Ok(states)
}

/// Advance one decoded symbol: reverse the rANS map and renormalize up.
#[inline]
fn step_state(x: &mut u32, model: &Model, s: u8, r: &mut Reader<'_>) -> Result<()> {
    let f = u32::from(model.freq[s as usize]);
    let c = u32::from(model.cum[s as usize]);
    let mut v = f * (*x >> SCALE_BITS) + (*x & (TOTFREQ - 1)) - c;
    while v < RANS_L {
        v = (v << 16) | u32::from(r.u16()?);
    }
    *x = v;
    Ok(())
}

/// Minimal forward reader over the encoded stream. Shared with the SIMD backend.
pub(crate) struct Reader<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    pub(crate) fn new(buf: &'a [u8]) -> Self {
        Reader { buf, pos: 0 }
    }
    /// Byte offset of the next unread byte (start of the renorm region once the
    /// header has been parsed).
    pub(crate) fn pos(&self) -> usize {
        self.pos
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
    pub(crate) fn u8(&mut self) -> Result<u8> {
        Ok(self.take::<1>()?[0])
    }
    pub(crate) fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.take::<2>()?))
    }
    pub(crate) fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take::<4>()?))
    }
    pub(crate) fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.take::<8>()?))
    }
}
