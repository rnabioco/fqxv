//! The never-worse selection rule, stated once (issue #203).
//!
//! Several codec paths choose between two encodings of the *same* data and keep
//! whichever is smaller — the per-block overlap codec vs order-k, the level-8+
//! hashed tier vs plain order-k, the shared long-read reference vs the plain
//! layout, and single-end reorder vs the plain layout. Each such choice used to be
//! written inline, and three of them shipped as bug fixes after a regression was
//! measured (#184, #192, #196): an alternative was adopted that was not actually
//! smaller. Routing every choice through these two helpers makes the rule — and
//! its one hard precondition — impossible to get subtly wrong the next time.
//!
//! **Precondition (the callers' invariant, not checked here): the two candidates
//! MUST decode to identical content.** Only then is "keep the smaller" a size
//! optimization rather than a correctness change. Every caller pairs encodings
//! that are byte-exact round trips of the same reads.
//!
//! **Ties keep the incumbent.** A challenger is adopted only on a *strict* win, so
//! an equal-size alternative — especially one that also stores a fixed cost like a
//! reference frame — is never taken for nothing.

/// Keep the smaller of two already-coded candidates that decode to identical
/// content; a tie keeps `incumbent`. See the module docs for the precondition.
pub(crate) fn keep_smaller(challenger: Vec<u8>, incumbent: Vec<u8>) -> Vec<u8> {
    if challenger.len() < incumbent.len() {
        challenger
    } else {
        incumbent
    }
}

/// Whether to adopt `candidate` over `baseline` when `candidate` also incurs a
/// one-time `fixed_cost` that `baseline` does not — a stored reference frame, a
/// permutation stream. The byte-count companion of [`keep_smaller`] for whole-file
/// layout decisions, where the caller acts on the boolean rather than swapping a
/// buffer. Strict, so the fixed cost is never paid without a net win.
pub(crate) fn adopt_over(candidate: usize, fixed_cost: usize, baseline: usize) -> bool {
    candidate + fixed_cost < baseline
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keep_smaller_breaks_ties_toward_the_incumbent() {
        assert_eq!(
            keep_smaller(vec![1, 2], vec![3]),
            vec![3],
            "challenger larger"
        );
        assert_eq!(
            keep_smaller(vec![1], vec![2, 3]),
            vec![1],
            "challenger smaller"
        );
        // Equal length: the incumbent stays, so nothing is swapped for no gain.
        assert_eq!(
            keep_smaller(vec![9, 9], vec![1, 2]),
            vec![1, 2],
            "tie keeps incumbent"
        );
    }

    #[test]
    fn adopt_over_requires_a_strict_win_including_the_fixed_cost() {
        assert!(adopt_over(10, 0, 11), "smaller with no fixed cost");
        assert!(!adopt_over(10, 0, 10), "a tie is not a win");
        assert!(!adopt_over(10, 5, 12), "fixed cost erases the saving");
        assert!(
            adopt_over(10, 1, 12),
            "still wins once the fixed cost is counted"
        );
    }
}
