//! Emit the linker arguments an extension module needs.
//!
//! An extension module must not link libpython — the interpreter that loads it
//! already supplies those symbols. `pyo3/extension-module` suppresses the link
//! line, but nothing in PyO3 tells the linker to expect the resulting undefined
//! symbols; `pyo3-ffi`'s own build script does not, and the helper below is
//! documented as the downstream crate's job.
//!
//! On ELF targets this is a no-op: a shared object may carry undefined symbols,
//! so the link succeeds either way. Mach-O resolves everything at link time, so
//! without `-undefined dynamic_lookup` every `_Py_*` reference is a hard error
//! and `cargo build` fails on macOS while CI (Linux) stays green.
//!
//! `add_extension_module_link_args` emits `rustc-cdylib-link-arg`, so the flags
//! reach the `_fqxv` cdylib and nothing else — no binary or test target has its
//! undefined-symbol checking weakened.

fn main() {
    pyo3_build_config::add_extension_module_link_args();
}
