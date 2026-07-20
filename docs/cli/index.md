# CLI Reference

The `fqxv` binary has four subcommands:

| Command | Purpose |
| --- | --- |
| [`compress`](compress.md) | Compress one or more FASTQ files to a `.fqxv` archive |
| [`decompress`](decompress.md) | Restore FASTQ — to files, split mates, or stdout |
| [`verify`](verify.md) | Validate a container's CRC-32C checksums without decompressing |
| [`info`](info.md) | Print container metadata and per-stream sizes |

![fqxv --help](../images/help.gif)

## Global options

| Option | Description |
| --- | --- |
| `--threads <N>` | Worker threads (0 = all available cores). Default: 16, capped at available cores. |
| `-v, --verbose` | Increase log verbosity (`-v` debug, `-vv`/`-vvv` trace); overridden by `RUST_LOG`. |
| `-q, --quiet` | Silence all output except warnings and errors (also suppresses the summary). |
| `-h, --help` | Print help. |
| `-V, --version` | Print version. A build that is not a clean release tag also reports its git description (`0.3.0 (v0.3.0-7-gab12cd34-dirty)`). |

`compress` and `decompress` fan blocks out across threads with `rayon`; the
output is byte-identical regardless of thread count.
