# CLI Reference

The `fqxv` binary has four subcommands:

| Command | Purpose |
| --- | --- |
| [`compress`](compress.md) | Compress one or more FASTQ files to a `.fqxv` archive |
| [`decompress`](decompress.md) | Restore FASTQ — to files, split mates, or stdout |
| [`verify`](verify.md) | Validate a container's CRC-32C checksums without decompressing |
| [`info`](info.md) | Print container metadata and per-stream sizes |

`verify` and `info` also take several archives, or a directory scanned
recursively for `*.fqxv`, and report one entry per archive.

![fqxv --help](../images/help.gif)

## Global options

| Option | Description |
| --- | --- |
| `--threads <N>` | Worker threads (0 = all available cores). Default: 16, capped at available cores. |
| `-v, --verbose` | Increase log verbosity, repeatable (`-v` info, `-vv` debug, `-vvv` trace with targets, thread ids, and span timing); overridden by `RUST_LOG`. Logs go to stderr, so piped FASTQ stays clean. |
| `-q, --quiet` | Silence all output except warnings and errors (also suppresses the progress indicator and the summary). |
| `-h, --help` | Print help (`-h` for the summary, `--help` for the long form). |
| `-V, --version` | Print version. A build that is not a clean release tag also reports its git description (`0.7.0 (v0.7.0-5-gb2c6fee)`). |

`compress` and `decompress` fan blocks out across threads with `rayon`; the
output is byte-identical regardless of thread count.
