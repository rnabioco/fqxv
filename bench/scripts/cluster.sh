#!/usr/bin/env bash
# Single source of truth for Slurm submission parameters, sourced by the bench
# drivers (submit_parallel.sh, corpus.sh). Partition/QoS/account used to be
# hardcoded across a dozen .sbatch files that had drifted apart — five defaulted
# to Alpine's `amilan`, three to Bodhi's `rna` — so a default `submit_parallel.sh`
# on Bodhi queued against a partition that does not exist there.
#
# The harness runs on two clusters:
#   bodhi  (default) — `rna` partition, 88+ cores / 754G, no per-core memory cap
#   alpine           — `amilan` partition, hard-caps memory at 3840M/core, and
#                      may need an explicit --account (e.g. amc-general)
#
# Detection probes `sinfo` for the partition, so it normally just works. Override
# any piece:
#   FQXV_CLUSTER=bodhi|alpine                pick the profile, skipping detection
#   FQXV_PARTITION / FQXV_QOS / FQXV_ACCOUNT override one field directly
#
# Provides:
#   FQXV_SBATCH_OPTS  array of flags to splat into a submission:
#                       sbatch "${FQXV_SBATCH_OPTS[@]}" some.sbatch
# and exports SBATCH_PARTITION / SBATCH_QOS / SBATCH_ACCOUNT so that a plain
# `sbatch slurm/foo.sbatch` from the same shell also lands on the right
# partition — Slurm environment variables override a script's own #SBATCH
# directives, which is what makes the .sbatch defaults safe to keep as a
# fallback for direct submission.
#
# Note: every statement below must succeed, because callers source this under
# `set -e`. Hence `if` blocks rather than `cond && action`, whose non-zero
# "false" would abort the caller.

fqxv_has_partition() {  # name -> 0 if that partition exists on this cluster
  sinfo -h -o '%P' 2>/dev/null | sed 's/\*$//' | grep -qx "$1"
}

fqxv_cluster() {  # -> echoes the resolved cluster name
  if [[ -n "${FQXV_CLUSTER:-}" ]]; then echo "$FQXV_CLUSTER"; return 0; fi
  if fqxv_has_partition rna; then echo bodhi; return 0; fi
  if fqxv_has_partition amilan; then echo alpine; return 0; fi
  echo bodhi   # neither probe matched (no sinfo?); an explicit FQXV_PARTITION still wins
}

FQXV_CLUSTER_RESOLVED="$(fqxv_cluster)"
case "$FQXV_CLUSTER_RESOLVED" in
  alpine) _fqxv_def_partition=amilan ;;
  *)      _fqxv_def_partition=rna ;;
esac

FQXV_PARTITION="${FQXV_PARTITION:-$_fqxv_def_partition}"
FQXV_QOS="${FQXV_QOS:-normal}"
# Empty by default: on both clusters the user's default association works (this
# is how the 2026-07-20 Bodhi run was submitted). Alpine associations are
# sometimes unset — export FQXV_ACCOUNT=amc-general there if sbatch complains.
FQXV_ACCOUNT="${FQXV_ACCOUNT:-}"

FQXV_SBATCH_OPTS=( --partition="$FQXV_PARTITION" --qos="$FQXV_QOS" )
export SBATCH_PARTITION="$FQXV_PARTITION"
export SBATCH_QOS="$FQXV_QOS"
if [[ -n "$FQXV_ACCOUNT" ]]; then
  FQXV_SBATCH_OPTS+=( --account="$FQXV_ACCOUNT" )
  export SBATCH_ACCOUNT="$FQXV_ACCOUNT"
fi

unset _fqxv_def_partition
