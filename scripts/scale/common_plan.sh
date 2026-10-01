#!/usr/bin/env bash
# Sourced by scale-study planning presets.
set -Eeuo pipefail

link_or_copy() {
  local source="$1"
  local destination="$2"
  mkdir -p "$(dirname "$destination")"
  rm -f "$destination"
  if ln "$source" "$destination" 2>/dev/null; then
    return 0
  fi
  cp "$source" "$destination"
}

stage_scale_reference() {
  local reference_root="$1"
  local study_root="$2"
  mkdir -p "$study_root"/{run_input_truth,run_realized_truth,run_stats,run_provenance}
  local name
  for name in \
    OpenSwathTimSim.high_signal_selection.tsv \
    OpenSwathTimSim.reference_truth.tsv \
    OpenSwathTimSim.target_only.transitions.tsv \
    OpenSwathTimSim.transitions.tsv; do
    [[ -s "$reference_root/$name" ]] || { echo "ERROR: scale reference artifact missing: $reference_root/$name" >&2; return 1; }
    link_or_copy "$reference_root/$name" "$study_root/$name"
  done
  if [[ -s "$reference_root/OpenSwathTimSim.entrapment_truth.tsv" ]]; then
    link_or_copy "$reference_root/OpenSwathTimSim.entrapment_truth.tsv" "$study_root/OpenSwathTimSim.entrapment_truth.tsv"
  fi
}
