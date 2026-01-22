#!/usr/bin/env bash
set -euo pipefail

DEFAULT_BASE="/eos/cms/store/group/phys_smp/ec/zhye/pp6j_15GeV_MINIAOD"

usage() {
  echo "Usage:"
  echo "  $0 <partNumber> [pattern]"
  echo "  $0 </path/to/PartX> [pattern]"
  echo "Examples:"
  echo "  $0 5"
  echo "  $0 5 'JME-RunIISummer20UL17MiniAODv2-*.root'"
  echo "  $0 /eos/cms/store/group/phys_smp/ec/zhye/pp6j_15GeV_MINIAOD/Part5"
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage
  exit 1
fi

ARG="$1"
PATTERN="${2:-*.root}"

if [[ "$ARG" =~ ^[0-9]+$ ]]; then
  PART_DIR="${DEFAULT_BASE}/Part${ARG}"
else
  PART_DIR="$ARG"
fi

if [[ ! -d "$PART_DIR" ]]; then
  echo "Error: directory not found: $PART_DIR" >&2
  exit 1
fi

mapfile -t FILES < <(find "$PART_DIR" -maxdepth 1 -type f -name "$PATTERN" | sort)

if [[ ${#FILES[@]} -eq 0 ]]; then
  echo "No files matched in $PART_DIR with pattern $PATTERN" >&2
  exit 1
fi

total=0
count=0

read_events() {
  local path="$1"
  local out=""

  out=$( { das-cmssw el8 edmFileUtil -P -f "$path" 2>/dev/null | awk '/^Events/{print $2}'; } || true )
  if [[ -z "$out" && "$path" == /eos/* ]]; then
    local xrootd="root://eosuser.cern.ch/${path}"
    out=$( { das-cmssw el8 edmFileUtil -P -f "$xrootd" 2>/dev/null | awk '/^Events/{print $2}'; } || true )
  fi

  echo "$out"
}

for f in "${FILES[@]}"; do
  # edmFileUtil prints a line like: "Events: 12345"
  n=$(read_events "$f")
  if [[ -z "${n:-}" ]]; then
    echo "Warning: could not read Events from $f" >&2
    continue
  fi
  total=$((total + n))
  count=$((count + 1))
done

echo "Files counted: $count"
echo "Total events: $total"
