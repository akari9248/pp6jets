#!/usr/bin/env bash
set -euo pipefail

DEFAULT_BASE_MINIAOD="/eos/cms/store/group/phys_smp/ec/zhye/pp6j_20GeV_MINIAOD"
DEFAULT_BASE_LHE="/eos/cms/store/group/phys_smp/ec/zhye/pp6j_20GeV"

usage() {
  echo "Usage:"
  echo "  $0 <partNumber> [pattern]"
  echo "  $0 </path/to/PartX> [pattern]"
  echo ""
  echo "Examples:"
  echo "  $0 5                           # MINIAOD (default *.root)"
  echo "  $0 5 '*.lhe'                   # LHE files"
  echo "  $0 /path/to/Part5 '*.root'"
}

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage
  exit 1
fi

ARG="$1"
PATTERN="${2:-*.root}"

# 判断文件类型
if [[ "$PATTERN" == *".lhe"* ]]; then
  FILE_TYPE="lhe"
  DEFAULT_BASE="$DEFAULT_BASE_LHE"
else
  FILE_TYPE="root"
  DEFAULT_BASE="$DEFAULT_BASE_MINIAOD"
fi

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

echo "Processing ${#FILES[@]} $FILE_TYPE files from $PART_DIR ..."

total=0
count=0
failed=0

if [[ "$FILE_TYPE" == "lhe" ]]; then
  idx=0
  nfiles=${#FILES[@]}
  for f in "${FILES[@]}"; do
    ((idx++)) || true
    printf "\r[%d/%d] %d events so far..." "$idx" "$nfiles" "$total" >&2
    line=$(grep -m1 "! unwtd events" "$f" 2>/dev/null || true)
    if [[ "$line" =~ ^[[:space:]]*([0-9]+) ]]; then
      n="${BASH_REMATCH[1]}"
      total=$((total + n))
      ((count++)) || true
    else
      ((failed++)) || true
    fi
  done
  echo "" >&2  # 换行
else
  output=$(das-cmssw el8 edmFileUtil -P -f "${FILES[@]}" 2>/dev/null || true)
  while IFS= read -r line; do
    if [[ "$line" =~ ([0-9]+)\ events ]]; then
      n="${BASH_REMATCH[1]}"
      total=$((total + n))
      ((count++)) || true
    fi
  done <<< "$output"
  failed=$((${#FILES[@]} - count))
fi

echo "Files counted: $count / ${#FILES[@]}"
[[ $failed -gt 0 ]] && echo "Failed files: $failed"
echo "Total events: $total"
