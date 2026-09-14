#!/usr/bin/env bash
set -euo pipefail

PRODUCTION_BASE="/eos/cms/store/group/phys_smp/ec/zhye/ALPGEN/Run2026C_13p6TeV/pp6j_25GeV"
if [[ "$(id -un)" != zhye ]]; then
  PRODUCTION_BASE="$PRODUCTION_BASE/$(id -un)"
fi
DEFAULT_BASE_MINIAOD="$PRODUCTION_BASE/MiniAOD_CP5_AQCDUP_v1"
DEFAULT_BASE_LHE="$PRODUCTION_BASE/LHE"

NJOBS="$(nproc 2>/dev/null || echo 8)"
if [[ "$NJOBS" -gt 16 ]]; then
  NJOBS=16
fi

GNU_PARALLEL=""
find_gnu_parallel() {
  local p
  for p in /usr/bin/parallel /bin/parallel; do
    if [[ -x "$p" ]] && "$p" --version 2>&1 | grep -q GNU; then
      GNU_PARALLEL="$p"
      return 0
    fi
  done
  return 1
}

usage() {
  echo "Usage:"
  echo "  $0 [-j N] <partNumber|all|/path/to/PartX> [pattern]"
  echo ""
  echo "Examples:"
  echo "  $0 -j 16 all '*.lhe'          # all parts under default LHE base"
  echo "  $0 all '*.lhe'                # same, auto job count (<=16)"
  echo "  $0 51 '*.lhe'                 # single part"
  echo "  $0 5                          # MINIAOD single part (*.root)"
  echo "  $0 /path/to/pp6j_25GeV all '*.lhe'"
}

# Count events in one LHE file.
# ALPGEN puts "! unwtd events" in the header comment at the top (not at EOF); scan head only on EOS.
count_lhe_file() {
  local f="$1"
  local line n

  line=$(head -c 131072 "$f" 2>/dev/null | grep -m1 '! unwtd events' || true)
  if [[ "$line" =~ ^[[:space:]]*([0-9]+) ]]; then
    echo "${BASH_REMATCH[1]}"
    return 0
  fi

  n=$(grep -c '<event>' "$f" 2>/dev/null) || n=0
  echo "$n"
  return 0
}

count_lhe_file_tagged() {
  local f="$1"
  local p n
  p=$(basename "$(dirname "$f")")
  n=$(count_lhe_file "$f" || true)
  printf '%s\t%s\n' "$p" "$n"
}

discover_part_dirs() {
  local base="$1"
  find "$base" -maxdepth 1 -mindepth 1 -type d -iname 'part*' -printf '%f\t%p\n' \
    | sort -t$'\t' -k1,1 -V \
    | cut -f2-
}

collect_files() {
  local dir pattern
  dir="$1"
  pattern="$2"
  find "$dir" -maxdepth 1 -type f -name "$pattern" -printf '%p\n' | sort
}

run_parallel_lhe() {
  local filelist="$1"
  local nfiles="$2"
  local tmp results joblog
  tmp=$(mktemp -d)
  results="$tmp/counts.tsv"
  joblog="$tmp/joblog.txt"
  trap 'rm -rf "$tmp"' RETURN

  if ! find_gnu_parallel; then
    echo "Error: GNU parallel not found (need /usr/bin/parallel, not DAS parallel)" >&2
    return 1
  fi

  # Wrapper script so GNU parallel always invokes a real executable (exported -f is flaky with bash -c on some nodes)
  cat >"$tmp/count_one.sh" <<'EOS'
#!/usr/bin/env bash
f="$1"
p=$(basename "$(dirname "$f")")
line=$(head -c 131072 "$f" 2>/dev/null | grep -m1 '! unwtd events' || true)
if [[ "$line" =~ ^[[:space:]]*([0-9]+) ]]; then
  n="${BASH_REMATCH[1]}"
else
  n=$(grep -c '<event>' "$f" 2>/dev/null) || n=0
fi
printf '%s\t%s\n' "$p" "$n"
EOS
  chmod +x "$tmp/count_one.sh"

  echo "Processing $nfiles lhe files with $NJOBS jobs ($GNU_PARALLEL) ..." >&2
  echo "(Tip: run 'parallel --citation' once to silence citation notice)" >&2

  local -a par_opts=(--will-cite -j"$NJOBS" --line-buffer --joblog "$joblog")
  # Progress bar on stderr only when attached to a terminal (stdout goes to $results)
  if [[ -t 2 ]]; then
    par_opts+=(--bar --eta)
  fi

  "$GNU_PARALLEL" "${par_opts[@]}" \
    "$tmp/count_one.sh" :::: "$filelist" >"$results"

  if [[ -s "$joblog" ]]; then
    local nfail
    nfail=$(awk 'NR > 1 && $7 != 0 {c++} END {print c+0}' "$joblog")
    if [[ "$nfail" -gt 0 ]]; then
      echo "Warning: $nfail parallel jobs failed (see $joblog)" >&2
    fi
  fi

  awk -F'\t' '
    {
      part[$1] += $2
      files[$1]++
      total += $2
      if ($2 > 0) ok++
      else fail++
    }
    END {
      print "FILES_COUNTED", ok+0, ok+fail+0
      print "FILES_FAILED", fail+0
      print "TOTAL", total+0
      n = asorti(part, sorted)
      for (i = 1; i <= n; i++) {
        p = sorted[i]
        printf "PART\t%s\t%d\t%d\n", p, part[p], files[p]
      }
    }
  ' "$results"
}

run_serial_lhe() {
  local dir="$1"
  local pattern="$2"
  local -a files=()
  mapfile -t files < <(collect_files "$dir" "$pattern")

  if [[ ${#files[@]} -eq 0 ]]; then
    echo "No files matched in $dir with pattern $pattern" >&2
    exit 1
  fi

  echo "Processing ${#files[@]} lhe files from $dir ..." >&2
  local total=0 count=0 failed=0 idx=0 nfiles line n
  nfiles=${#files[@]}
  for f in "${files[@]}"; do
    ((idx++)) || true
    printf '\r[%d/%d] %d events so far...' "$idx" "$nfiles" "$total" >&2
    n=$(count_lhe_file "$f" || echo 0)
    if [[ "$n" -gt 0 ]]; then
      total=$((total + n))
      ((count++)) || true
    else
      ((failed++)) || true
    fi
  done
  echo "" >&2
  echo "FILES_COUNTED $count ${#files[@]}"
  echo "FILES_FAILED $failed"
  echo "TOTAL $total"
}

run_parallel_root_batch() {
  local filelist="$1"
  local batch_size=100
  local tmp results total=0 count=0 failed=0 nfiles
  tmp=$(mktemp -d)
  results="$tmp/root_counts.txt"
  trap 'rm -rf "$tmp"' RETURN

  nfiles=$(wc -l <"$filelist")
  echo "Processing $nfiles root files (batched edmFileUtil) ..." >&2

  local -a batch=()
  local i=0 nb=0
  while IFS= read -r f; do
    batch+=("$f")
    if [[ ${#batch[@]} -ge $batch_size ]]; then
      ((nb++)) || true
      printf '\r[batch %d] %d / %d files submitted...' "$nb" "$((nb * batch_size))" "$nfiles" >&2
      das-cmssw el8 edmFileUtil -P -f "${batch[@]}" 2>/dev/null >>"$results" || true
      batch=()
    fi
    ((i++)) || true
  done <"$filelist"
  if [[ ${#batch[@]} -gt 0 ]]; then
    ((nb++)) || true
    printf '\r[batch %d] %d / %d files submitted...' "$nb" "$nfiles" "$nfiles" >&2
    das-cmssw el8 edmFileUtil -P -f "${batch[@]}" 2>/dev/null >>"$results" || true
  fi
  echo "" >&2

  while IFS= read -r line; do
    if [[ "$line" =~ ([0-9]+)\ events ]]; then
      n="${BASH_REMATCH[1]}"
      total=$((total + n))
      ((count++)) || true
    fi
  done <"$results"
  failed=$((nfiles - count))

  echo "FILES_COUNTED $count $nfiles"
  echo "FILES_FAILED $failed"
  echo "TOTAL $total"
}

# --- parse args ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    -j)
      NJOBS="$2"
      shift 2
      ;;
    -j*)
      NJOBS="${1#-j}"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      break
      ;;
  esac
done

if [[ $# -lt 1 || $# -gt 2 ]]; then
  usage
  exit 1
fi

ARG="$1"
PATTERN="${2:-*.root}"

if [[ "$PATTERN" == *".lhe"* ]]; then
  FILE_TYPE="lhe"
  DEFAULT_BASE="$DEFAULT_BASE_LHE"
else
  FILE_TYPE="root"
  DEFAULT_BASE="$DEFAULT_BASE_MINIAOD"
fi

MODE="single"
BASE="$DEFAULT_BASE"
PART_DIR=""

if [[ "$ARG" == "all" ]]; then
  MODE="all"
elif [[ -d "$ARG" && "$(basename "$ARG")" =~ ^[Pp]art ]]; then
  PART_DIR="$ARG"
  BASE="$(dirname "$ARG")"
elif [[ -d "$ARG" ]]; then
  if [[ -n "$(find "$ARG" -maxdepth 1 -type f -name "$PATTERN" -print -quit)" ]]; then
    PART_DIR="$ARG"
    BASE="$(dirname "$ARG")"
  else
    BASE="$ARG"
    MODE="all"
  fi
elif [[ "$ARG" =~ ^[0-9]+$ ]]; then
  PART_DIR="${DEFAULT_BASE}/Part${ARG}"
  if [[ ! -d "$PART_DIR" ]]; then
    PART_DIR="${DEFAULT_BASE}/part${ARG}"
  fi
else
  echo "Error: invalid argument: $ARG" >&2
  exit 1
fi

if [[ "$MODE" == "single" ]]; then
  if [[ ! -d "$PART_DIR" ]]; then
    echo "Error: directory not found: $PART_DIR" >&2
    exit 1
  fi
  if [[ "$FILE_TYPE" == "lhe" ]]; then
    if find_gnu_parallel; then
      tmp=$(mktemp -d)
      trap 'rm -rf "$tmp"' EXIT
      filelist="$tmp/files.lst"
      collect_files "$PART_DIR" "$PATTERN" >"$filelist"
      nfiles=$(wc -l <"$filelist")
      [[ "$nfiles" -eq 0 ]] && { echo "No files matched in $PART_DIR with pattern $PATTERN" >&2; exit 1; }
      summary=$(run_parallel_lhe "$filelist" "$nfiles")
    else
      run_serial_lhe "$PART_DIR" "$PATTERN" | {
        read -r _ c t
        read -r _ f _
        read -r _ total _
        echo "Files counted: $c / $t"
        [[ "${f:-0}" -gt 0 ]] && echo "Failed files: $f"
        echo "Total events: $total"
      }
      exit 0
    fi
  else
    tmp=$(mktemp -d)
    trap 'rm -rf "$tmp"' EXIT
    filelist="$tmp/files.lst"
    collect_files "$PART_DIR" "$PATTERN" >"$filelist"
    nfiles=$(wc -l <"$filelist")
    [[ "$nfiles" -eq 0 ]] && { echo "No files matched in $PART_DIR with pattern $PATTERN" >&2; exit 1; }
    echo "Processing $nfiles root files from $PART_DIR ..."
    summary=$(run_parallel_root_batch "$filelist")
  fi
else
  if [[ ! -d "$BASE" ]]; then
    echo "Error: base directory not found: $BASE" >&2
    exit 1
  fi
  mapfile -t PART_DIRS < <(discover_part_dirs "$BASE")
  if [[ ${#PART_DIRS[@]} -eq 0 ]]; then
    echo "Error: no Part* directories under $BASE" >&2
    exit 1
  fi

  echo "Found ${#PART_DIRS[@]} parts under $BASE" >&2
  tmp=$(mktemp -d)
  trap 'rm -rf "$tmp"' EXIT
  filelist="$tmp/files.lst"
  for pdir in "${PART_DIRS[@]}"; do
    collect_files "$pdir" "$PATTERN" >>"$filelist"
  done
  nfiles=$(wc -l <"$filelist")

  if [[ "$nfiles" -eq 0 ]]; then
    echo "No files matched under $BASE with pattern $PATTERN" >&2
    exit 1
  fi

  if [[ "$FILE_TYPE" == "lhe" ]]; then
    if ! find_gnu_parallel; then
      echo "Error: GNU parallel required for 'all' mode (/usr/bin/parallel)" >&2
      exit 1
    fi
    summary=$(run_parallel_lhe "$filelist" "$nfiles")
  else
    echo "Processing $nfiles root files across ${#PART_DIRS[@]} parts ..."
    summary=$(run_parallel_root_batch "$filelist")
  fi
fi

# --- print summary ---
files_counted=0 files_total=0 files_failed=0 total_events=0
declare -A part_events part_files
while IFS= read -r line; do
  case "$line" in
    FILES_COUNTED*)
      read -r _ files_counted files_total <<<"$line"
      ;;
    FILES_FAILED*)
      read -r _ files_failed <<<"$line"
      ;;
    TOTAL*)
      read -r _ total_events <<<"$line"
      ;;
    PART*)
      read -r _ pname pe pf <<<"$line"
      part_events["$pname"]=$pe
      part_files["$pname"]=$pf
      ;;
  esac
done <<<"$summary"

if [[ "$MODE" == "all" ]]; then
  echo ""
  echo "========== Per-part summary =========="
  printf '%-8s %12s %8s\n' "Part" "Events" "Files"
  for pdir in "${PART_DIRS[@]}"; do
    pname=$(basename "$pdir")
    # match PartN / partN case-insensitively
    pe="${part_events[$pname]:-0}"
    pf="${part_files[$pname]:-0}"
    # Case-insensitive match (Part4 vs part4)
    pe="${part_events[$pname]:-${part_events[$pname]}}"
    pf="${part_files[$pname]:-${part_files[$pname]}}"
    if [[ "$pe" == "0" && "$pf" == "0" ]]; then
      for k in "${!part_events[@]}"; do
        if [[ "${k,,}" == "${pname,,}" ]]; then
          pe="${part_events[$k]}"
          pf="${part_files[$k]}"
          break
        fi
      done
    fi
    printf '%-8s %12d %8d\n' "$pname" "$pe" "$pf"
  done
  echo "======================================"
fi

echo "Files counted: $files_counted / $files_total"
[[ "$files_failed" -gt 0 ]] && echo "Failed files: $files_failed"
echo "Total events: $total_events"
