#!/usr/bin/env bash

set -euo pipefail

# Detect schedd nodes where current user has queued/running/held jobs.
# Usage:
#   ./detect_active_schedds.sh
#   ./detect_active_schedds.sh bigbird13.cern.ch bigbird19.cern.ch

USER_NAME="${USER:-$(id -un)}"
NODE_PREFIX="${NODE_PREFIX:-bigbird}"
NODE_START="${NODE_START:-1}"
NODE_END="${NODE_END:-101}"
NODE_SUFFIX="${NODE_SUFFIX:-.cern.ch}"
QUERY_TIMEOUT="${QUERY_TIMEOUT:-8}"
MAX_PARALLEL="${MAX_PARALLEL:-12}"

build_default_nodes() {
  local i
  for ((i = NODE_START; i <= NODE_END; i++)); do
    printf "%s%s%s\n" "${NODE_PREFIX}" "${i}" "${NODE_SUFFIX}"
  done
}

query_node_once() {
  local node="$1"
  local out

  # Query this schedd directly; avoid changing global myschedd assignment.
  if ! out="$(timeout "${QUERY_TIMEOUT}" condor_q -name "${node}" "${USER_NAME}" 2>/dev/null)"; then
    return 0
  fi

  awk -v user="${USER_NAME}" -v node="${node}" '
    BEGIN {
      printed=0
      jobs=0
      owner_re = "^" user "[[:space:]]"
    }
    # For `condor_q -name <schedd> <user>`, "Total for query" equals this user jobs.
    # Some environments omit "Total for <user>", so parse query first, then user as fallback.
    $1 == "Total" && $2 == "for" && $3 == "query:" {
      val=$4
      gsub(/[^0-9]/, "", val)
      if (val != "") jobs = val + 0
      next
    }
    $1 == "Total" && $2 == "for" && $3 == user ":" {
      val=$4
      gsub(/[^0-9]/, "", val)
      if (val != "") jobs = val + 0
      next
    }
    /^OWNER[[:space:]]+/ {
      header=$0
      next
    }
    $0 ~ owner_re {
      if (!printed) {
        print "[" node "]"
        if (header != "") print header
        printed=1
      }
      print
    }
    END {
      if (!printed && jobs > 0) {
        print "[" node "]"
        if (header != "") print header
      }
    }
  ' <<< "${out}"
}

run_one_node() {
  local idx="$1"
  local node="$2"
  local out_dir="$3"
  local out_file="${out_dir}/${idx}.out"

  query_node_once "${node}" > "${out_file}" || true
}

main() {
  local -a nodes
  local node
  local out_dir
  local -a pids
  local i pid
  local running=0
  local found_any=0

  if (($# > 0)); then
    nodes=("$@")
  else
    mapfile -t nodes < <(build_default_nodes)
  fi

  out_dir="$(mktemp -d)"
  trap '[[ -n "${out_dir:-}" ]] && rm -rf "${out_dir}"' EXIT

  for i in "${!nodes[@]}"; do
    node="${nodes[$i]}"
    run_one_node "$i" "${node}" "${out_dir}" &
    pid=$!
    pids+=("${pid}")
    ((running += 1))

    if ((running >= MAX_PARALLEL)); then
      wait -n || true
      ((running -= 1))
    fi
  done

  for pid in "${pids[@]}"; do
    wait "${pid}" || true
  done

  for i in "${!nodes[@]}"; do
    if [[ -s "${out_dir}/${i}.out" ]]; then
      while IFS= read -r line; do
        [[ -z "${line}" ]] && continue
        echo "${line}"
      done < "${out_dir}/${i}.out"
      echo
      found_any=1
    fi
  done

  if ((found_any == 0)); then
    echo "No active jobs found on scanned nodes."
    exit 0
  fi
}

main "$@"
