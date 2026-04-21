#!/usr/bin/env bash

set -euo pipefail

# Detect schedd nodes where current user has queued/running/held jobs.
# Usage:
#   ./detect_active_schedds.sh
#   ./detect_active_schedds.sh bigbird13.cern.ch bigbird19.cern.ch
#
# Environment variables:
#   NODE_PREFIX   default: bigbird
#   NODE_START    default: 1
#   NODE_END      default: 40
#   NODE_SUFFIX   default: .cern.ch
#   QUERY_TIMEOUT default: 8 (seconds per node)

USER_NAME="${USER:-$(id -un)}"
NODE_PREFIX="${NODE_PREFIX:-bigbird}"
NODE_START="${NODE_START:-1}"
NODE_END="${NODE_END:-40}"
NODE_SUFFIX="${NODE_SUFFIX:-.cern.ch}"
QUERY_TIMEOUT="${QUERY_TIMEOUT:-8}"

build_default_nodes() {
  local i
  for ((i = NODE_START; i <= NODE_END; i++)); do
    printf "%s%s%s\n" "${NODE_PREFIX}" "${i}" "${NODE_SUFFIX}"
  done
}

get_user_jobs_on_node() {
  local node="$1"
  local out

  # Query this schedd directly; avoid changing global myschedd assignment.
  if ! out="$(timeout "${QUERY_TIMEOUT}" condor_q -name "${node}" "${USER_NAME}" 2>/dev/null)"; then
    return 1
  fi

  # For `condor_q -name <schedd> <user>`, "Total for query" equals this user's jobs.
  # Some environments omit "Total for <user>", so parse query first, then user as fallback.
  awk -v user="${USER_NAME}" '
    $1 == "Total" && $2 == "for" && $3 == "query:" {
      gsub(/[^0-9]/, "", $4)
      print $4
      found=1
      exit
    }
    $1 == "Total" && $2 == "for" && $3 == user ":" {
      gsub(/[^0-9]/, "", $4)
      print $4
      found=1
      exit
    }
    END {
      if (!found) print 0
    }
  ' <<< "${out}"
}

print_node_rows() {
  local node="$1"
  local out

  if ! out="$(timeout "${QUERY_TIMEOUT}" condor_q -name "${node}" "${USER_NAME}" 2>/dev/null)"; then
    return 1
  fi

  awk -v user="${USER_NAME}" -v node="${node}" '
    BEGIN {
      printed=0
      owner_re = "^" user "[[:space:]]"
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
  ' <<< "${out}"
}

main() {
  local -a nodes
  local node jobs

  if (($# > 0)); then
    nodes=("$@")
  else
    mapfile -t nodes < <(build_default_nodes)
  fi

  local found_any=0
  for node in "${nodes[@]}"; do
    jobs="$(get_user_jobs_on_node "${node}" || true)"
    if [[ -n "${jobs}" ]] && [[ "${jobs}" =~ ^[0-9]+$ ]] && ((jobs > 0)); then
      print_node_rows "${node}" || true
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
