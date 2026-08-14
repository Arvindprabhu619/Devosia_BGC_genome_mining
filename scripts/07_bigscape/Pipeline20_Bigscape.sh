#!/usr/bin/env bash

set -euo pipefail

# ============================================================
# BiG-SCAPE clustering
#
# Repository-relative, portable version
#
# Required environment variable:
#   PFAM_DB
#
# Example:
#   export PFAM_DB=/path/to/Pfam-A.hmm
#
# ============================================================

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

ANTISMASH_DIR="${ANTISMASH_DIR:-${REPO_ROOT}/antismash_output}"
OUTPUT_DIR="${OUTPUT_DIR:-${REPO_ROOT}/bigscape_output_full}"
PFAM_DB="${PFAM_DB:-}"

if [[ -z "${PFAM_DB}" ]]; then
    echo "ERROR: PFAM_DB is not set."
    echo
    echo "Set it before running BiG-SCAPE, for example:"
    echo "  export PFAM_DB=/path/to/Pfam-A.hmm"
    exit 1
fi

if [[ ! -f "${PFAM_DB}" ]]; then
    echo "ERROR: Pfam database not found:"
    echo "  ${PFAM_DB}"
    exit 1
fi

if [[ ! -d "${ANTISMASH_DIR}" ]]; then
    echo "ERROR: antiSMASH directory not found:"
    echo "  ${ANTISMASH_DIR}"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

echo "============================================================"
echo "BiG-SCAPE"
echo "============================================================"
echo "Repository : ${REPO_ROOT}"
echo "Input      : ${ANTISMASH_DIR}"
echo "Output     : ${OUTPUT_DIR}"
echo "Pfam       : ${PFAM_DB}"
echo "============================================================"

bigscape cluster \
    -i "${ANTISMASH_DIR}" \
    -o "${OUTPUT_DIR}" \
    --input-mode recursive \
    -m 3.1 \
    -p "${PFAM_DB}" \
    --gcf-cutoffs 0.3 \
    --classify category \
    --record-type region \
    --include-singletons \
    -c 8 \
    -l Devosia_110_full

echo
echo "BiG-SCAPE completed successfully."
