#!/usr/bin/env bash
# Run transcriptomics DE/volcano scripts with the /work_space R env (limma installed).
set -euo pipefail
RSCRIPT="/work_space/envs/transcriptomics2/bin/Rscript"
if [[ ! -x "$RSCRIPT" ]]; then
  echo "Missing $RSCRIPT — falling back to system Rscript" >&2
  RSCRIPT="$(command -v Rscript)"
fi
cd "$(dirname "$0")"
exec "$RSCRIPT" "$@"
