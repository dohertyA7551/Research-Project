#!/usr/bin/env bash
# Full TMA composition pipeline: assign groups → proteomics → transcriptomics.
set -euo pipefail
R=/work_space/envs/transcriptomics2/bin/Rscript
echo "=== Step 1: TMA composition group assignment ==="
cd /work_space/files/proteomics
"$R" run_tma_composition_groups.R
echo "=== Step 2: Proteomics ==="
"$R" run_tma_comp_proteomics.R
echo "=== Step 3: Transcriptomics ==="
cd /work_space/files/Transcriptomics
"$R" run_tma_comp_transcriptomics.R
echo "=== Done ==="
