dataset_pubchemqc_opv_17458.csv (12 MB) is archived on Zenodo, not mirrored here.

  1. download it from DOI 10.5281/zenodo.18201813
  2. place it in this directory
  3. run  bash ../run_all.sh

Everything else the pipeline needs is in this repository. Two large regenerable
outputs are also Zenodo-only — outputs/recompute/screen_full.csv (4 MB, the
per-molecule results) and genuine_osc_ranked.csv (1.2 MB) — because run_all.sh
rewrites both from the property table. The small result tables, the JSON summaries
and the full run log ARE in the repository, so the headline numbers can be checked
without downloading anything. The PNG copies of the figures are Zenodo-only; the
PDFs here are the versions used in the paper.
