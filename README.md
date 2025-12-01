# ImproperGeneProfiles

Code and data accompanying the manuscript “Enhancing Differential Expression Analysis by ROC-Based Detection of Improper Expression Profiles.”

## What this repository does
- Implements ROC-based indices (generalized AUC, gAUC, and Length of ROC, LROC) to flag non-monotonic “improper” gene/miRNA expression profiles that standard DE pipelines often miss.
- Benchmarks gAUC and LROC against classical AUC (cAUC) and DESeq2 in extensive negative-binomial simulations spanning sample size (100/300/500), dispersion (very low/moderate/very high), and differential expression prevalence (5% vs 30%).
- Applies the workflow to a cervical cancer miRNA-Seq dataset (29 tumors vs 29 normals), showing that gAUC robustly recovers improper miRNAs and LROC provides complementary signal at low–moderate dispersion.
- Produces figures and saved R objects for both simulation and real-data analyses for downstream visualization and inspection.

## Repository layout
- `R/`
  - `simulation.R`: Full simulation study comparing DESeq2, cAUC, gAUC, and LROC; generates detection-rate summaries and figures.
  - `cervical.R`: Re-analysis of the cervical cancer miRNA-Seq dataset; runs DESeq2 and ROC-based screens and saves results/plots.
  - `helper_functions.R`: Utilities for data generation, preprocessing (filtering, VST), ROC-based calculations, and plotting helpers.
- `data/`
  - `count.csv`: Raw miRNA count matrix (genes x samples) for the cervical cancer study.
  - `condition.csv`: Sample metadata with class labels used by `cervical.R`.
- `figure/`
  - `Simulation/`: PNGs of simulation performance summaries and improper-gene inspections.
  - `Cervical/`: Plots for the cervical cohort (volcano, joint ROC-based diagnostics).
  - `volcano_plot.pdf`: Additional DE visualization.
- `saved/`
  - `Simulation/simRes.Rda`: Serialized simulation results (true detection rates across scenarios).
  - `Cervical/Cervical_results.Rda`: Serialized DE/ROC results for the cervical dataset.
- `LROC.Rproj`: RStudio project file for convenience.
- `LICENSE`: Repository license.

## Methodological summary (from the manuscript)
- **Problem**: Standard count-based DE tools (e.g., DESeq2) are tuned to monotonic shifts and can miss “Goldilocks” profiles where both very low and very high expression are disease-associated.
- **Approach**: Extend ROC analysis beyond classical AUC.
  - **gAUC** generalizes AUC by allowing dual thresholds to capture tail-driven signals.
  - **LROC** scores diagnostic ability via ROC curve length, highlighting curvature that cAUC underestimates.
  - **cAUC** and **DESeq2** serve as baselines.
- **Simulations**: Negative-binomial counts with improper profiles induced via high–low mixtures; 18 scenarios (sample sizes 100/300/500; dispersion 0.01/0.1/1; DE proportion 0.05/0.30; 25% of DEGs improper; 3,000 features; 1,000 replicates per scenario in the manuscript, `nSim` set smaller in the shared script for quick runs).
  - **Key finding**: gAUC delivers the highest and most stable detection of improper features; LROC adds complementary power at low–moderate dispersion but degrades under extreme overdispersion; cAUC performs worst for improper profiles.
- **Real data (cervical miRNA-Seq, 29 vs 29)**: gAUC and LROC highlight additional biologically plausible improper miRNAs beyond DESeq2 hits; joint gAUC–LROC–cAUC plots visually separate improper from proper profiles.

## How to run
1. Open the project in R or RStudio (`LROC.Rproj`).
2. Ensure required packages: `DESeq2`, `tidyverse`, `nsROC`, `parallel`, `patchwork`, `ggrepel`, `ggtext` (see `R/simulation.R` and `R/cervical.R` for full list).
3. Simulations: run `Rscript R/simulation.R` to regenerate detection-rate summaries and figures (adjust `nSim` and cluster settings as needed).
4. Cervical analysis: run `Rscript R/cervical.R` to reproduce the miRNA study outputs; results saved to `saved/Cervical/Cervical_results.Rda` and plots in `figure/Cervical/`.

## Outputs to expect
- Saved R objects in `saved/` capturing DESeq2 and ROC-based statistics for downstream exploration.
- Figures in `figure/` summarizing simulation performance, improper-profile inspections, and cervical miRNA volcano/ROC diagnostics.

## Notes
- Scripts assume working directory is the repository root.
- The shared scripts load pre-saved results for plotting; uncomment the `save`/`ggsave` lines inside the scripts to regenerate from scratch.
