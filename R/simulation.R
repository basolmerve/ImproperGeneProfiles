# Simulation study
library(tidyverse)
library(parallel)
library(DESeq2)
library(nsROC)
library(patchwork)
library(ggrepel)
library(ggtext)
source("R/helper_functions.R")

# 1 - Simulations ----
## 1.1 Scenarios ----
# We will run 27 combinations in total, each repeated 1000 times.
simulation_scenarios <- expand_grid(
  simID = NA_integer_,
  n = c(100, 300, 500),
  p = 3000,
  propDE = c(.05, .30),
  propIR = .25,
  phi = c(.01, .1, 1),
  sdsignal = 1.5,
  propNZ = .10
) |> 
  mutate(
    simID = 1:n()
  )

# Number of replicated simulations for each scenario. (e.g., 250, 1000, etc.)
nSim <- 10

# Register parallel backend via "parallel" package.
cl <- makeCluster(10)   # Number of CPUs activated.

# Export required elements to parallel nodes.
clusterExport(
  cl, 
  c(
    "createDDSobject", "filterCounts", "diffExp", "selectDEfeatures_DESeq",
    "preProcessCounts", "diffExp_ROC", "selectDEfeatures_ROC", "nLhat",
    "generateCountData"
  )
)

# Select specific scenarios for testing purpose. Uncomment as needed.
# simulation_scenarios <- simulation_scenarios[1, ]
# simulation_scenarios <- simulation_scenarios[15:18, ]
simRes <- list()

## 1.2 Computations ----
for (i in 1:nrow(simulation_scenarios)){
  scen.i <- simulation_scenarios[i, ]
  nFeat_select <- ceiling(scen.i$propDE * scen.i$p)
  
  printStatus(idx = i, scen.i, nrow(simulation_scenarios))
  
  # Perform analysis in parallel cores for data list (DF_List)
  clusterExport(cl, c("nFeat_select", "scen.i"))
  
  # Set seed within clusters
  clusterSetRNGStream(cl, iseed = 3627)
  simRes.i <- parLapply(cl, X = 1:nSim, function(idx = X, nFeat = nFeat_select, scen = scen.i, ...){
    library(tidyverse)
    library(DESeq2)
    
    # Step 1: Generate RNA-Seq Data and store in a list
    DF <- generateCountData(
      n = scen$n, p = scen$p, K = 2, 
      param = 1 / scen$phi, sdsignal = scen$sdsignal, 
      DE = scen$propDE, IR = scen$propIR, 
      nonzero_prop = 0, min_count = 1,
      allZero.rm = TRUE, tag.samples = TRUE
    )
    
    # Create DESeqDataSet object.
    dds <- createDDSobject(DF)
    
    # Step 2: Filter data (Remove low quality genes with near zero variances.
    dds_processed <- filterCounts(dds)
    
    # We will force each method to select all DE features (one-shot). nFeat is the number of DEGs available in 
    # the generated RNA-Sequence dataset.
    nFeat <- pmin(nFeat, length(dds_processed$DE_Genes))
    
    # DE Analysis via DESeq2
    dds_diffExp <- diffExp(dds_processed, nonzero = TRUE)
    selectedGenes_DESeq <- selectDEfeatures_DESeq(dds_diffExp, nFeatures = nFeat)
    
    # Step 3: Transform filtered counts using VST method.
    dds_processed <- preProcessCounts(dds_processed, normalize = TRUE, transform = TRUE, 
                                      transformationMethod = "vst", nonzero = TRUE)
    
    # DE Analysis via ROC-based methods.
    dds_diffExp_ROC <- diffExp_ROC(.object = dds_processed)
    selectedGenes_ROC <- selectDEfeatures_ROC(dds_diffExp_ROC, nFeat = nFeat)
    
    selectedGenes_all <- bind_cols(
      selectedGenes_ROC$selectedFeatures, 
      selectedGenes_DESeq$selectedFeatures
    )
    
    simParams <- bind_cols(scen, tibble(simReplication = idx))
    
    # True Detection Rates for Improper Genes
    TDR_IGs <- bind_cols(
      simParams,
      tibble(
        group = "IGs",
        LROC = sum(selectedGenes_all$LROC %in% dds_diffExp_ROC$improper_Genes) / length(dds_diffExp_ROC$improper_Genes),
        gAUC = sum(selectedGenes_all$gAUC %in% dds_diffExp_ROC$improper_Genes) / length(dds_diffExp_ROC$improper_Genes),
        AUC = sum(selectedGenes_all$AUC %in% dds_diffExp_ROC$improper_Genes) / length(dds_diffExp_ROC$improper_Genes),
        DESeq = sum(selectedGenes_all$DESeq %in% dds_diffExp_ROC$improper_Genes) / length(dds_diffExp_ROC$improper_Genes)
      )
    ) 
    
    # Overall True Detection Rates for DE Genes (including proper and improper profiles)
    TDR_DEGs <- bind_cols(
      simParams,
      tibble(
        group = "DEGs",
        LROC = sum(selectedGenes_all$LROC %in% dds_diffExp_ROC$DE_Genes) / length(dds_diffExp_ROC$DE_Genes),
        gAUC = sum(selectedGenes_all$gAUC %in% dds_diffExp_ROC$DE_Genes) / length(dds_diffExp_ROC$DE_Genes),
        AUC = sum(selectedGenes_all$AUC %in% dds_diffExp_ROC$DE_Genes) / length(dds_diffExp_ROC$DE_Genes),
        DESeq = sum(selectedGenes_all$DESeq %in% dds_diffExp_ROC$DE_Genes) / length(dds_diffExp_ROC$DE_Genes)
      )
    )
    
    TDR <- bind_rows(TDR_IGs, TDR_DEGs)
    return(TDR)
  })

  simRes[[i]] <- bind_rows(simRes.i)

  # Uncomment the line below to save intermediate results.
  # save(simRes, file = "saved/Simulation/simRes.Rda")
}

# Save final results of simulations. Uncomment if needed.
# dir.create("saved/Simulation", recursive = TRUE, showWarnings = FALSE)
# save(simRes, file = "saved/Simulation/simRes.Rda")

# 2. Plots ----
# 2.1 Performances ----
load("saved/Simulation/simRes.Rda")
simRes <- bind_rows(simRes)

# Convert to long format for ggplot visualizations.
simRes_long <- simRes |> 
  select(-c(p, propIR)) |> 
  pivot_longer(cols = LROC:DESeq, names_to = "Method", values_to = "Value") |> 
  mutate(
    Method = factor(Method, levels = c("DESeq", "AUC", "gAUC", "LROC")),
    group = factor(group),
    phi = factor(phi, levels = c(.01, .1, 1), labels = c("Very Low", "Moderate", "Very High")),
    n = factor(n, levels = c(100, 300, 500), labels = c("100", "300", "500")),
    propDE = factor(propDE, levels = c(.05, .30), labels = c("Low", "High"))
  )

# Overall performances for low DE scenarios.
fig_low <- ggplot(filter(simRes_long, propDE == "Low"), aes(x = Method, y = Value, fill = n)) + 
  geom_boxplot(width = .6, outlier.colour = rgb(0, 0, 0, .3)) + 
  theme_bw() + 
  facet_grid(cols = vars(group), rows = vars(phi), scales = "free_y") + 
  theme(legend.position = "top",
        axis.text.x = element_text(margin = margin(t = 5, b = 5)),
        axis.text.y = element_text(margin = margin(r = 5, l = 5))) + 
  guides(fill = guide_legend(title = "Sample size (n) ")) + 
  labs(x = "Differential Expression Methods", y = "True Detection Rate")

# ggsave(filename = "figure/Simulation/simRes_overall_low.png", plot = fig_low, device = png, 
#        width = 18, height = 22, units = "cm", dpi = 320, create.dir = TRUE)
# ggsave(filename = "document/manuscript/figure/simRes_overall_low.png", plot = fig_low, device = png, 
#        width = 18, height = 22, units = "cm", dpi = 320, create.dir = TRUE)

# Overall performances for high DE scenarios.
fig_high <- ggplot(filter(simRes_long, propDE != "Low"), aes(x = Method, y = Value, fill = n)) + 
  geom_boxplot(width = .6, outlier.colour = rgb(0, 0, 0, .3)) + 
  theme_bw() + 
  facet_grid(cols = vars(group), rows = vars(phi), scales = "free_y") + 
  theme(legend.position = "top",
        axis.text.x = element_text(margin = margin(t = 5, b = 5)),
        axis.text.y = element_text(margin = margin(r = 5, l = 5))) + 
  guides(fill = guide_legend(title = "Sample size (n) ")) + 
  labs(x = "Differential Expression Methods", y = "True Detection Rate")

# Save figure
# ggsave(filename = "figure/Simulation/simRes_overall_high.png", plot = fig_high, device = png, 
#        width = 18, height = 22, units = "cm", dpi = 320, create.dir = TRUE)

# ggsave(filename = "document/manuscript/figure/simRes_overall_high.png", plot = fig_high, device = png, 
#        width = 18, height = 22, units = "cm", dpi = 320, create.dir = TRUE)

# ggsave(filename = "figure/Simulation/simRes_overall.svg", plot = fig, device = svg, 
#        width = 18, height = 26, units = "cm", dpi = 320, create.dir = TRUE)

# 2.2 IG Inspection ----
# Here, we focused on a specific scenario for illustration.
set.seed(3627)
dat <- generateCountData(
  n = 300, p = 1000, K = 2, param = 1, sdsignal = 1.5, 
  tag.samples = TRUE, allZero.rm = TRUE, nonzero_prop = .10
)

dds <- createDDSobject(dat)

# Step 2: Filter data (Remove low quality genes with near zero variances.)
dds_processed <- filterCounts(dds)

# 'nFeat' is the number of DE features selected from complete dataset.
nFeat <- 300

# DE Analysis via DESeq2. Select DE features via DESeq2. Number of features is 'nFeat'.
dds_diffExp <- diffExp(dds_processed)
selectedGenes_DESeq <- selectDEfeatures_DESeq(dds_diffExp, nFeatures = nFeat)

# Step 3: Transform filtered counts using VST method available in DESeq2. 
# VST transformed counts will be used for ROC-based DE analysis.
dds_processed <- preProcessCounts(dds_processed, normalize = TRUE, transform = TRUE, transformationMethod = "vst")

# Extract sample information from DESeqDataSet object.
col_data <- colData(dds_processed$DESeqObject) |> 
  as.data.frame() |> 
  mutate(
    condition01 = if_else(condition == "C2", 1, 0),
    condition = factor(condition, levels = c("C1", "C2"))
  ) |> 
  DataFrame()

# Prepare data for ROC analysis.
rocData <- bind_cols(
  as_tibble(t(dds_processed$transformedCounts)), 
  tibble(response = col_data$condition01)
)

# DE Analysis via ROC-based methods.
# For computational efficiency, we perform parallel processing here. 
# If 'cluster' is NULL or not provided, computations will be done sequentially.
dds_diffExp_ROC <- diffExp_ROC(.object = dds_processed, cluster = cl)

# Select DE features via ROC-based methods. Number of features is 'nFeat'.
selectedGenes_ROC <- selectDEfeatures_ROC(dds_diffExp_ROC, nFeat = nFeat)

# Combine selected genes from both methods (i.e., DESeq2 and ROC-based).
selectedGenes <- bind_cols(
  selectedGenes_ROC$selectedFeatures, 
  selectedGenes_DESeq$selectedFeatures
)

# Results for each method with the values of metrics. 
# Metrics are AUC values for ROC-based methods and log-fold changes for DESeq2.
results <- bind_rows(
  selectedGenes_DESeq$results,
  selectedGenes_ROC$results
)

# Combine all results in a list.
results_DiffExp <- list(
  results = results,
  selectedFeatures = selectedGenes,
  rocData = rocData
)

# True Detection Rates for Improper Genes
improperGenes <- dds_diffExp_ROC$improper_Genes
DEGs <- dds_diffExp_ROC$DE_Genes

TDR_IGs <- tibble(
  LROC = sum(selectedGenes$LROC %in% improperGenes) / length(improperGenes),
  gAUC = sum(selectedGenes$gAUC %in% improperGenes) / length(improperGenes),
  AUC = sum(selectedGenes$AUC %in% improperGenes) / length(improperGenes),
  DESeq = sum(selectedGenes$DESeq %in% improperGenes) / length(improperGenes)
)

DE_Results <- results_DiffExp$results |> 
  arrange(Method) |>
  pivot_wider(id_cols = Gene, names_from = Method, values_from = Value) |>
  mutate(
    improper = factor(if_else(Gene %in% improperGenes, 1, 0)),
    DEGs = factor(if_else(Gene %in% DEGs, 1, 0)),
    type = factor(
      case_when(
        (DEGs == 1 & improper == 1) ~ "IG",
        (DEGs == 1 & improper == 0) ~ "DEG",
        .default = "Not significant"
      )
    )
  )

axisTheme <- theme(
  axis.text.x = element_text(margin = margin(t = 5, b = 5)),
  axis.text.y = element_text(margin = margin(r = 5, l = 5))
)
baseFontSize <- 12

LROC_AUC <- ggplot(DE_Results, aes(x = AUC, y = LROC, colour = type)) +
  geom_point(size = 2) + 
  theme_bw(base_size = baseFontSize) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_markdown(lineheight = 3)
  ) + 
  axisTheme + 
  scale_colour_manual(values = c(rgb(1, 0, 0, .2), rgb(0, 0.3, 1, .3), "#0000004D"), labels = c("Diff. Exp. Genes", "Improper Genes", "Not significant")) + 
  guides(colour = "none") + 
  labs(
    x = "Area Under Curve (cAUC)<br><span style='display:block; text-align:center; font-size:18pt;'>(a)</span>",
    y = "Length of ROC Curve (LROC)"
  )
# guides(colour = guide_legend(override.aes = list(colour = c(rgb(1, 0, 0), rgb(0, 0.3, 1), "#000000")), position = "top", title = NULL))

gAUC_AUC <- ggplot(DE_Results, aes(x = AUC, y = gAUC, colour = type)) +
  geom_point(size = 2) + 
  geom_abline(slope = 1, intercept = 0, colour = "gray30", linetype = 2) + 
  geom_hline(yintercept = .5, linetype = 2, colour = "gray60") + 
  geom_vline(xintercept = .5, linetype = 2, colour = "gray60") + 
  geom_label(inherit.aes = FALSE, x = .65, y = .45, label = "gAUC < 0.5") + 
  theme_bw(base_size = baseFontSize) +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_text(margin = margin(l = 10)),
    axis.title.x = element_markdown(lineheight = 3)
  ) + 
  axisTheme + 
  scale_colour_manual(values = c(rgb(1, 0, 0, .2), rgb(0, 0.3, 1, .3), "#0000004D"), labels = c("Diff. Exp. Genes", "Improper Genes", "Not significant")) + 
  guides(colour = guide_legend(override.aes = list(colour = c(rgb(1, 0, 0), rgb(0, 0.3, 1), "#000000")), position = "top", title = NULL)) + 
  labs(
    x = "Area Under Curve (cAUC)<br><span style='display:block; text-align:center; font-size:18pt;'>(b)</span>",
    y = "Generalized Area Under Curve (gAUC)"
  )

gAUC_LROC <- ggplot(DE_Results, aes(x = LROC, y = gAUC, colour = type)) +
  geom_point(size = 2) + 
  theme_bw(base_size = baseFontSize) +
  theme(
    panel.grid = element_blank(),
    axis.title.y = element_text(margin = margin(l = 10)),
    axis.title.x = element_markdown(lineheight = 3)
  ) + 
  axisTheme + 
  scale_colour_manual(values = c(rgb(1, 0, 0, .2), rgb(0, 0.3, 1, .3), "#0000004D"), labels = c("Diff. Exp. Genes", "Improper Genes", "Not significant")) + 
  guides(colour = "none") +
  labs(
    y = "Generalized Area Under Curve (gAUC)",
    # Başlığı iki satır yapıyoruz; (a) için ayrı stil verebilirsin
    x = "Length of ROC Curve (LROC)<br><span style='display:block; text-align:center; font-size:18pt'>(c)</span>"
  )
# guides(colour = guide_legend(override.aes = list(colour = c(rgb(1, 0, 0), rgb(0, 0.3, 1), "#000000")), position = "top", title = NULL))
gAUC_LROC

smoothData <- DE_Results |> 
  filter(type != "Not significant")

gAUC_LROC <- gAUC_LROC + 
  geom_smooth(data = smoothData, mapping = aes(x = LROC, y = gAUC, colour = type), se = FALSE) + 
  geom_hline(yintercept = .5, linetype = 2, colour = "gray60")

figALL <- (LROC_AUC + gAUC_AUC + gAUC_LROC) + 
  plot_layout() &
  theme(plot.margin = margin(t = 0, b = 0, r = 5, l = 5))

# Save the figure
# ggsave(filename = "figure/Simulation/IG_visual_inspect.png", plot = figALL, device = png, 
#        width = 28, height = 12, units = "cm", dpi = 320)

# ggsave(filename = "document/manuscript/figure/IG_visual_inspect.png", plot = figALL, device = png, 
#        width = 28, height = 12, units = "cm", dpi = 320)