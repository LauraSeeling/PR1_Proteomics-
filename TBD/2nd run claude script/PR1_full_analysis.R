# ==================================================================
# PR1: heart valve disease (PXD023800) - full analysis
#   Part A: transform + filter (valid values per disease x sex cell)
#   Part B: sex-specific detection test (Fisher), catches Y-linked proteins
#   Part C: limma contrasts incl. disease x sex interaction
# Needs: readr, dplyr, tidyr, tibble, ggplot2, limma
# Before running, set the working directory to the folder with the two CSVs
# (basic R GUI: File -> Change dir...)
# ==================================================================

library(readr); library(dplyr); library(tidyr); library(tibble); library(ggplot2)
library(limma)

# ---- Settings ------------------------------------------------------
pg_file        <- "PR1_heart_valve_disease_proteinGroups.csv"
samples_file   <- "PR1_heart_valve_disease_samples.csv"
min_valid_frac <- 0.70   # keep protein if >= 70% quantified in at least ONE
                         # disease x sex cell
min_per_cell   <- 3      # limma only: needs >= this many values in EVERY cell
fdr_cut        <- 0.05
fc_cut         <- 0.58   # |log2FC| ~ 1.5-fold (0 = ignore)

# ==================================================================
# PART A: transform and filter
# ==================================================================

# ---- A1. Read ------------------------------------------------------
pg <- read_csv(pg_file, show_col_types = FALSE,
               col_types = cols(.default = col_guess(),
                                Reverse = col_character(),
                                `Potential contaminant` = col_character(),
                                `Only identified by site` = col_character()))
samples <- read_csv(samples_file, show_col_types = FALSE)

samples <- samples %>%
  mutate(group = recode(group,
                        "control (donor heart)" = "CON",
                        "aortic stenosis"       = "AS",
                        "mitral regurgitation"  = "MR"),
         group = factor(group, levels = c("CON", "AS", "MR")),
         sex   = factor(sex,   levels = c("female", "male")),
         group_sex = interaction(group, sex, sep = "_"))

# ---- A2. Remove flagged rows ---------------------------------------
n0 <- nrow(pg)
pg <- pg %>%
  filter(is.na(Reverse), is.na(`Potential contaminant`),
         is.na(`Only identified by site`))
cat("Protein groups:", n0, "->", nrow(pg), "after removing flagged rows\n")

# ---- A3. Intensity matrix: 0 -> NA, log2 ---------------------------
lfq_cols <- samples$column
stopifnot(all(lfq_cols %in% colnames(pg)))

mat <- pg %>% select(all_of(lfq_cols)) %>% as.matrix()
mat[mat == 0] <- NA
mat <- log2(mat)
colnames(mat) <- samples$sample

# Readable unique row IDs: first gene name, else first UniProt ID
gene <- pg$`Gene names`
gene_first <- sapply(strsplit(ifelse(is.na(gene), "", gene), ";"),
                     function(x) if (length(x) == 0) "" else x[1])
no_gene <- is.na(gene_first) | gene_first == ""
gene_first[no_gene] <- sapply(strsplit(pg$`Protein IDs`[no_gene], ";"), `[`, 1)
rownames(mat) <- make.unique(gene_first)

annot <- pg %>%
  transmute(row_id = rownames(mat),
            protein_ids = `Protein IDs`,
            gene_names  = `Gene names`,
            protein_names = `Protein names`,
            unique_peptides = `Unique peptides`)

cat("Missing values overall:", round(100 * mean(is.na(mat)), 1), "%\n")

# ---- A4. Valid-value filter per disease x sex cell -----------------
grp <- samples$group_sex                       # <-- the fix
valid_frac <- sapply(levels(grp), function(g)
  rowMeans(!is.na(mat[, grp == g, drop = FALSE])))
keep <- rowSums(valid_frac >= min_valid_frac) >= 1

mat_f   <- mat[keep, ]
annot_f <- annot[keep, ]
cat("Proteins kept after valid-value filter:", sum(keep), "of", length(keep), "\n")

# Y-linked proteins: are they kept now?
y_genes <- c("DDX3Y", "EIF1AY", "UTY", "KDM5D", "RPS4Y1", "USP9Y", "ZFY", "NLGN4Y", "TMSB4Y")
cat("Y-linked proteins in filtered data:",
    paste(intersect(y_genes, rownames(mat_f)), collapse = ", "), "\n")

# ---- A5. QC plots --------------------------------------------------
p_box <- as.data.frame(mat_f) %>%
  rownames_to_column("protein") %>%
  pivot_longer(-protein, names_to = "sample", values_to = "log2LFQ") %>%
  left_join(samples, by = "sample") %>%
  ggplot(aes(sample, log2LFQ, fill = group)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6)) +
  labs(title = "log2 LFQ per sample (after filtering)")
print(p_box)
ggsave("QC_boxplot.png", p_box, width = 10, height = 5, dpi = 200)

cc  <- mat_f[complete.cases(mat_f), ]
pca <- prcomp(t(cc), scale. = TRUE)
pve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
p_pca <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample") %>%
  left_join(samples, by = "sample") %>%
  ggplot(aes(PC1, PC2, colour = group, shape = sex)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(x = paste0("PC1 (", pve[1], "%)"), y = paste0("PC2 (", pve[2], "%)"),
       title = paste("PCA on", nrow(cc), "complete proteins"))
print(p_pca)
ggsave("QC_pca.png", p_pca, width = 6, height = 5, dpi = 200)

saveRDS(list(mat = mat_f, annot = annot_f, samples = samples), "PR1_filtered.rds")
write.csv(cbind(annot_f, mat_f), "PR1_filtered_log2LFQ.csv", row.names = FALSE)

# ==================================================================
# PART B: detection test by sex (for proteins limma cannot handle)
# ==================================================================
sex <- samples$sex
det <- data.frame(
  row_id     = rownames(mat_f),
  f_detected = rowSums(!is.na(mat_f[, sex == "female"])),
  m_detected = rowSums(!is.na(mat_f[, sex == "male"])),
  n_f = sum(sex == "female"),
  n_m = sum(sex == "male")
)
det$p <- mapply(function(a, b, nf, nm)
                  fisher.test(matrix(c(a, nf - a, b, nm - b), 2))$p.value,
                det$f_detected, det$m_detected, det$n_f, det$n_m)
det$adj.p <- p.adjust(det$p, "BH")
det <- det %>% left_join(annot_f, by = "row_id") %>% arrange(p)
write.csv(det, "detection_by_sex.csv", row.names = FALSE)

cat("\nProteins with different detection rate between sexes (FDR < 0.05):",
    sum(det$adj.p < fdr_cut), "\n")
print(head(det[, c("row_id", "f_detected", "m_detected", "adj.p")], 10))

# ==================================================================
# PART C: limma (disease x sex)
# ==================================================================

# limma needs every cell to have data; proteins limited to some cells
# are covered by Part B instead
cell_counts <- sapply(levels(samples$group_sex), function(g)
  rowSums(!is.na(mat_f[, samples$group_sex == g, drop = FALSE])))
testable <- rowSums(cell_counts >= min_per_cell) == ncol(cell_counts)
mat_l   <- mat_f[testable, ]
annot_l <- annot_f[testable, ]
cat("\nProteins testable with limma (>=", min_per_cell, "values in every cell):",
    sum(testable), "of", nrow(mat_f), "\n")

# ---- C1. Design ----------------------------------------------------
design <- model.matrix(~ 0 + group_sex, data = samples)
colnames(design) <- sub("group_sex", "", colnames(design))
print(colSums(design))

# ---- C2. Contrasts -------------------------------------------------
contrasts <- makeContrasts(
  AS_vs_CON_F = AS_female - CON_female,
  AS_vs_CON_M = AS_male   - CON_male,
  MR_vs_CON_F = MR_female - CON_female,
  MR_vs_CON_M = MR_male   - CON_male,
  AS_vs_CON   = (AS_female + AS_male)/2 - (CON_female + CON_male)/2,
  MR_vs_CON   = (MR_female + MR_male)/2 - (CON_female + CON_male)/2,
  AS_vs_MR    = (AS_female + AS_male)/2 - (MR_female + MR_male)/2,
  M_vs_F      = (CON_male + AS_male + MR_male)/3 -
                (CON_female + AS_female + MR_female)/3,
  AS_x_sex    = (AS_male - CON_male) - (AS_female - CON_female),
  MR_x_sex    = (MR_male - CON_male) - (MR_female - CON_female),
  levels = design
)

# ---- C3. Fit -------------------------------------------------------
fit <- lmFit(mat_l, design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)

# ---- C4. Results ---------------------------------------------------
results <- lapply(colnames(contrasts), function(cn) {
  topTable(fit, coef = cn, number = Inf, adjust.method = "BH", sort.by = "P") %>%
    rownames_to_column("row_id") %>%
    left_join(annot_l, by = "row_id") %>%
    mutate(contrast = cn,
           significant = adj.P.Val < fdr_cut & abs(logFC) >= fc_cut)
})
names(results) <- colnames(contrasts)

for (cn in names(results))
  write.csv(results[[cn]], paste0("limma_", cn, ".csv"), row.names = FALSE)

summary_tab <- bind_rows(results) %>%
  group_by(contrast) %>%
  summarise(tested = n(),
            sig_FDR05 = sum(adj.P.Val < fdr_cut, na.rm = TRUE),
            sig_FDR05_FC = sum(significant, na.rm = TRUE),
            up = sum(significant & logFC > 0, na.rm = TRUE),
            down = sum(significant & logFC < 0, na.rm = TRUE))
print(as.data.frame(summary_tab))
write.csv(summary_tab, "limma_summary.csv", row.names = FALSE)

# ---- C5. Volcano plots ---------------------------------------------
volcano <- function(res, title, n_label = 10) {
  res <- res %>%
    mutate(status = case_when(significant & logFC > 0 ~ "up",
                              significant & logFC < 0 ~ "down",
                              TRUE ~ "ns"))
  lab <- res %>% filter(significant) %>%
    slice_min(adj.P.Val, n = n_label, with_ties = FALSE)
  ggplot(res, aes(logFC, -log10(P.Value), colour = status)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text(data = lab, aes(label = row_id), size = 2.5, vjust = -0.5,
              check_overlap = TRUE, show.legend = FALSE) +
    scale_colour_manual(values = c(up = "firebrick", down = "steelblue", ns = "grey70")) +
    theme_bw() +
    labs(title = title, x = "log2 fold change", y = "-log10(p-value)")
}
for (cn in names(results)) {
  p <- volcano(results[[cn]], cn)
  print(p)
  ggsave(paste0("volcano_", cn, ".png"), p, width = 6, height = 5, dpi = 200)
}

# ---- C6. Sanity checks ---------------------------------------------
print(head(results$M_vs_F[, c("row_id", "logFC", "adj.P.Val")], 10))
hist(results$AS_x_sex$P.Value, breaks = 30,
     main = "AS x sex: p-value distribution", xlab = "p-value")

saveRDS(list(fit = fit, results = results, design = design,
             contrasts = contrasts, detection = det),
        "PR1_limma_results.rds")
cat("\nDone. Results written to the working directory.\n")
