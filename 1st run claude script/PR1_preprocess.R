# ------------------------------------------------------------------
# PR1: heart valve disease (PXD023800) - transform and filter
# Input : proteinGroups CSV + samples CSV
# Output: log2 matrix (NA = not quantified) filtered on valid values,
#         plus sample table, ready for limma
# ------------------------------------------------------------------

library(tidyverse)

# ---- 0. Settings --------------------------------------------------
pg_file      <- "PR1_heart_valve_disease_proteinGroups.csv"
samples_file <- "PR1_heart_valve_disease_samples.csv"
min_valid_frac <- 0.70   # protein must be quantified in >= 70% of samples
                         # of at least one group (disease x sex groups
                         # can be used instead, see step 5)

# ---- 1. Read data -------------------------------------------------
pg      <- read_csv(pg_file, show_col_types = FALSE,
                    col_types = cols(.default = col_guess(),
                                     Reverse = col_character(),
                                     `Potential contaminant` = col_character(),
                                     `Only identified by site` = col_character()))
samples <- read_csv(samples_file, show_col_types = FALSE)

# Short group labels and factors (control as reference level)
samples <- samples %>%
  mutate(group = recode(group,
                        "control (donor heart)" = "CON",
                        "aortic stenosis"       = "AS",
                        "mitral regurgitation"  = "MR"),
         group = factor(group, levels = c("CON", "AS", "MR")),
         sex   = factor(sex,   levels = c("female", "male")),
         group_sex = interaction(group, sex, sep = "_"))

# ---- 2. Remove reverse, contaminants, only-identified-by-site -----
n0 <- nrow(pg)
pg <- pg %>%
  filter(is.na(Reverse), is.na(`Potential contaminant`),
         is.na(`Only identified by site`))
cat("Protein groups:", n0, "->", nrow(pg), "after removing flagged rows\n")

# ---- 3. Build intensity matrix; 0 -> NA; log2 ---------------------
lfq_cols <- samples$column                       # same order as samples table
stopifnot(all(lfq_cols %in% colnames(pg)))

mat <- pg %>% select(all_of(lfq_cols)) %>% as.matrix()
mat[mat == 0] <- NA                              # 0 = not quantified
mat <- log2(mat)
colnames(mat) <- samples$sample

# Unique, readable row IDs (gene names can be missing or duplicated)
gene <- pg$`Gene names`
gene_first <- sapply(strsplit(ifelse(is.na(gene), "", gene), ";"), `[`, 1)
gene_first[is.na(gene_first) | gene_first == ""] <-
  sapply(strsplit(pg$`Protein IDs`[is.na(gene_first) | gene_first == ""], ";"), `[`, 1)
rownames(mat) <- make.unique(gene_first)

annot <- pg %>%
  transmute(row_id = rownames(mat),
            protein_ids = `Protein IDs`,
            gene_names  = `Gene names`,
            protein_names = `Protein names`,
            unique_peptides = `Unique peptides`)

# ---- 4. Missingness overview (before filtering) -------------------
cat("Missing values overall:", round(100 * mean(is.na(mat)), 1), "%\n")
miss_by_sample <- colMeans(is.na(mat))
print(summary(miss_by_sample))

# ---- 5. Valid-value filter ----------------------------------------
# Keep a protein if, in at least one group, >= min_valid_frac of samples
# have a quantified value. Use `group` (CON/AS/MR); swap for `group_sex`
# if you want the filter to respect the disease x sex cells.
grp <- samples$group
valid_frac <- sapply(levels(grp), function(g) {
  rowMeans(!is.na(mat[, grp == g, drop = FALSE]))
})
keep <- rowSums(valid_frac >= min_valid_frac) >= 1

mat_f   <- mat[keep, ]
annot_f <- annot[keep, ]
cat("Proteins kept after valid-value filter:", sum(keep), "of", length(keep), "\n")

# ---- 6. Optional: median normalisation per sample -----------------
# LFQ values are already normalised by MaxQuant; check before adding more.
# Uncomment to centre every sample on the global median:
# med <- apply(mat_f, 2, median, na.rm = TRUE)
# mat_f <- sweep(mat_f, 2, med - mean(med), "-")

# ---- 7. Quick QC plots --------------------------------------------
# Intensity distributions per sample
as.data.frame(mat_f) %>%
  rownames_to_column("protein") %>%
  pivot_longer(-protein, names_to = "sample", values_to = "log2LFQ") %>%
  left_join(samples, by = "sample") %>%
  ggplot(aes(sample, log2LFQ, fill = group)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6)) +
  labs(title = "log2 LFQ per sample (after filtering)")

# Sex sanity check: Y-linked proteins should separate males from females
y_genes <- intersect(c("DDX3Y", "EIF1AY", "UTY", "KDM5D", "RPS4Y1", "USP9Y"),
                     rownames(mat_f))
print(y_genes)
if (length(y_genes) > 0) {
  as.data.frame(t(mat_f[y_genes, , drop = FALSE])) %>%
    rownames_to_column("sample") %>%
    left_join(samples, by = "sample") %>%
    pivot_longer(all_of(y_genes), names_to = "protein", values_to = "log2LFQ") %>%
    ggplot(aes(sex, log2LFQ, colour = group)) +
    geom_jitter(width = 0.15) +
    facet_wrap(~protein) +
    theme_bw() +
    labs(title = "Y-linked proteins by sex (check for sample swaps)")
}

# PCA on complete cases only (no imputation) to look for outliers / batch
cc  <- mat_f[complete.cases(mat_f), ]
pca <- prcomp(t(cc), scale. = TRUE)
pve <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("sample") %>%
  left_join(samples, by = "sample") %>%
  ggplot(aes(PC1, PC2, colour = group, shape = sex)) +
  geom_point(size = 3) +
  theme_bw() +
  labs(x = paste0("PC1 (", pve[1], "%)"), y = paste0("PC2 (", pve[2], "%)"),
       title = paste("PCA on", nrow(cc), "complete proteins"))

# ---- 8. Save ------------------------------------------------------
saveRDS(list(mat = mat_f, annot = annot_f, samples = samples),
        "PR1_filtered.rds")
write.csv(cbind(annot_f, mat_f), "PR1_filtered_log2LFQ.csv", row.names = FALSE)

# ---- 9. Next step (limma) -----------------------------------------
# library(limma)
# design <- model.matrix(~ group * sex, data = samples)
# fit <- lmFit(mat_f, design) |> eBayes()
# topTable(fit, coef = "groupAS:sexmale")   # disease x sex interaction
