library(TCGAbiolinks)
library(tidyverse)

dir.create("results", showWarnings = FALSE)

# 1. Ask the GDC which breast cancer 450k methylation files exist
query <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "DNA Methylation",
  platform      = "Illumina Human Methylation 450",
  data.type     = "Methylation Beta Value"
)
files <- getResults(query)

cat("Files:", nrow(files), "\n")
print(table(files$sample_type))
cat("Unique patients:", n_distinct(files$cases.submitter_id), "\n")

set.seed(42)

# 1. Keep primary tumours and normals; one file per patient per tissue type
files <- files %>%
  filter(sample_type %in% c("Primary Tumor", "Solid Tissue Normal")) %>%
  distinct(cases.submitter_id, sample_type, .keep_all = TRUE)

normals <- files %>% filter(sample_type == "Solid Tissue Normal")
tumours <- files %>% filter(sample_type == "Primary Tumor")

# 2. Patients with both a tumour and a normal sample
paired_ids <- intersect(normals$cases.submitter_id, tumours$cases.submitter_id)
cat("Patients with matched tumour and normal:", length(paired_ids), "\n")

# 3. Select: all normals, all paired tumours, then random tumours up to 200
paired_tumours <- tumours %>% filter(cases.submitter_id %in% paired_ids)
other_tumours <- tumours %>%
  filter(!cases.submitter_id %in% paired_ids) %>%
  slice_sample(n = max(0, 200 - nrow(paired_tumours)))
selected <- bind_rows(normals, paired_tumours, other_tumours)

print(table(selected$sample_type))
cat("Barcode example:", head(selected$cases, 2), "\n")

dir.create("data", showWarnings = FALSE)

# 1. Save the list of selected samples (small file, goes on GitHub)
write.csv(selected[, c("cases", "cases.submitter_id", "sample_type")],
          "results/selected_samples.csv", row.names = FALSE)

# 2. Query only the selected samples
query_sel <- GDCquery(
  project       = "TCGA-BRCA",
  data.category = "DNA Methylation",
  platform      = "Illumina Human Methylation 450",
  data.type     = "Methylation Beta Value",
  barcode       = selected$cases
)

# 3. Download, then build the data object
GDCdownload(query_sel, method = "api", files.per.chunk = 20)
met <- GDCprepare(query_sel, save = TRUE,
                  save.filename = "data/brca_methylation_selected.rda")

# 4. Checks (do not print the full matrix)
cat("Probes x samples:", dim(met), "\n")
print(table(colData(met)$sample_type))
beta <- assay(met)
cat("Missing values:", sum(is.na(beta)), "\n")
cat("Beta range:", range(beta, na.rm = TRUE), "\n")

library(TCGAbiolinks)
sesameData::sesameDataCache()   # small annotation files, needs internet

met <- GDCprepare(query_sel, save = TRUE,
                  save.filename = "data/brca_methylation_selected.rda")

cat("Probes x samples:", dim(met), "\n")
print(table(colData(met)$sample_type))
beta <- assay(met)
cat("Missing values:", sum(is.na(beta)), "\n")
cat("Beta range:", range(beta, na.rm = TRUE), "\n")

met <- GDCprepare(query_sel, save = TRUE,
                  save.filename = "data/brca_methylation_selected.rda")

cat("Probes x samples:", dim(met), "\n")
print(table(colData(met)$sample_type))
beta <- assay(met)
cat("Missing values:", sum(is.na(beta)), "\n")
cat("Beta range:", range(beta, na.rm = TRUE), "\n")

library(SummarizedExperiment)

print(table(colData(met)$sample_type))
beta <- assay(met)
cat("Missing values:", sum(is.na(beta)), "\n")
cat("Probes with any missing value:", sum(rowSums(is.na(beta)) > 0), "\n")
cat("Beta range:", range(beta, na.rm = TRUE), "\n")
print(head(rownames(beta), 3))
print(colnames(rowData(met)))

library(SummarizedExperiment)

rd <- as.data.frame(rowData(met))
cat("Probe types:\n");  print(table(rd$probeType, useNA = "ifany"))
cat("MASK_general:\n"); print(table(rd$MASK_general, useNA = "ifany"))
cat("Chromosomes:\n");  print(head(sort(table(rd$chrm_A), decreasing = TRUE), 30))

# Filters
na_frac <- rowMeans(is.na(beta))
keep <- rd$probeType == "cg" &
  !as.logical(rd$MASK_general) &
  !(rd$chrm_A %in% c("chrX", "chrY", "chrM")) &
  !is.na(rd$chrm_A) &
  na_frac <= 0.05
keep[is.na(keep)] <- FALSE

cat("Probes before filtering:", nrow(beta), "\n")
cat("Probes kept:", sum(keep), "\n")

beta_f <- beta[keep, ]
rd_f   <- rd[keep, ]
cat("Missing values left:", sum(is.na(beta_f)), "\n")

library(limma)
library(tidyverse)
dir.create("results", showWarnings = FALSE)

# 1. Groups (normal is the reference level)
group <- factor(ifelse(colData(met)$sample_type == "Primary Tumor", "Tumour", "Normal"),
                levels = c("Normal", "Tumour"))
print(table(group))

# 2. Beta -> M-values (clip to avoid infinite values)
b <- pmin(pmax(beta_f, 0.001), 0.999)
m_vals <- log2(b / (1 - b))

# 3. limma on M-values
design <- model.matrix(~ group)
fit <- eBayes(lmFit(m_vals, design))
res <- topTable(fit, coef = "groupTumour", number = Inf, sort.by = "none")

# 4. Effect size on the beta scale
res$mean_normal <- rowMeans(beta_f[, group == "Normal"], na.rm = TRUE)
res$mean_tumour <- rowMeans(beta_f[, group == "Tumour"], na.rm = TRUE)
res$delta_beta  <- res$mean_tumour - res$mean_normal
res$probe <- rownames(res)
res$gene  <- rd_f$gene_HGNC
res$chr   <- rd_f$chrm_A

res$status <- case_when(
  res$adj.P.Val < 0.05 & res$delta_beta >=  0.2 ~ "Hypermethylated in tumour",
  res$adj.P.Val < 0.05 & res$delta_beta <= -0.2 ~ "Hypomethylated in tumour",
  TRUE ~ "Not significant / small change"
)

# 5. Checks
cat("Probes tested:", nrow(res), "\n")
print(table(res$status))
cat("Significant by adjusted p only:", sum(res$adj.P.Val < 0.05), "\n")
print(head(res %>% arrange(adj.P.Val) %>%
             dplyr::select(probe, gene, delta_beta, adj.P.Val), 10))

# 6. Save (the full table is large, so it stays in data/ and off GitHub)
saveRDS(res, "data/methylation_results_full.rds")
top <- res %>% filter(status != "Not significant / small change") %>%
  arrange(desc(abs(delta_beta))) %>% slice_head(n = 1000)
write.csv(top, "results/top1000_differential_cpgs.csv", row.names = FALSE)

# 7. Volcano
volcano <- ggplot(res, aes(delta_beta, -log10(adj.P.Val), colour = status)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_vline(xintercept = c(-0.2, 0.2), linetype = "dashed") +
  scale_colour_manual(values = c("Hypermethylated in tumour" = "firebrick",
                                 "Hypomethylated in tumour" = "steelblue",
                                 "Not significant / small change" = "grey75")) +
  labs(title = "Differential DNA methylation: primary tumour vs normal (TCGA-BRCA)",
       x = "Difference in mean beta (tumour - normal)",
       y = "-log10 adjusted p-value", colour = NULL) +
  theme_classic() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))
print(volcano)
ggsave("results/volcano_methylation.png", volcano, width = 8, height = 6, dpi = 300, bg = "white")

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
res$island <- ann[res$probe, "Relation_to_Island"]
res$island[res$island == ""] <- NA
res$island <- factor(res$island,
                     levels = c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea"))
cat("Probes without island annotation:", sum(is.na(res$island)), "\n")

# Share of each island category within each status group
island_tab <- res %>%
  filter(!is.na(island)) %>%
  count(status, island) %>%
  group_by(status) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()
print(island_tab %>% dplyr::select(status, island, pct) %>%
        pivot_wider(names_from = island, values_from = pct) %>%
        mutate(across(where(is.numeric), ~ round(.x, 1))))

island_plot <- ggplot(island_tab, aes(island, pct, fill = status)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Hypermethylated in tumour" = "firebrick",
                               "Hypomethylated in tumour" = "steelblue",
                               "Not significant / small change" = "grey70")) +
  labs(title = "CpG island context of differentially methylated probes",
       x = "Position relative to CpG island", y = "% of probes in group", fill = NULL) +
  theme_classic()
print(island_plot)
ggsave("results/island_context.png", island_plot, width = 8, height = 5, dpi = 300, bg = "white")
write.csv(island_tab, "results/island_context_table.csv", row.names = FALSE)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

ann <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
res$island <- ann[res$probe, "Relation_to_Island"]
res$island[res$island == ""] <- NA
res$island <- factor(res$island,
                     levels = c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea"))
cat("Probes without island annotation:", sum(is.na(res$island)), "\n")

# Share of each island category within each status group
island_tab <- res %>%
  dplyr::filter(!is.na(island)) %>%
  dplyr::count(status, island) %>%
  dplyr::group_by(status) %>%
  dplyr::mutate(pct = 100 * n / sum(n)) %>%
  dplyr::ungroup()
print(island_tab %>% dplyr::select(status, island, pct) %>%
        tidyr::pivot_wider(names_from = island, values_from = pct) %>%
        dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 1))))

island_plot <- ggplot(island_tab, aes(island, pct, fill = status)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Hypermethylated in tumour" = "firebrick",
                               "Hypomethylated in tumour" = "steelblue",
                               "Not significant / small change" = "grey70")) +
  labs(title = "CpG island context of differentially methylated probes",
       x = "Position relative to CpG island", y = "% of probes in group", fill = NULL) +
  theme_classic()
print(island_plot)
ggsave("results/island_context.png", island_plot, width = 8, height = 5, dpi = 300, bg = "white")
write.csv(island_tab, "results/island_context_table.csv", row.names = FALSE)

library(cBioPortalData)
library(MultiAssayExperiment)
library(tidyverse)

# 0. Free memory from objects we no longer need
rm(b, m_vals, fit, beta); invisible(gc())

# 1. Expression data (log2 RSEM + 1)
brca <- cBioDataPack("brca_tcga_pan_can_atlas_2018", ask = FALSE)
expr_raw <- brca[["mrna_seq_v2_rsem"]]
expr_mat <- if (is.matrix(expr_raw)) expr_raw else assay(expr_raw)
expr_mat <- expr_mat[!is.na(rownames(expr_mat)) & rownames(expr_mat) != "", ]
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]
expr_log <- log2(expr_mat + 1)

# 2. Match tumours that have both methylation and expression
tum_cols <- which(group == "Tumour")
tum_ids  <- substr(colnames(beta_f)[tum_cols], 1, 15)
use <- tum_ids %in% colnames(expr_log)
cat("Tumours with methylation and expression:", sum(use), "\n")
meth_idx <- tum_cols[use]
expr_sub <- expr_log[, tum_ids[use]]

# 3. Probe -> gene and promoter annotation
ann_sub <- ann[res$probe, c("UCSC_RefGene_Name", "UCSC_RefGene_Group")]
res$ref_gene <- sub(";.*", "", ann_sub$UCSC_RefGene_Name)
res$ref_gene[res$ref_gene == ""] <- NA
res$promoter <- grepl("TSS200|TSS1500", ann_sub$UCSC_RefGene_Group)

ok <- !is.na(res$ref_gene) & res$promoter & res$ref_gene %in% rownames(expr_log)
hyper_probes <- res$probe[ok & res$status == "Hypermethylated in tumour"]
bg_all <- res$probe[ok & res$status == "Not significant / small change"]
set.seed(42)
bg_probes <- sample(bg_all, min(5000, length(bg_all)))
cat("Hypermethylated promoter probes:", length(hyper_probes),
    " Background probes:", length(bg_probes), "\n")

# 4. Spearman correlation of methylation vs expression across tumours
gene_of <- setNames(res$ref_gene, res$probe)
cor_one <- function(p) {
  x <- beta_f[p, meth_idx]
  y <- expr_sub[unname(gene_of[p]), ]
  if (sum(complete.cases(x, y)) < 30 || sd(y) == 0) return(NA_real_)
  suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
}

cors <- tibble(
  probe = c(hyper_probes, bg_probes),
  group = rep(c("Hypermethylated promoter probes", "Background promoter probes"),
              c(length(hyper_probes), length(bg_probes)))
) %>%
  dplyr::mutate(gene = unname(gene_of[probe]),
                rho = unname(vapply(probe, cor_one, numeric(1)))) %>%
  dplyr::filter(!is.na(rho))

print(cors %>% dplyr::group_by(group) %>%
        dplyr::summarise(n = dplyr::n(),
                         median_rho = round(median(rho), 3),
                         pct_rho_below_neg03 = round(100 * mean(rho < -0.3), 1)))
print(wilcox.test(rho ~ group, data = cors))

# 5. Plot and save
cor_plot <- ggplot(cors, aes(rho, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Hypermethylated promoter probes" = "firebrick",
                               "Background promoter probes" = "grey60")) +
  labs(title = "Promoter methylation vs expression of the same gene (tumours)",
       x = "Spearman correlation (methylation vs expression)",
       y = "Density", fill = NULL) +
  theme_classic()
print(cor_plot)
ggsave("results/methylation_vs_expression.png", cor_plot, width = 8, height = 5, dpi = 300, bg = "white")
write.csv(cors, "results/methylation_expression_correlations.csv", row.names = FALSE)

print(head(cors %>% dplyr::filter(group == "Hypermethylated promoter probes") %>%
             dplyr::arrange(rho) %>% dplyr::select(probe, gene, rho), 15))

