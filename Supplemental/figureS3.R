suppressPackageStartupMessages({
  library(tidyverse)
})
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Paths and colors
root_dir <- file.path("..","Data")
print(root_dir)
getwd()
list.files(root_dir)

flow_path <- file.path(root_dir, "markov_clinical.csv")
out_dir <- file.path(root_dir, "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

pal_fill <- c(BM = "#CBF2CE", PB = "#DDC0C2", RM = "#B410D1", RL = "#37A749")

# Load and reshape
flow_raw <- read.csv(flow_path, check.names = FALSE, stringsAsFactors = FALSE)
rn <- flow_raw[[1]]
mat <- flow_raw[, -1, drop = FALSE]
flow_t <- as.data.frame(t(as.matrix(mat)), stringsAsFactors = FALSE)
colnames(flow_t) <- rn
flow_t <- tibble::rownames_to_column(flow_t, var = "Sample")

# Numeric features
M_cols <- grep("^M\\d+$", names(flow_t), value = TRUE)
flow_t[M_cols] <- lapply(flow_t[M_cols], function(v) suppressWarnings(as.numeric(v)))


# Specimen: BM vs PB
spec_long <- flow_t |>
  select(Sample, Specimen, all_of(M_cols)) |>
  filter(!is.na(Specimen), Specimen %in% c("BM", "PB")) |>
  pivot_longer(cols = all_of(M_cols), names_to = "Feature", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(
    Specimen = factor(Specimen, levels = c("BM", "PB")),
    Feature  = factor(Feature, levels = M_cols)
  )


# ----- Four separate figures: BM, PB, RM, RL -----
# BM only
bm_long <- spec_long %>% filter(Specimen == "BM")
p_bm <- ggplot(bm_long, aes(x = Feature, y = value)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.85,
               fill = pal_fill["BM"], color = "gray20") +
  geom_point(color = pal_fill["BM"],
             position = position_jitter(width = 0.15),
             alpha = 0.6, size = 1) +
  labs(x = "", y = "", title = "") +
  scale_x_discrete(labels = NULL) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

print(p_bm)

ggsave(file.path(out_dir, "Markov_box_BM_only.pdf"), p_bm, width = 3, height = 3.35, units = "in")

# PB only
pb_long <- spec_long %>% filter(Specimen == "PB")
p_pb <- ggplot(pb_long, aes(x = Feature, y = value)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.85,
               fill = pal_fill["PB"], color = "gray20") +
  geom_point(color = pal_fill["PB"],
             position = position_jitter(width = 0.15),
             alpha = 0.6, size = 1) +
  labs(x = "", y = "", title = "") +
  scale_x_discrete(labels = NULL) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

print(p_pb)
ggsave(file.path(out_dir, "flow_box_PB_only.pdf"), p_pb, width = 3, height = 3.35, units = "in")

# RM (Remission) from DiseaseState
rm_long <- flow_t |>
  select(Sample, DiseaseState, all_of(M_cols)) |>
  filter(!is.na(DiseaseState), DiseaseState == "Remission") |>
  pivot_longer(cols = all_of(M_cols), names_to = "Feature", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(
    Feature = factor(Feature, levels = M_cols)
  )

p_rm <- ggplot(rm_long, aes(x = Feature, y = value)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.85,
               fill = pal_fill["RM"], color = "gray20") +
  geom_point(color = pal_fill["RM"],
             position = position_jitter(width = 0.15),
             alpha = 0.6, size = 1) +
  labs(x = "", y = "", title = "") +
  scale_x_discrete(labels = NULL) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.y  = element_text(size = 12, face = "bold"),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none"
  )
print(p_rm)
ggsave(file.path(out_dir, "flow_box_RM_only.pdf"), p_rm, width = 3, height = 3.35, units = "in")


# Relapse from Relapse == 1
rl_long <- flow_t |>
  mutate(Relapse = suppressWarnings(as.numeric(Relapse))) |>
  select(Sample, Relapse, all_of(M_cols)) |>
  filter(!is.na(Relapse), Relapse == 1) |>
  pivot_longer(cols = all_of(M_cols), names_to = "Feature", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(
    Feature = factor(Feature, levels = M_cols)
  )

p_rl <- ggplot(rl_long, aes(x = Feature, y = value)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.85,
               fill = pal_fill["RL"], color = "gray20") +
  geom_point(color = pal_fill["RL"],
             position = position_jitter(width = 0.15),
             alpha = 0.6, size = 1) +
  labs(x = "", y = "", title = "") +
  scale_x_discrete(labels = NULL) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )
print(p_rl)
ggsave(file.path(out_dir, "flow_box_RL_only.pdf"), p_rl, width = 3, height = 3.35, units = "in")


# Genetics: BCR::ABL1 vs Negative (exclude BCR::ABL1-like)
bcrneg_long <- flow_t |>
  select(Sample, `BCR::ABL1/like`, all_of(M_cols)) |>
  filter(!is.na(`BCR::ABL1/like`), `BCR::ABL1/like` %in% c("BCR::ABL1", "Negative")) |>
  pivot_longer(cols = all_of(M_cols), names_to = "Feature", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(
    Feature = factor(Feature, levels = M_cols),
    Genetics = factor(`BCR::ABL1/like`, levels = c("BCR::ABL1", "Negative"))
  )

# BCR::ABL1 only
bcr_only <- bcrneg_long %>% filter(Genetics == "BCR::ABL1")
p_bcr <- ggplot(bcr_only, aes(x = Feature, y = value)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.85,
               fill = "#4A85BD", color = "gray20") +
  geom_point(color = "#4A85BD",
             position = position_jitter(width = 0.15),
             alpha = 0.6, size = 1) +
  labs(x = "", y = "", title = "") +
  scale_x_discrete(labels = NULL) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.y  = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )
print(p_bcr)
ggsave(file.path(out_dir, "flow_box_BCR_ABL1_only.pdf"), p_bcr, width = 3, height = 3.35, units = "in")

# Negative only
neg_only <- bcrneg_long %>% filter(Genetics == "Negative")
p_neg <- ggplot(neg_only, aes(x = Feature, y = value)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.85,
               fill = "#DC362C", color = "gray20") +
  geom_point(color = "#DC362C",
             position = position_jitter(width = 0.15),
             alpha = 0.6, size = 1) +
  labs(x = "", y = "", title = "") +
  scale_x_discrete(labels = NULL) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.y  = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )
print(p_neg)
ggsave(file.path(out_dir, "flow_box_BCRlike_Negative_only.pdf"), p_neg, width = 3, height = 3.35, units = "in")

# --- Genetics x Specimen subsets ---

# BCR::ABL1 & BM
bcr_bm_long <- flow_t |>
  select(Sample, Specimen, `BCR::ABL1/like`, all_of(M_cols)) |>
  filter(!is.na(Specimen), Specimen == "BM",
         !is.na(`BCR::ABL1/like`), `BCR::ABL1/like` == "BCR::ABL1") |>
  pivot_longer(cols = all_of(M_cols), names_to = "Feature", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(Feature = factor(Feature, levels = M_cols))

p_bcr_bm <- ggplot(bcr_bm_long, aes(x = Feature, y = value)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.85,
               fill = "#4A85BD", color = "gray20") +
  geom_point(color = "#4A85BD",
             position = position_jitter(width = 0.15),
             alpha = 0.6, size = 1) +
  labs(x = "", y = "", title = "") +
  scale_x_discrete(labels = NULL) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )
print(p_bcr_bm)
ggsave(file.path(out_dir, "flow_box_BCR_ABL1_BM.pdf"), p_bcr_bm, width = 3, height = 3.35, units = "in")

# BCR::ABL1 & PB
bcr_pb_long <- flow_t |>
  select(Sample, Specimen, `BCR::ABL1/like`, all_of(M_cols)) |>
  filter(!is.na(Specimen), Specimen == "PB",
         !is.na(`BCR::ABL1/like`), `BCR::ABL1/like` == "BCR::ABL1") |>
  pivot_longer(cols = all_of(M_cols), names_to = "Feature", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(Feature = factor(Feature, levels = M_cols))

p_bcr_pb <- ggplot(bcr_pb_long, aes(x = Feature, y = value)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.85,
               fill = "#4A85BD", color = "gray20") +
  geom_point(color = "#4A85BD",
             position = position_jitter(width = 0.15),
             alpha = 0.6, size = 1) +
  labs(x = "", y = "", title = "") +
  scale_x_discrete(labels = NULL) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )
print(p_bcr_pb)
ggsave(file.path(out_dir, "flow_box_BCR_ABL1_PB.pdf"), p_bcr_pb, width = 3, height = 3.35, units = "in")

# Negative & BM
neg_bm_long <- flow_t |>
  select(Sample, Specimen, `BCR::ABL1/like`, all_of(M_cols)) |>
  filter(!is.na(Specimen), Specimen == "BM",
         !is.na(`BCR::ABL1/like`), `BCR::ABL1/like` == "Negative") |>
  pivot_longer(cols = all_of(M_cols), names_to = "Feature", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(Feature = factor(Feature, levels = M_cols))

p_neg_bm <- ggplot(neg_bm_long, aes(x = Feature, y = value)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.85,
               fill = "#DC362C", color = "gray20") +
  geom_point(color = "#DC362C",
             position = position_jitter(width = 0.15),
             alpha = 0.6, size = 1) +
  labs(x = "", y = "", title = "") +
  scale_x_discrete(labels = NULL) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )
print(p_neg_bm)
ggsave(file.path(out_dir, "flow_box_Negative_BM.pdf"), p_neg_bm, width = 3, height = 3.35, units = "in")

# Negative & PB
neg_pb_long <- flow_t |>
  select(Sample, Specimen, `BCR::ABL1/like`, all_of(M_cols)) |>
  filter(!is.na(Specimen), Specimen == "PB",
         !is.na(`BCR::ABL1/like`), `BCR::ABL1/like` == "Negative") |>
  pivot_longer(cols = all_of(M_cols), names_to = "Feature", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(Feature = factor(Feature, levels = M_cols))

p_neg_pb <- ggplot(neg_pb_long, aes(x = Feature, y = value)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, alpha = 0.85,
               fill = "#DC362C", color = "gray20") +
  geom_point(color = "#DC362C",
             position = position_jitter(width = 0.15),
             alpha = 0.6, size = 1) +
  labs(x = "", y = "", title = "") +
  scale_x_discrete(labels = NULL) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )
print(p_neg_pb)
ggsave(file.path(out_dir, "flow_box_Negative_PB.pdf"), p_neg_pb, width = 3, height = 3.35, units = "in")

# Per-Feature summary by Specimen (BM/PB) on the same data used for the plots
summary_tbl <- spec_long %>%
  group_by(Specimen, Feature) %>%
  summarize(
    n      = dplyr::n(),
    mean   = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd     = sd(value, na.rm = TRUE),
    IQR    = IQR(value, na.rm = TRUE),
    q1     = quantile(value, 0.25, na.rm = TRUE),
    q3     = quantile(value, 0.75, na.rm = TRUE),
    min    = min(value, na.rm = TRUE),
    max    = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Feature, Specimen)

print(summary_tbl, n = 50)  # adjust n as needed
write.csv(summary_tbl, file.path(out_dir, "specimen_feature_summary.csv"), row.names = FALSE)

# If you also want BM-only and PB-only quick tables:
bm_summary <- bm_long %>%
  group_by(Feature) %>%
  summarize(n = dplyr::n(),
            mean = mean(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            IQR = IQR(value, na.rm = TRUE),
            q1 = quantile(value, 0.25, na.rm = TRUE),
            q3 = quantile(value, 0.75, na.rm = TRUE),
            min = min(value, na.rm = TRUE),
            max = max(value, na.rm = TRUE),
            .groups = "drop")

pb_summary <- pb_long %>%
  group_by(Feature) %>%
  summarize(n = dplyr::n(),
            mean = mean(value, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            IQR = IQR(value, na.rm = TRUE),
            q1 = quantile(value, 0.25, na.rm = TRUE),
            q3 = quantile(value, 0.75, na.rm = TRUE),
            min = min(value, na.rm = TRUE),
            max = max(value, na.rm = TRUE),
            .groups = "drop")

# --- Stat summaries for RM, RL, BCR::ABL1, Negative ---
rm_summary <- rm_long %>%
  group_by(Feature) %>%
  summarize(
    n      = dplyr::n(),
    mean   = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd     = sd(value, na.rm = TRUE),
    IQR    = IQR(value, na.rm = TRUE),
    q1     = quantile(value, 0.25, na.rm = TRUE),
    q3     = quantile(value, 0.75, na.rm = TRUE),
    min    = min(value, na.rm = TRUE),
    max    = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(Feature)

rl_summary <- rl_long %>%
  group_by(Feature) %>%
  summarize(
    n      = dplyr::n(),
    mean   = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd     = sd(value, na.rm = TRUE),
    IQR    = IQR(value, na.rm = TRUE),
    q1     = quantile(value, 0.25, na.rm = TRUE),
    q3     = quantile(value, 0.75, na.rm = TRUE),
    min    = min(value, na.rm = TRUE),
    max    = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(Feature)

bcr_summary <- bcrneg_long %>%
  filter(Genetics == "BCR::ABL1") %>%
  group_by(Feature) %>%
  summarize(
    n      = dplyr::n(),
    mean   = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd     = sd(value, na.rm = TRUE),
    IQR    = IQR(value, na.rm = TRUE),
    q1     = quantile(value, 0.25, na.rm = TRUE),
    q3     = quantile(value, 0.75, na.rm = TRUE),
    min    = min(value, na.rm = TRUE),
    max    = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(Feature)

neg_summary <- bcrneg_long %>%
  filter(Genetics == "Negative") %>%
  group_by(Feature) %>%
  summarize(
    n      = dplyr::n(),
    mean   = mean(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    sd     = sd(value, na.rm = TRUE),
    IQR    = IQR(value, na.rm = TRUE),
    q1     = quantile(value, 0.25, na.rm = TRUE),
    q3     = quantile(value, 0.75, na.rm = TRUE),
    min    = min(value, na.rm = TRUE),
    max    = max(value, na.rm = TRUE),
    .groups = "drop"
  ) %>% arrange(Feature)

# Combined table for convenience
rm_summary$Group  <- "RM"
rl_summary$Group  <- "RL"
bcr_summary$Group <- "BCR::ABL1"
neg_summary$Group <- "Negative"

all_fig_summary <- bind_rows(rm_summary, rl_summary, bcr_summary, neg_summary) %>%
  relocate(Group)

print(all_fig_summary, n = 50)

# Save CSVs
write.csv(all_fig_summary, file.path(out_dir, "figure_groups_feature_summary.csv"), row.names = FALSE)
write.csv(rm_summary,  file.path(out_dir, "rm_feature_summary.csv"),        row.names = FALSE)
write.csv(rl_summary,  file.path(out_dir, "rl_feature_summary.csv"),        row.names = FALSE)
write.csv(bcr_summary, file.path(out_dir, "bcrabl1_feature_summary.csv"),   row.names = FALSE)
write.csv(neg_summary, file.path(out_dir, "negative_feature_summary.csv"),  row.names = FALSE)

# --- Stat summaries for Genetics x Specimen (BCR::ABL1/Negative x BM/PB) ---

summarize_long <- function(df) {
  df %>%
    group_by(Feature) %>%
    summarize(
      n      = dplyr::n(),
      mean   = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      sd     = sd(value, na.rm = TRUE),
      IQR    = IQR(value, na.rm = TRUE),
      q1     = quantile(value, 0.25, na.rm = TRUE),
      q3     = quantile(value, 0.75, na.rm = TRUE),
      min    = min(value, na.rm = TRUE),
      max    = max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>% arrange(Feature)
}

bcr_bm_summary <- summarize_long(bcr_bm_long); bcr_bm_summary$Group <- "BCR::ABL1_BM"
bcr_pb_summary <- summarize_long(bcr_pb_long); bcr_pb_summary$Group <- "BCR::ABL1_PB"
neg_bm_summary <- summarize_long(neg_bm_long); neg_bm_summary$Group <- "Negative_BM"
neg_pb_summary <- summarize_long(neg_pb_long); neg_pb_summary$Group <- "Negative_PB"

gen_spec_summary <- bind_rows(bcr_bm_summary, bcr_pb_summary, neg_bm_summary, neg_pb_summary) %>%
  relocate(Group)

print(gen_spec_summary, n = 50)

# Save CSVs
write.csv(gen_spec_summary, file.path(out_dir, "genetics_specimen_feature_summary.csv"), row.names = FALSE)
write.csv(bcr_bm_summary, file.path(out_dir, "bcrabl1_bm_feature_summary.csv"), row.names = FALSE)
write.csv(bcr_pb_summary, file.path(out_dir, "bcrabl1_pb_feature_summary.csv"), row.names = FALSE)
write.csv(neg_bm_summary, file.path(out_dir, "negative_bm_feature_summary.csv"), row.names = FALSE)
write.csv(neg_pb_summary, file.path(out_dir, "negative_pb_feature_summary.csv"), row.names = FALSE)
