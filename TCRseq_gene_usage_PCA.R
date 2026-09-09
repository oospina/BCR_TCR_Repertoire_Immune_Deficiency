# TCR metrics

library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)

setwd("")

# List all csv files
files <- list.files(pattern = "\\.csv$", full.names = TRUE)

# Read and combine them
df_all <- bind_rows(
  lapply(files, function(x) {
    df <- read.csv(x)
    df$sample <- tools::file_path_sans_ext(basename(x))
    return(df)
  })
)


clone_ratio <- df_all %>%
  group_by(sample) %>%
  summarise(
    unique_clones = n_distinct(clone_id),
    total_sequences = sum(duplicate_count),
    ratio_unique_total = unique_clones / total_sequences
  )

clone_ratio

write.csv(clone_ratio, "clone_ratio_TRB.csv", row.names = T)


# Gene usages


df_all$V_gene <- ifelse(
  grepl("TRBV", df_all$v_call),
  sub(".*?(TRBV[0-9]+(-[0-9]+)?).*", "\\1", df_all$v_call),
  NA
)
df_all$D_gene <- sub(".* (TRBD[^* ]+).*", "\\1", df_all$d_call)

# There are 2 clusters J gene segments


df_all$J_gene <- gsub(".* (TRBJ[0-9]+-[0-9]+).*", "\\1", df_all$j_call)

# Heavy gene usages (unique clones makes more sense for TCR)


V_usage <- df_all %>%
  filter(!is.na(V_gene)) %>%
  group_by(sample, V_gene) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

D_usage <- df_all %>%
  filter(!is.na(D_gene)) %>%
  group_by(sample, D_gene) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

J_usage <- df_all %>%
  filter(!is.na(J_gene)) %>%
  group_by(sample, J_gene) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

# Combined

heavy_usage <- bind_rows(
  V_usage %>% mutate(chain="V", gene = V_gene) %>% select(sample, chain, gene, count, freq),
  D_usage %>% mutate(chain="D", gene = D_gene) %>% select(sample, chain, gene, count, freq),
  J_usage %>% mutate(chain="J", gene = J_gene) %>% select(sample, chain, gene, count, freq)
)

write.csv(heavy_usage, "heavy_segments_TRB.csv", row.names = T)

write.csv(df_all, "combined_df_allclones_TRB.csv", row.names = T)


colnames(df_all)

trbv_order <- c(
  "TRBV1",
  "TRBV2",
  "TRBV3-1",
  
  "TRBV4-1", "TRBV4-2", "TRBV4-3",
  
  "TRBV5-1", "TRBV5-3", "TRBV5-4", "TRBV5-5", "TRBV5-6", "TRBV5-7", "TRBV5-8",
  
  "TRBV6-1", "TRBV6-2", "TRBV6-3", "TRBV6-4", "TRBV6-5", "TRBV6-6", "TRBV6-7", "TRBV6-8", "TRBV6-9",
  
  "TRBV7-1", "TRBV7-2", "TRBV7-3", "TRBV7-4", "TRBV7-5", "TRBV7-6", "TRBV7-7", "TRBV7-8", "TRBV7-9",
  
  "TRBV8-1", "TRBV8-2", "TRBV8-3",
  
  "TRBV9",
  
  "TRBV10-1", "TRBV10-2", "TRBV10-3",
  
  "TRBV11-1", "TRBV11-2", "TRBV11-3",
  
  "TRBV12-1", "TRBV12-2", "TRBV12-3", "TRBV12-4", "TRBV12-5",
  
  "TRBV13", "TRBV14", "TRBV15", "TRBV16", "TRBV17", "TRBV18", "TRBV19",
  
  "TRBV20-1",
  
  "TRBV21-1",
  
  "TRBV22-1",
  
  "TRBV23-1", "TRBV23-2",
  
  "TRBV24-1",
  
  "TRBV25-1",
  
  "TRBV26",
  
  "TRBV27", "TRBV28", "TRBV29-1",
  
  "TRBV30"
)



V_usage <- V_usage %>%
  complete(sample, V_gene = trbv_order, fill = list(count = 0, freq = 0))


V_usage$V_gene <- factor(V_usage$V_gene, levels = trbv_order)


ggplot(V_usage, aes(x = V_gene, y = sample, fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(
    colours = c("white", "blue", "yellow"),
    values = scales::rescale(c(0, 0.01, max(V_usage$freq, na.rm = TRUE))),
    na.value = "white",
    name = "Frequency"
  ) +
  labs(x = "V gene (distal → proximal)", y = "Sample") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Exclusion of unused genes

unused_V <- V_usage %>%
  group_by(V_gene) %>%
  summarise(total = sum(freq, na.rm = TRUE)) %>%
  filter(total == 0) %>%
  pull(V_gene)

unused_V
# Remove them

V_usage_filtered <- V_usage %>%
  filter(!V_gene %in% unused_V)

# update the order vector

V_order_filtered <- trbv_order[!trbv_order %in% unused_V]


V_usage_filtered <- V_usage_filtered[!is.na(V_usage_filtered$V_gene), ]
V_usage_filtered$V_gene <- factor(
  V_usage_filtered$V_gene,
  levels = V_order_filtered
)

unique(V_usage_filtered$V_gene)
range(V_usage_filtered$freq, na.rm = TRUE)


ggplot(V_usage_filtered, aes(x = V_gene, y = sample, fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    na.value = "white",
    name = "TRBV usage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Reordering and renaming the samples
V_usage_filtered$sample<- as.character(V_usage_filtered$sample)
unique(V_usage_filtered$sample)


# Define your desired order
sample_order <- c("P2_T2", "P6", "P1_T1", "P1_T2", 
                  "P10", "P11", "P4", "P15", "P13", 
                  "P12", "P2_T1", "Het1", "Het2", "Het4",
                  "Het5", "HD1", "HD9", "HD12")

V_usage_filtered <- V_usage_filtered %>%
  mutate(
    # Remove "_1indexXX_TRB_allclone" from all samples
    sample = sub("_\\d+index\\d+_TRB_allclone$", "", sample),
    
    # Fix HD samples: optional underscore + leading zeros
    sample = gsub("HD_?0*(\\d+)", "HD\\1", sample),
    
    # Fix P4 and P6 to remove _T1
    sample = case_when(
      sample == "P4_T1" ~ "P4",
      sample == "P6_T1" ~ "P6",
      TRUE ~ sample
    ),
    
    # Convert to factor with your exact order
    sample = factor(sample, levels = sample_order)
  )


# Check levels
levels(V_usage_filtered$sample)
unique(V_usage_filtered$sample)


p <- ggplot(V_usage_filtered, aes(x = V_gene, y = forcats::fct_rev(sample), fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#ffffcc", "#a1dab4", "#41b6c4", "#225ea8"),
    values = scales::rescale(c(0, 0.001, 0.01, 0.05, max(V_usage_filtered$freq, na.rm = TRUE))),
    na.value = "white",
    name = "TRBV usage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed(ratio = 1)

p
ggsave("V_gene_usage_unique_clones_TRB.png", p, bg = "white")


write.csv(V_usage_filtered, "V_usage_TRB_unique_clones.csv", row.names = T)

# D genes

D_order <- c("TRBD1", "TRBD2")

D_usage <- D_usage %>%
  complete(sample, D_gene = D_order, fill = list(count = 0, freq = 0))

D_usage$D_gene <- factor(D_usage$D_gene, levels = D_order)


ggplot(D_usage, aes(x = D_gene, y = sample, fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradient(
    low = "blue",
    high = "yellow",
    na.value = "white",
    name = "Frequency"
  ) +
  labs(x = "D gene (distal → proximal)", y = "Sample") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
unique(D_usage$D_gene)

range(D_usage$freq, na.rm = TRUE)


ggplot(D_usage, aes(x = D_gene, y = sample, fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    na.value = "white",
    name = "TRBD usage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Clean up the names

D_usage <- D_usage %>%
  mutate(
    # Remove "_1indexXX_TRB_allclone" from all samples
    sample = sub("_\\d+index\\d+_TRB_allclone$", "", sample),
    
    # Fix HD samples: optional underscore + leading zeros
    sample = gsub("HD_?0*(\\d+)", "HD\\1", sample),
    
    # Fix P4 and P6 to remove _T1
    sample = case_when(
      sample == "P4_T1" ~ "P4",
      sample == "P6_T1" ~ "P6",
      TRUE ~ sample
    ),
    
    # Convert to factor with your exact order
    sample = factor(sample, levels = sample_order)
  )

# Check levels
levels(D_usage$sample)
unique(D_usage$sample)


ggplot(D_usage, aes(x = D_gene, y = forcats::fct_rev(sample), fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#ffffcc", "#a1dab4", "#41b6c4", "#225ea8"),
    values = scales::rescale(c(0, 0.001, 0.01, 0.05, max(D_usage$freq, na.rm = TRUE))),
    na.value = "white",
    name = "TRBD usage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed(ratio = 1)


# Bigger fonts
p_d <- ggplot(D_usage, aes(x = D_gene, y = forcats::fct_rev(sample), fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#ffffcc", "#a1dab4", "#41b6c4", "#225ea8"),
    values = scales::rescale(c(0, 0.001, 0.01, 0.05, max(D_usage$freq, na.rm = TRUE))),
    na.value = "white",
    name = "TRBD usage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 11)) + 
  coord_fixed(ratio = 1)

p_d
D_usage
ggsave("D_gene_usage_unique_clones_TRB.png", p_d, bg = "white")

write.csv(D_usage, "D_usage_unique_clones_TRB.csv", row.names = T)

#J usage

J_order <- c("TRBJ1-1", "TRBJ1-2", "TRBJ1-3", "TRBJ1-4", "TRBJ1-5", "TRBJ1-6", "TRBJ2-1", "TRBJ2-2", "TRBJ2-3", "TRBJ2-4", "TRBJ2-5", "TRBJ2-6", "TRBJ2-7")

J_usage <- J_usage %>%
  complete(sample, J_gene = J_order, fill = list(count = 0, freq = 0))

J_usage$J_gene <- factor(J_usage$J_gene, levels = J_order)


ggplot(J_usage, aes(x = J_gene, y = sample, fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradient(
    low = "blue",
    high = "yellow",
    na.value = "white",
    name = "Frequency"
  ) +
  labs(x = "J gene (distal → proximal)", y = "Sample") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggplot(J_usage, aes(x = J_gene, y = sample, fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    na.value = "white",
    name = "TRBJ usage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Clean up the names

J_usage <- J_usage %>%
  mutate(
    # Remove "_1indexXX_TRB_allclone" from all samples
    sample = sub("_\\d+index\\d+_TRB_allclone$", "", sample),
    
    # Fix HD samples: optional underscore + leading zeros
    sample = gsub("HD_?0*(\\d+)", "HD\\1", sample),
    
    # Fix P4 and P6 to remove _T1
    sample = case_when(
      sample == "P4_T1" ~ "P4",
      sample == "P6_T1" ~ "P6",
      TRUE ~ sample
    ),
    
    # Convert to factor with your exact order
    sample = factor(sample, levels = sample_order)
  )


sample_order <- c("P2_T2", "P6", "P1_T1", "P1_T2", 
                  "P10", "P11", "P4", "P15", "P13", 
                  "P12", "P2_T1", "Het1", "Het2", "Het4",
                  "Het5", "HD1", "HD9", "HD12")



p_j <- ggplot(J_usage, aes(x = J_gene, y = forcats::fct_rev(sample), fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#ffffcc", "#a1dab4", "#41b6c4", "#225ea8"),
    values = scales::rescale(c(0, 0.001, 0.01, 0.05, max(J_usage$freq, na.rm = TRUE))),
    na.value = "white",
    name = "TRBJ usage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14)) + 
  coord_fixed(ratio = 1)

p_j
ggsave("J_gene_usage_unique_clones_TRB.png", p_j, bg = "white")

write.csv(J_usage, "J_usage_unique_clones_TRB.csv", row.names = T)

# Top 10 clones per sample

clone_freq <- df_all %>%
  group_by(sample, clone_id) %>%
  summarise(
    clone_size = sum(duplicate_count),
    .groups = "drop"
  ) %>%
  group_by(sample) %>%
  mutate(
    total_sequences = sum(clone_size),
    clone_freq_percent = 100 * clone_size / total_sequences
  ) %>%
  arrange(sample, desc(clone_freq_percent))

top10_clones <- clone_freq %>%
  group_by(sample) %>%
  slice_max(order_by = clone_freq_percent, n = 10) %>%
  ungroup()

top10_clones


write.csv(top10_clones, "top10_clones_frequency_per_patient_TRB.csv", row.names = T)

# cumulative frequency

cumulative_freq <- top10_clones %>%
  group_by(sample) %>%
  summarise(top10_cum_freq = sum(clone_freq_percent))

write.csv(cumulative_freq, "cumulative_top10_freq_TRB.csv", row.names = T)

# PCA from gene usages (unique clones)

library(tidyr)
library(ggplot2)

df_all <- read.csv("combined_df_allclones_TRB.csv")

V_usage <- df_all %>%
  filter(!is.na(V_gene)) %>%
  group_by(sample, V_gene) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

D_usage <- df_all %>%
  filter(!is.na(D_gene)) %>%
  group_by(sample, D_gene) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

J_usage <- df_all %>%
  filter(!is.na(J_gene)) %>%
  group_by(sample, J_gene) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

V_matrix <- V_usage %>%
  select(sample, V_gene, freq) %>%
  pivot_wider(names_from = V_gene, values_from = freq, values_fill = 0)

# Remove the sample column
V_pca_input <- as.data.frame(V_matrix)
rownames(V_pca_input) <- V_pca_input$sample
V_pca_input$sample <- NULL

#Hint: this code ensures pca only sees numeric values: df_numeric <- dplyr::select(df, where(is.numeric))

V_pca <- prcomp(V_pca_input, scale. = TRUE)

V_pca_df <- as.data.frame(V_pca$x)
V_pca_df$sample <- rownames(V_pca_df)


ggplot(V_pca_df, aes(PC1, PC2, label = sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.7) +
  theme_classic() +
  labs(title = "PCA of V gene usage")

# Adding categories
V_pca_df$sample

library(stringr)
library(RColorBrewer)
# Extracting the group names from the sample names, merging HD + Het


V_pca_df <- as.data.frame(V_pca$x)
V_pca_df$sample <- rownames(V_pca_df)

V_pca_df <- V_pca_df %>%
  mutate(group = case_when(
    grepl("^HD", sample) ~ "HD+Het",
    grepl("^Het", sample) ~ "HD+Het",
    TRUE ~ "pRD"
  ))

colors <- c("pRD" = "#E69F00", "HD+Het" = "#48D1CC")
plot_pca <- ggplot(V_pca_df, aes(PC1, PC2, color = group, fill = group)) +
  geom_point(size = 4) +
  stat_ellipse(type = "t", geom = "polygon", alpha = 0.2) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)+
  theme_classic()

ggsave("2D_PCA_TRBV_HD_Het_merged.png", plot_pca, bg = "white")

# disease categories
#Remove HD+Het
disease_df <- V_pca_df %>%
  dplyr::filter(!group %in% "HD+Het")

# Assign the categories
rownames(disease_df)

# Define which samples belong to each category
mild_samples     <- c( "P15_1index20_TRB_allclone","P12_1index16_TRB_allclone", "P2_T1_1index13_TRB_allclone", "P13_1index19_TRB_allclone")
moderate_samples <- c("P11_1index18_TRB_allclone", "P4_T1_2index9_TRB_allclone", "P10_1index15_TRB_allclone")
severe_samples   <- c( "P1_T2_2index2_TRB_allclone", "P2_T2_2index3_TRB_allclone","P6_T1_2index1_TRB_allclone", "P1_T1_1index12_TRB_allclone")

# Assign categories to the disease dataframe
disease_df <- disease_df %>%
  mutate(group = case_when(
    sample %in% mild_samples     ~ "mild",
    sample %in% moderate_samples ~ "moderate",
    sample %in% severe_samples   ~ "severe"
  ))

# Check assignment
table(disease_df$group)


ggplot(disease_df, aes(x = PC1, y = PC2, color = group, fill = group)) +
  geom_point(size = 4) +
  stat_ellipse(geom = "polygon", alpha = 0.2) +
  scale_color_manual(values = c(
    mild = "#FFB6C1",     # pink
    moderate = "#FFA500", # orange
    severe = "#FF0000"    # red
  )) +
  scale_fill_manual(values = c(
    mild = "#FFB6C1",     # pink
    moderate = "#FFA500", # orange
    severe = "#FF0000"    # red
  )) +
  theme_minimal()

# Top 10 genes for PC1
top10_PC1 <- sort(abs(V_pca$rotation[, "PC1"]), decreasing = TRUE)[1:10]

# Top 10 genes for PC2
top10_PC2 <- sort(abs(V_pca$rotation[, "PC2"]), decreasing = TRUE)[1:10]

top10_PC1
top10_PC2

write.csv(top10_PC1, "top10_PC1_TRB.csv", row.names = T)
write.csv(top10_PC2, "top10_PC2_TRB.csv", row.names = T)
write.csv(V_pca_df, "PCA_df_unique_TRB_clones.csv", row.names = T)

# Doing the same with weighted clone size
V_usage_clone <- df_all %>%
  filter(!is.na(V_gene)) %>%
  group_by(sample, V_gene) %>%
  summarise(count = sum(duplicate_count), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

write.csv(V_usage_clone, "TRBV_usage_weighted.csv", row.names = T)


V_matrix_cl <- V_usage_clone %>%
  select(sample, V_gene, freq) %>%
  pivot_wider(names_from = V_gene, values_from = freq, values_fill = 0)

# Remove the sample column
V_pca_input_cl <- as.data.frame(V_matrix_cl)
rownames(V_pca_input_cl) <- V_pca_input_cl$sample
V_pca_input_cl$sample <- NULL

#Hint: this code ensures pca only sees numeric values: df_numeric <- dplyr::select(df, where(is.numeric))

V_pca_cl <- prcomp(V_pca_input_cl, scale. = TRUE)

V_pca_df_cl <- as.data.frame(V_pca_cl$x)
V_pca_df_cl$sample <- rownames(V_pca_df_cl)

ggplot(V_pca_df_cl, aes(PC1, PC2, label = sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.7) +
  theme_classic() 

# Adding categories
V_pca_df_cl$sample

# Extracting the group names from the sample names, merging HD + Het

V_pca_df_cl$sample <- rownames(V_pca_df_cl)

V_pca_df_cl <- V_pca_df_cl %>%
  mutate(group = case_when(
    grepl("^HD", sample) ~ "HD+Het",
    grepl("^Het", sample) ~ "HD+Het",
    TRUE ~ "pRD"
  ))


colors <- c("pRD" = "#E69F00", "HD+Het" = "#48D1CC")
plot_pca <- ggplot(V_pca_df_cl, aes(PC1, PC2, color = group, fill = group)) +
  geom_point(size = 4) +
  stat_ellipse(type = "t", geom = "polygon", alpha = 0.2) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors)+
  theme_classic()
plot_pca
ggsave("2D_PCA_TRBV_HD_Het_merged_weighted_to_clonesize.png", plot_pca, bg = "white")

# disease categories
#Remove HD+Het
disease_df_cl <- V_pca_df_cl %>%
  dplyr::filter(!group %in% "HD+Het")

# Assign the categories
rownames(disease_df_cl)

# Define which samples belong to each category
mild_samples     <- c( "P15_1index20_TRB_allclone","P12_1index16_TRB_allclone", "P2_T1_1index13_TRB_allclone", "P13_1index19_TRB_allclone")
moderate_samples <- c("P11_1index18_TRB_allclone", "P4_T1_2index9_TRB_allclone", "P10_1index15_TRB_allclone")
severe_samples   <- c( "P1_T2_2index2_TRB_allclone", "P2_T2_2index3_TRB_allclone","P6_T1_2index1_TRB_allclone", "P1_T1_1index12_TRB_allclone")

# Assign categories to the disease dataframe
disease_df_cl <- disease_df_cl %>%
  mutate(group = case_when(
    sample %in% mild_samples     ~ "mild",
    sample %in% moderate_samples ~ "moderate",
    sample %in% severe_samples   ~ "severe"
  ))

# Check assignment
table(disease_df$group)


ggplot(disease_df_cl, aes(x = PC1, y = PC2, color = group, fill = group)) +
  geom_point(size = 4) +
  stat_ellipse(geom = "polygon", alpha = 0.2) +
  scale_color_manual(values = c(
    mild = "#FFB6C1",     # pink
    moderate = "#FFA500", # orange
    severe = "#FF0000"    # red
  )) +
  scale_fill_manual(values = c(
    mild = "#FFB6C1",     # pink
    moderate = "#FFA500", # orange
    severe = "#FF0000"    # red
  )) +
  theme_minimal()

# Top 10 genes for PC1
top10_PC1_cl <- sort(abs(V_pca_cl$rotation[, "PC1"]), decreasing = TRUE)[1:10]

# Top 10 genes for PC2
top10_PC2_cl <- sort(abs(V_pca_cl$rotation[, "PC2"]), decreasing = TRUE)[1:10]

top10_PC1_cl
top10_PC2_cl

write.csv(top10_PC1_cl, "top10_PC1_TRB_weigthed_to_clonesize.csv", row.names = T)
write.csv(top10_PC2_cl, "top10_PC2_TRB_weighted_to_clonesize.csv", row.names = T)
write.csv(V_pca_df_cl, "PCA_df_weighted_TRB_clones.csv", row.names = T)
write.csv(V_usage, "V_usage_unique_clones_TRB.csv", row.names = T)
write.csv(V_usage_clone, "V_usage_weighted_TRB.csv", row.names = T)
write.csv(D_usage, "D_usage_unique_clones_TRB.csv", row.names = T)
write.csv(J_usage, "J_usage_unique_clones_TRB.csv", row.names = T)

# Top 10 PC1 loadings from unique clones

genes_of_interest <- c("TRBV7-9", "TRBV6-6", "TRBV7-2", "TRBV29-1", "TRBV11-2",
                       "TRBV9", "TRBV4-1", "TRBV25-1", "TRBV7-3", "TRBV12-3")

top10_PC1_unique <- V_usage %>%
  filter(V_gene %in% genes_of_interest) %>%
  select(sample, V_gene, freq)

write.csv(top10_PC1_unique, "top10_PC1_loadings_per_sample_unique.csv", row.names = T)

# PCA only for the disease patients from the gene usage values

V_usage_patients <- V_usage %>%
  filter(!grepl("^HD|^Het", sample))

# Make it tidy
V_matrix_patients <- V_usage_patients %>%
  select(sample, V_gene, freq) %>%
  pivot_wider(names_from = V_gene, values_from = freq, values_fill = 0)

# Remove the sample column
sample_names <- V_matrix_patients$sample
V_matrix_patients <- V_matrix_patients %>% select(-sample)

# Run PCA
V_pca_patients <- prcomp(V_matrix_patients, scale. = TRUE)

V_pca_df_patients <- as.data.frame(V_pca_patients$x)
V_pca_df_patients$sample <- sample_names


ggplot(V_pca_df_patients, aes(PC1, PC2)) +
  geom_point(size = 4) +
  theme_minimal()

# Adding the groups
V_pca_df_patients <- V_pca_df_patients %>%
  mutate(group = case_when(
    sample %in% mild_samples ~ "mild",
    sample %in% moderate_samples ~ "moderate",
    sample %in% severe_samples ~ "severe"
  ))

table(V_pca_df_patients$group)


p <- ggplot(V_pca_df_patients, aes(PC1, PC2, color = group, fill = group)) +
  geom_point(size = 4) +
  stat_ellipse(geom = "polygon", alpha = 0.2) +
  scale_color_manual(values = c(
    mild = "#FFB6C1",
    moderate = "#FFA500",
    severe = "#FF0000"
  )) +
  scale_fill_manual(values = c(
    mild = "#FFB6C1",
    moderate = "#FFA500",
    severe = "#FF0000"
  )) +
  theme_minimal()

ggsave("categories_PCA_unique_TRBV.png", p, bg = "white")

# Same with weighted clones

V_usage_patients_cl <- V_usage_clone %>%
  filter(!grepl("^HD|^Het", sample))

# Make it tidy
V_matrix_patients_cl <- V_usage_patients_cl %>%
  select(sample, V_gene, freq) %>%
  pivot_wider(names_from = V_gene, values_from = freq, values_fill = 0)

# Remove the sample column
sample_names_cl <- V_matrix_patients_cl$sample
V_matrix_patients_cl <- V_matrix_patients_cl %>% select(-sample)

# Run PCA
V_pca_patients_cl <- prcomp(V_matrix_patients_cl, scale. = TRUE)

V_pca_df_patients_cl <- as.data.frame(V_pca_patients_cl$x)
V_pca_df_patients_cl$sample <- sample_names


ggplot(V_pca_df_patients_cl, aes(PC1, PC2)) +
  geom_point(size = 4) +
  theme_minimal()

# Adding the groups
V_pca_df_patients_cl <- V_pca_df_patients_cl %>%
  mutate(group = case_when(
    sample %in% mild_samples ~ "mild",
    sample %in% moderate_samples ~ "moderate",
    sample %in% severe_samples ~ "severe"
  ))

table(V_pca_df_patients_cl$group)


p2 <- ggplot(V_pca_df_patients_cl, aes(PC1, PC2, color = group, fill = group)) +
  geom_point(size = 4) +
  stat_ellipse(geom = "polygon", alpha = 0.2) +
  scale_color_manual(values = c(
    mild = "#FFB6C1",
    moderate = "#FFA500",
    severe = "#FF0000"
  )) +
  scale_fill_manual(values = c(
    mild = "#FFB6C1",
    moderate = "#FFA500",
    severe = "#FF0000"
  )) +
  theme_minimal()

ggsave("categories_PCA_weighted_TRBV.png", p2, bg = "white")

