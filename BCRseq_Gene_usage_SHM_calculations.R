

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


clone_ratio_2 <- df_all %>%
  group_by(sample) %>%
  summarise(
    unique_clones = n_distinct(clone_id),
    total_sequences = sum(duplicate_count),
    ratio_unique_total = unique_clones / total_sequences
  )

clone_ratio_2

write.csv(clone_ratio_2, "clone_ratio.csv", row.names = T)

# J5 and J6 usage

df_all <- df_all %>%
  mutate(Jgene = str_extract(j_call, "IGHJ[0-9]+"))

# weighted to the clone counts

J_usage <- df_all %>%
  group_by(sample) %>%
  summarise(
    J5_count = sum(duplicate_count[Jgene == "IGHJ5"]),
    J6_count = sum(duplicate_count[Jgene == "IGHJ6"]),
    total_J = sum(duplicate_count),
    J5_freq = J5_count / total_J,
    J6_freq = J6_count / total_J
  )

write.csv(J_usage, "J_gene_usage.csv", row.names = T)

# Heavy gene usages (per sample)

df_all <- df_all %>%
  mutate(
    Vgene = str_extract(v_call, "IGHV[0-9]+-[0-9]+"),
    Dgene = str_extract(d_call, "IGHD[0-9]+-[0-9]+"),
    Jgene = str_extract(j_call, "IGHJ[0-9]+")
  )

V_usage <- df_all %>%
  filter(!is.na(Vgene)) %>%
  group_by(sample, Vgene) %>%
  summarise(count = sum(duplicate_count), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

D_usage <- df_all %>%
  filter(!is.na(Dgene)) %>%
  group_by(sample, Dgene) %>%
  summarise(count = sum(duplicate_count), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

J_usage <- df_all %>%
  filter(!is.na(Jgene)) %>%
  group_by(sample, Jgene) %>%
  summarise(count = sum(duplicate_count), .groups = "drop") %>%
  group_by(sample) %>%
  mutate(freq = count / sum(count)) %>%
  ungroup()

# Combined

heavy_usage <- bind_rows(
  V_usage %>% mutate(chain="V", gene = Vgene) %>% select(sample, chain, gene, count, freq),
  D_usage %>% mutate(chain="D", gene = Dgene) %>% select(sample, chain, gene, count, freq),
  J_usage %>% mutate(chain="J", gene = Jgene) %>% select(sample, chain, gene, count, freq)
)
write.csv(heavy_usage, "heavy_chain_usages_freq.csv", row.names = T)
write.csv(df_all, "combined_df_allclones.csv", row.names = T)
write.csv(V_usage, "IGHV_usage,csv", row.names = T)
write.csv(D_usage, "IGHD_usage,csv", row.names = T)
write.csv(J_usage, "IGHJ_usage,csv", row.names = T)


# making a heatmap 

V_order <- c(
  "IGHV7-81","IGHV3-74","IGHV3-73","IGHV3-72","IGHV2-70","IGHV1-69",
  "IGHV3-66","IGHV3-64","IGHV4-61","IGHV4-59","IGHV1-58","IGHV3-53",
  "IGHV5-51","IGHV3-49","IGHV3-48","IGHV1-46","IGHV1-45","IGHV3-43",
  "IGHV4-39","IGHV3-38","IGHV3-35","IGHV4-34","IGHV3-33","IGHV4-31",
  "IGHV3-30","IGHV4-28","IGHV2-26","IGHV1-24","IGHV3-23","IGHV3-21",
  "IGHV3-20","IGHV1-18","IGHV3-16","IGHV3-15","IGHV3-13","IGHV3-11",
  "IGHV3-9","IGHV1-8","IGHV3-7","IGHV2-5","IGHV4-4","IGHV1-3",
  "IGHV1-2","IGHV6-1"
)


V_usage <- V_usage %>%
  complete(sample, Vgene = V_order, fill = list(count = 0, freq = 0))

V_usage$Vgene <- factor(V_usage$Vgene, levels = V_order)

ggplot(V_usage, aes(x = Vgene, y = sample, fill = freq)) +
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

# there are 3 non-used segments
unused_V <- V_usage %>%
  group_by(Vgene) %>%
  summarise(total = sum(freq, na.rm = TRUE)) %>%
  filter(total == 0) %>%
  pull(Vgene)

unused_V
# Remove them

V_usage_filtered <- V_usage %>%
  filter(!Vgene %in% unused_V)

# update the order vector

V_order_filtered <- V_order[!V_order %in% unused_V]

V_usage_filtered <- V_usage_filtered[!is.na(V_usage_filtered$Vgene), ]
V_usage_filtered$Vgene <- factor(
  V_usage_filtered$Vgene,
  levels = V_order_filtered
)

unique(V_usage_filtered$Vgene)
range(V_usage_filtered$freq, na.rm = TRUE)


ggplot(V_usage_filtered, aes(x = Vgene, y = sample, fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    na.value = "white",
    name = "IGHV usage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Reordering and renaming the samples
V_usage_filtered$sample <- as.character(V_usage_filtered$sample)
unique(V_usage_filtered$sample)


V_usage_filtered$sample <- dplyr::recode(
  V_usage_filtered$sample,
  "P6_T1_2index1_IGH_allclone" = "P6",
  "P2_T2_2index3_IGH_allclone" = "P2_T2", 
  "P1_T1_1index12_IGH_allclone" = "P1_T1", 
  "P1_T2_2index2_IGH_allclone" = "P1_T2", 
  "P10_1index15_IGH_allclone" = "P10", 
  "P4_T1_2index9_IGH_allclone" = "P4", 
  "P11_1index18_IGH_allclone" = "P11", 
  "P15_1index20_IGH_allclone" = "P15", 
  "P12_1index16_IGH_allclone" = "P12", 
  "P13_1index19_IGH_allclone" = "P13", 
  "P2_T1_1index13_IGH_allclone" = "P2_T1", 
  "Het1_2index15_IGH_allclone" = "Het1", 
  "Het2_2index16_IGH_allclone" = "Het2", 
  "Het4_2index18_IGH_allclone" = "Het4", 
  "Het5_2index8_IGH_allclone" = "Het5", 
  "HD1_1index27_IGH_allclone" = "HD1", 
  "HD9_1index25_IGH_allclone" = "HD9", 
  "HD12_1index23_IGH_allclone" = "HD12"
)

sample_order <- c("P2_T2", "P6", "P1_T1", "P1_T2", 
                  "P10", "P11", "P4", "P15", "P13", 
                  "P12", "P2_T1", "Het1", "Het2", "Het4",
                  "Het5", "HD1", "HD9", "HD12")



V_usage_filtered$sample <- factor(
  V_usage_filtered$sample,
  levels = sample_order)


write.csv(V_usage_filtered, "V_usage_filtered_IGH.csv", row.names = T)


# Then plot as usual
p <- ggplot(V_usage_filtered, aes(x = Vgene, y = forcats::fct_rev(sample), fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#ffffcc", "#a1dab4", "#41b6c4", "#225ea8"),
    values = scales::rescale(c(0, 0.001, 0.01, 0.05, max(V_usage_filtered$freq, na.rm = TRUE))),
    na.value = "white",
    name = "IGHV usage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_fixed(ratio = 1)

p
ggsave("V_gene_usage_heatmap_IGH.png", p, bg = "white")

# D gene segments


D_order <- c("IGHD1-1", "IGHD2-2", "IGHD3-3", "IGHD4-4", "IGHD5-5", "IGHD6-6", "IGHD1-7", "IGHD2-8", "IGHD3-9", "IGHD3-10", 
             "IGHD4-11", "IGHD5-12", "IGHD6-13", "IGHD1-14", "IGHD2-15", "IGHD3-16", "IGHD4-17", "IGHD5-18", "IGHD6-19", "IGHD1-20", 
             "IGHD2-21", "IGHD3-22", "IGHD4-23", "IGHD5-24", "IGHD6-25", "IGHD1-26", "IGHD7-27")

D_usage <- D_usage %>%
  complete(sample, Dgene = D_order, fill = list(count = 0, freq = 0))

D_usage$Dgene <- factor(D_usage$Dgene, levels = D_order)

# Making zeros to NA
# D_usage$freq[D_usage$freq == 0] <- NA

ggplot(D_usage, aes(x = Dgene, y = sample, fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(
    colours = c("white", "blue", "yellow"),
    values = scales::rescale(c(0, 0.01, max(V_usage$freq, na.rm = TRUE))),
    na.value = "white",
    name = "Frequency"
  ) +
  labs(x = "D gene (distal → proximal)", y = "Sample") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# there are 2 non-used segments
unused_D <- D_usage %>%
  group_by(Dgene) %>%
  summarise(total = sum(freq, na.rm = TRUE)) %>%
  filter(total == 0) %>%
  pull(Dgene)

unused_D
# Remove them

D_usage_filtered <- D_usage %>%
  filter(!Dgene %in% unused_D)

# update the order vector

D_order_filtered <- D_order[!D_order %in% unused_D]

D_usage_filtered <- D_usage_filtered[!is.na(D_usage_filtered$Dgene), ]
D_usage_filtered$Dgene <- factor(
  D_usage_filtered$Dgene,
  levels = D_order_filtered
)

unique(D_usage_filtered$Dgene)
range(D_usage_filtered$freq, na.rm = TRUE)


ggplot(D_usage_filtered, aes(x = Dgene, y = sample, fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    na.value = "white",
    name = "IGHD usage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Reordering and renaming the samples
D_usage_filtered$sample <- as.character(D_usage_filtered$sample)
unique(D_usage_filtered$sample)


D_usage_filtered$sample <- dplyr::recode(
  D_usage_filtered$sample,
  "P6_T1_2index1_IGH_allclone" = "P6",
  "P2_T2_2index3_IGH_allclone" = "P2_T2", 
  "P1_T1_1index12_IGH_allclone" = "P1_T1", 
  "P1_T2_2index2_IGH_allclone" = "P1_T2", 
  "P10_1index15_IGH_allclone" = "P10", 
  "P4_T1_2index9_IGH_allclone" = "P4", 
  "P11_1index18_IGH_allclone" = "P11", 
  "P15_1index20_IGH_allclone" = "P15", 
  "P12_1index16_IGH_allclone" = "P12", 
  "P13_1index19_IGH_allclone" = "P13", 
  "P2_T1_1index13_IGH_allclone" = "P2_T1", 
  "Het1_2index15_IGH_allclone" = "Het1", 
  "Het2_2index16_IGH_allclone" = "Het2", 
  "Het4_2index18_IGH_allclone" = "Het4", 
  "Het5_2index8_IGH_allclone" = "Het5", 
  "HD1_1index27_IGH_allclone" = "HD1", 
  "HD9_1index25_IGH_allclone" = "HD9", 
  "HD12_1index23_IGH_allclone" = "HD12"
)

sample_order <- c("P2_T2", "P6", "P1_T1", "P1_T2", 
                  "P10", "P11", "P4", "P15", "P13", 
                  "P12", "P2_T1", "Het1", "Het2", "Het4",
                  "Het5", "HD1", "HD9", "HD12")



D_usage_filtered$sample <- factor(
  D_usage_filtered$sample,
  levels = sample_order)


unique(D_usage_filtered$sample)

write.csv(D_usage_filtered, "D_usage_filtered_IGH.csv", row.names = T)

p_d <- ggplot(D_usage_filtered, aes(x = Dgene, y = forcats::fct_rev(sample), fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#ffffcc", "#a1dab4", "#41b6c4", "#225ea8"),
    values = scales::rescale(c(0, 0.001, 0.01, 0.05, max(D_usage_filtered$freq, na.rm = TRUE))),
    na.value = "white",
    name = "IGHD usage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14)) + 
  coord_fixed(ratio = 1)


p_d
ggsave("D_gene_usage_IGH.png", p_d, bg = "white")


# J segments

J_order <- c("IGHJ1", "IGHJ2", "IGHJ3", "IGHJ4", "IGHJ5", "IGHJ6")

J_usage <- J_usage %>%
  complete(sample, Jgene = J_order, fill = list(count = 0, freq = 0))

J_usage$Jgene <- factor(J_usage$Jgene, levels = J_order)

ggplot(J_usage, aes(x = Jgene, y = sample, fill = freq)) +
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


# no non-used J

unique(J_usage$Jgene)
range(J_usage$freq, na.rm = TRUE)

ggplot(J_usage, aes(x = Jgene, y = sample, fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(
    colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    na.value = "white",
    name = "IGHJ usage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Renaming the samples
J_usage$sample <- as.character(J_usage$sample)

J_usage$sample <- dplyr::recode(
  J_usage$sample,
  "P6_T1_2index1_IGH_allclone" = "P6",
  "P2_T2_2index3_IGH_allclone" = "P2_T2", 
  "P1_T1_1index12_IGH_allclone" = "P1_T1", 
  "P1_T2_2index2_IGH_allclone" = "P1_T2", 
  "P10_1index15_IGH_allclone" = "P10", 
  "P4_T1_2index9_IGH_allclone" = "P4", 
  "P11_1index18_IGH_allclone" = "P11", 
  "P15_1index20_IGH_allclone" = "P15", 
  "P12_1index16_IGH_allclone" = "P12", 
  "P13_1index19_IGH_allclone" = "P13", 
  "P2_T1_1index13_IGH_allclone" = "P2_T1", 
  "Het1_2index15_IGH_allclone" = "Het1", 
  "Het2_2index16_IGH_allclone" = "Het2", 
  "Het4_2index18_IGH_allclone" = "Het4", 
  "Het5_2index8_IGH_allclone" = "Het5", 
  "HD1_1index27_IGH_allclone" = "HD1", 
  "HD9_1index25_IGH_allclone" = "HD9", 
  "HD12_1index23_IGH_allclone" = "HD12"
)

sample_order <- c("P2_T2", "P6", "P1_T1", "P1_T2", 
                  "P10", "P11", "P4", "P15", "P13", 
                  "P12", "P2_T1", "Het1", "Het2", "Het4",
                  "Het5", "HD1", "HD9", "HD12")

J_usage$sample <- factor(
  J_usage$sample,
  levels = sample_order)

unique(J_usage$sample)

write.csv(J_usage, "J_usage_filtered_IGH.csv", row.names = T)

p_j <- ggplot(J_usage, aes(x = Jgene, y = forcats::fct_rev(sample), fill = freq)) +
  geom_tile(color = "grey90") +
  scale_fill_gradientn(
    colours = c("#FFFFFF", "#ffffcc", "#a1dab4", "#41b6c4", "#225ea8"),
    values = scales::rescale(c(0, 0.001, 0.01, 0.05, max(J_usage$freq, na.rm = TRUE))),
    na.value = "white",
    name = "IGHJ usage"
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
ggsave("J_gene_usage_IGH.png", p_j, bg = "white")


# IGHV4-34 usage

gene_of_interest <- "IGHV4-34"

V_single <- V_usage_filtered %>%
  filter(Vgene == gene_of_interest) %>%
  select(sample, freq)

write.csv(V_single, "IGHV4-34_usage.csv", row.names = T)

# SHM calculations (rsmut Excel files)

library(readxl)

# list all Excel files
files <- list.files(pattern = "\\.xlsx$", full.names = TRUE)

# read and combine
combined <- lapply(files, read_excel) %>%
  bind_rows(.id = "sample_file")


combined <- lapply(files, function(x) {
  df <- read_excel(x)
  df$sample <- tools::file_path_sans_ext(basename(x))
  return(df)
}) %>%
  bind_rows()

write.csv(combined, "combine_mu_calculations.csv", row.names = T)

# Renaming the samples
unique(combined$sample)
class(combined$sample)

combined$sample <- recode(combined$sample,
                          "1index12_IGH_rsmut" = "P1_T1", "1index13_IGH_rsmut" = "P2_T1", "1index15_IGH_rsmut" = "P10", 
                          "1index16_IGH_rsmut" = "P12", "1index18_IGH_rsmut" = "P11", "1index19_IGH_rsmut" = "P13", 
                          "1index20_IGH_rsmut" = "P15", "1index23_IGH_rsmut" = "HD12", "1index25_IGH_rsmut" = "HD9", 
                          "1index27_IGH_rsmut" = "HD1", "2index1_IGH_rsmut" = "P6", "2index15_IGH_rsmut" = "Het1", 
                          "2index16_IGH_rsmut" = "Het2", "2index18_IGH_rsmut" = "Het4", "2index2_IGH_rsmut" = "P1_T2", 
                          "2index3_IGH_rsmut" = "P2_T2", "2index8_IGH_rsmut" = "Het5", "2index9_IGH_rsmut" = "P4")

table(combined$sample)   

# Saving again
write.csv(combined, "combine_mu_calculations.csv", row.names = T)

# Define the order

combined$sample <- factor(combined$sample,
                          levels = c(
                            "P6", "P2_T2",
                            "P1_T1", "P1_T2",
                            "P10", "P4", "P11", "P15", "P12", "P13", "P2_T1",
                            "Het1", "Het2", "Het4", "Het5",
                            "HD1", "HD9", "HD12"
                          )
)

levels(combined$sample)

# Calculate the V region length for the frequency

combined$V_length <- combined$v_sequence_end - combined$v_sequence_start + 1                

combined$mutation_freq <- combined$mu_count / combined$V_length

summary(combined$mutation_freq)

# Mean SHM per sample


sample_SHM <- combined %>%
  group_by(sample) %>%
  summarise(mean_SHM = mean(mutation_freq, na.rm = TRUE))

write.csv(sample_SHM, "SHM_per_sample.csv", row.names = T)

colnames(combined)
# SHM frequency by clones

clone_SHM <- combined %>%
  group_by(sample, clone_id) %>%
  summarise(
    mean_SHM = mean(mutation_freq, na.rm = TRUE),   # SHM per clone
    clone_size = sum(duplicate_count),             # total sequences per clone
    n_sequences = n(),                               # number of rows collapsed
    .groups = "drop"
  )

weighted_SHM <- clone_SHM %>%
  group_by(sample) %>%
  summarise(
    weighted_SHM = sum(mean_SHM * clone_size) / sum(clone_size),
    n_clones = n()
  )
write.csv(weighted_SHM, "weighted_SHM_by_clone_size.csv", row.names = T)


weighted_SHM$group <- case_when(
  grepl("^P", weighted_SHM$sample) ~ "pRD",
  grepl("^Het", weighted_SHM$sample) ~ "Het",
  grepl("^HD", weighted_SHM$sample) ~ "HD"
)

ggplot(weighted_SHM, aes(x = group, y = weighted_SHM)) +
  geom_bar(stat = "identity", fill = "lightblue") +
  theme_minimal() +
  labs(y = "Weighted SHM per sample", x = "Group") 

ggplot(weighted_SHM, aes(x = sample, y = weighted_SHM, fill = group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(y = "Weighted SHM per sample", x = "Sample")

group_SHM <- weighted_SHM %>%
  group_by(group) %>%
  summarise(mean_weighted_SHM = mean(weighted_SHM))

ggplot(group_SHM, aes(x = group, y = mean_weighted_SHM, fill = group)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(y = "Mean weighted SHM per group", x = "Group")


# Top 10 clones per sample

clone_freq <- combined %>%
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


write.csv(top10_clones, "top10_clones_frequency_per_patient.csv", row.names = T)

# cumulative frequency

cumulative_freq <- top10_clones %>%
  group_by(sample) %>%
  summarise(top10_cum_freq = sum(clone_freq_percent))

write.csv(cumulative_freq, "cumulative_top10_freq.csv", row.names = T)


