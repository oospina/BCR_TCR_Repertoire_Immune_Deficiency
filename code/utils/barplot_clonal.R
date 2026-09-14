##
# barplot_clonal.R
# 
# By Oscar Ospina
# Created: Feb 03, 2026
#

library('tidyverse')
library('alakazam')

# size of random downsampling of data
sampsize = 2500

# input variables (mostly generated from parsed input filename)
args <- commandArgs(TRUE)  # read arguments into "args" vector
ccfile <- args[1]          # color-code file - 1st argument
#ccfile = './data/dest/IGH_v_gene_colorcodes.tsv'
idfile = args[2] # Key with sample index to sample name - 2nd argument
#idfile = './data/raw_data/tcrbcr_index_sample_key.txt' 
tmfile <- args[3]          # barplot file - 3rd argument
dbfile <- args[4:length(args)]          # List of AIRR files - 4th argument
#dbfile = list.files('./data/clone/', pattern='_IGH_allclone.tsv', full.names=TRUE)

# Get index-patient key
pat = read.table(idfile)

# Read files
airr_df = lapply(dbfile, function(i){
  sample_idx = str_extract(i, '[0-9]index[0-9]+')
  pat_id = pat[[2]][pat[[1]] == sample_idx]
  
  df_tmp = read_tsv(i) %>%
    mutate(sample_id=sample_idx,
           pat_id=pat_id,
           comp=str_extract(i, 'IGH|TRB'))
}) %>% bind_rows()

# read uniform color-code table
cc_df = read_tsv(ccfile)

# Get gene names
airr_df$v_gene = getGene(airr_df$v_call)
airr_df = airr_df %>% left_join(., cc_df, by=c("v_gene"="gene"))

# generate a clone size distribution table with V-genes and hues
clonesizes = airr_df %>% 
  select(c('pat_id', 'clone_id', 'duplicate_count', 'v_gene')) %>%
  group_by(pat_id, v_gene) %>%
  summarize(clone_ct=n()) %>%
  mutate(total_clones=sum(clone_ct),
         clone_perc=clone_ct/total_clones, .groups='drop') %>%
  arrange(desc(total_clones)) %>%
  mutate(v_gene=factor(v_gene, levels=cc_df[['gene']][ cc_df[['gene']] %in% airr_df[['v_gene']] ]))

# Make barplot
bp1 = ggplot(clonesizes) +
  geom_bar(aes(x=pat_id, y=clone_ct, fill=v_gene), stat='identity') +
  labs(title=paste0('Number of clones\n', unique(airr_df$comp), ' - Top ', sampsize, ' most abundant clones'),
       x=NULL, y='Number of clones', x=NULL) +
  scale_fill_manual(values=setNames(cc_df[['color']], nm=cc_df[['gene']])) +
  theme_minimal() +
  theme(panel.border=element_rect(fill=NA, color='black'),
        axis.text.x=element_text(angle=70, vjust=1, hjust=1))

bp2 = ggplot(clonesizes) +
  geom_bar(aes(x=pat_id, y=clone_perc, fill=v_gene), stat='identity') +
  labs(title=paste0('Percentage of clones within samples\n', unique(airr_df$comp), ' - Top ', sampsize, ' most abundant clones'), 
       y='Percentage of clones', x=NULL) +
  scale_fill_manual(values=setNames(cc_df[['color']], nm=cc_df[['gene']])) +
  theme_minimal() +
  theme(panel.border=element_rect(fill=NA, color='black'),
        axis.text.x=element_text(angle=70, vjust=1, hjust=1))

pdf(tmfile, width=16, height=8)
print(ggpubr::ggarrange(bp1, bp2, ncol=2, common.legend=TRUE, legend='right'))
dev.off()

