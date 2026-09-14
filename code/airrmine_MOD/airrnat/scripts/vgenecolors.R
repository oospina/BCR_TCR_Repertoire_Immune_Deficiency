# -----------------------------------------------------------------------------
# vgenecolors.R
# version: 1.0 (2020-09-21)
# author : Peter Blazso
#
# Collects V-genes from specified repertoires and assigns unified color codes
# -----------------------------------------------------------------------------

library( openxlsx )
library( dplyr )
library( readr )
library( colorspace )

# arguments
args    <- commandArgs(TRUE)
#args = list.files('./data/dest/', pattern='_TRB_genes.xlsx', full.names=TRUE)
outfile <- args[1]  # 1st argument is destination datafile

# all other arguments are input gene usage files
# read in and merge all data into one tibble
df <- tibble()
for( arg in args[2:length(args)] )
{
  # read sheet #1 (with V-gene usage list) from next XLSX file
  df_tmp = read.xlsx( arg, sheet=1 ) # To avoid merge of files without sequences (e.g., sample 1index14)
  if(all(unique(df_tmp[['seq_count']]) != "") ){
  	df <- bind_rows( df, df_tmp )
  } else{
    cat(paste0('NO SEQUENCES AVAILABLE IN THIS SAMPLE ( ', stringr::str_extract(arg, '[0-9]index[0-9]+_[IGHTRB]{3}'), ' )\n'))
  }
  # df <- bind_rows( df, read.xlsx( arg, sheet=1 ) )
}

# make a unique list of V-genes
df <- df %>% select( gene ) %>% group_by( gene ) %>% summarize() %>% ungroup()

# generate additional id to sort the genes correctly
## MODIFIED AS SOME TRB GENES DONT HAVE NUMERALS
# df$id <- sprintf( "%s%03d",
#                   unlist( lapply(strsplit(df$gene,"-"),getElement,1)),
#                   as.numeric( 
#                     unlist( lapply(strsplit(df$gene,"-"),getElement,2) ))
# )
parsed_tmp = strsplit(df$gene,"-")
parsed_tmp = lapply(parsed_tmp, function(i){if(length(i) == 1){i=c(i, ' ')} else{i}})
gene_tmp = unlist( lapply(parsed_tmp, getElement, 1) )
num_tmp = as.numeric( unlist( lapply(parsed_tmp, getElement, 2) ) )
df$id <- sprintf( "%s%03d", gene_tmp, num_tmp) %>% gsub(' NA', '', .)

# sort gene names
df <- df %>% arrange( id ) %>% select( -id )

# assign hue and sample color values to this sorted list
df$hue   <- round(seq(1, 360, length.out = nrow(df)))
df$color <- hex(HSV(df$hue, 0.8, 1))

# write out color-coded list
write_tsv( df, outfile )
