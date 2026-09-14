# -----------------------------------------------------------------------------
# treemap.R
# version: 1.0 (2020-09-28)
# author : Peter Blazso, Krisztian Csomos
#
# Generates uniformly color-coded treemap of V-gene usage/abundance
# -----------------------------------------------------------------------------

library( readr )
library( dplyr )
#library( treemap )
library( treemapify )
library( colorspace )
library( alakazam )

# size of random downsampling of data
sampsize = 2500

# input variables (mostly generated from parsed input filename)
args <- commandArgs(TRUE)  # read arguments into "args" vector
ccfile <- args[1]          # color-code file - 1st argument
#ccfile = './data/dest/IGH_v_gene_colorcodes.tsv'
dbfile <- args[2]          # AIRR data file  - 2nd argument
#dbfile = './data/coll/2index12_IGH_collapsed.tsv'
idfile = args[3] # Key with sample index to sample name - 3rd argumen
#idfile = './data/raw_data/tcrbcr_index_sample_key.txt' 
tmfile <- args[4]          # treemap file    - 4rd argument

# dissect input database filename
dbfile_parts <- strsplit(basename(dbfile), "_")[[1]]
subj <- dbfile_parts[1]    # subject ID
comp <- dbfile_parts[2]    # B-cell compartment

# read uniform color-code table
cc_df   <- read_tsv( ccfile )

# read AIRR data
airr_df <- read_tsv( dbfile )

if(nrow(airr_df) > 0){ # In case no sequences are present for a sample (e.g., 2index25 TRB)
  # Get index-patient key
  pat = read.table(idfile)
  pat = pat[[2]][pat[[1]] == subj]
  
  # downsample AIRR data if original sample size is greater than limit
  if( nrow (airr_df) > sampsize )
    airr_df <- airr_df %>% sample_n( sampsize )
  
  # assign uniform hue values to different V-genes
  airr_df$v_gene <- getGene( airr_df$v_call )
  airr_df   <- airr_df %>% left_join( cc_df, by=c("v_gene"="gene") )
  
  # mark distinct sequences in V-gene groups by (brightness) value levels
  vals <-  seq( 1, 0.2, by=-0.01 )
  vals = vals[round(seq(1, length(vals), length.out=length(airr_df$hue)))] # Select equally spaced hues because some samples have fewer sequences
  airr_df$color <- hex(HSV(airr_df$hue, 0.8, vals))
  
  # Plot treemap
  #tiff( filename=tmfile, width=3600, height=2880, res=600, compression="zip" )
  #pdf( filename=tmfile, width=10, height=10)
  #pdf(tmfile, width=10, height=10)
  # treemap( airr_df, index=c("v_call","sequence_id"),
  #          vSize="duplicate_count", vColor="color",
  #          type="color",
  #          #title="", # title=paste( "AIRR TreeMap","-",paste0( subj,"_",comp ) ),
  #          title=paste( "V-gene abundance\n", paste0( comp, ' - ', subj, "; Patient: ", pat ) ),
  #          algorithm="squarified", fontsize.labels = 0,
  #          border.lwds=c( 1,0.25 ),
  #          border.col = c("black","black") )
  # dev.off()
  
  # Create table with sequence x V-gene ID to assign color hues
  airr_df_tmp = airr_df %>% 
    arrange(desc(duplicate_count)) #%>%
    #mutate(idx=paste0(alakazam::getGene(v_call), '_', sequence_id))
  
  nclones_title = ifelse(nrow(airr_df_tmp) >= sampsize, 
                         paste0(sampsize, ' randomly-selected sequences'),
                         paste0(nrow(airr_df_tmp), ' sequences'))
  
  tm_p1 = ggplot(airr_df_tmp) +
    geom_treemap(aes(area=duplicate_count, subgroup=v_gene, fill=color), color='black', size=0.5, show.legend=FALSE) + 
    scale_fill_manual(values=setNames(airr_df_tmp[['color']], airr_df_tmp[['color']])) +
    guides(fill="none") +
    ggnewscale::new_scale_fill() +
    geom_treemap(aes(area=duplicate_count, subgroup=v_gene, fill=v_gene), size=0, alpha=0, show.legend=TRUE) +
    labs(fill='V-Gene', 
         title=paste0('V-gene abundance - ', comp, '\n', nclones_title, '\n', subj, '; Patient: ', pat)) +
    scale_fill_manual(values=setNames(cc_df[['color']], cc_df[['gene']])) +
    guides(fill=guide_legend(override.aes = list(alpha=1))) +
    coord_equal()
  
  pdf(tmfile, width=10, height=10)
  print(tm_p1)
  dev.off()
} else{
  cat(paste0('NO SEQUENCES AVAILABLE IN THIS SAMPLE ( ', stringr::str_extract(dbfile, '[0-9]index[0-9]+_[IGHTRB]{3}'), ' )\n'))
}

