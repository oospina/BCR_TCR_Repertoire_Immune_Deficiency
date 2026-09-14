# -----------------------------------------------------------------------------
# treemap_clonal.R
# version: 1.0 (2020-05-28)
# author : Peter Blazso, Krisztian Csomos
#
# Generates uniformly color-coded treemap of V-gene usage/abundance 
# granularity: unique clonotypes (instead of unique sequences)
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
#dbfile = './data/clone/1index12_IGH_allclone.tsv'
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

if(nrow(airr_df) > 0){ # In case no clones are present for a sample (e.g., 2index25 TRB)
  # Get index-patient key
  pat = read.table(idfile)
  pat = pat[[2]][pat[[1]] == subj]
  
  # downsample AIRR data if original sample size is greater than limit
  #if( nrow (airr_df) > sampsize )
  #  airr_df <- airr_df %>% sample_n( sampsize )
  
  # assign uniform hue values to different V-genes
  airr_df$v_gene <- getGene( airr_df$v_call )
  airr_df   <- airr_df %>% left_join( cc_df, by=c("v_gene"="gene") )
  
  # generate a clone size distribution table with V-genes and hues
  clonesizes <- airr_df %>% select( clone_id, duplicate_count, v_gene, hue ) %>%
    group_by( clone_id, v_gene, hue ) %>%
    summarize( sumcount=sum(duplicate_count)) %>%
    arrange( desc(sumcount) ) %>%
    head( sampsize )
  
  # mark distinct sequences in clonal groups by (brightness) value levels
  vals <-  seq( 1, 0.2, by=-0.01 )
  vals = vals[round(seq(1, length(vals), length.out=length(clonesizes$hue)))] # Select equally spaced hues because some samples have fewer sequences
  clonesizes$color <- hex(HSV(clonesizes$hue, 0.8, vals))
  
  # plot treemap
  # tiff( filename=tmfile, width=3600, height=2880, res=600, compression="zip" )
  # treemap( clonesizes, index=c("v_gene","clone_id"),
  #          vSize="sumcount", vColor="color",
  #          type="color",
  #          title="", # title=paste( "AIRR TreeMap","-",paste0( subj,"_",comp ) ),
  #          algorithm="squarified", fontsize.labels = 0,
  #          border.lwds=c( 1,0.25 ),
  #          border.col = c("black","black") )
  # dev.off()
  
  # Create table with sequence x V-gene ID to assign color hues
  clonesizes_tmp = clonesizes %>% 
    arrange(desc(sumcount)) #%>%
    #mutate(idx=paste0(alakazam::getGene(v_gene), '_', clone_id))
  
  
  nclones_title = ifelse(nrow(clonesizes_tmp) >= sampsize, 
                         paste0(sampsize, ' top most abundant clones'),
                         paste0(nrow(clonesizes_tmp), ' clones'))
  
  # Create treemap
  tm_p1 = ggplot(clonesizes_tmp) +
    geom_treemap(aes(area=sumcount, subgroup=v_gene, fill=color), color='black', size=0.5, show.legend=FALSE) + 
    scale_fill_manual(values=setNames(clonesizes_tmp[['color']], clonesizes_tmp[['color']])) +
    guides(fill="none") +
    ggnewscale::new_scale_fill() +
    geom_treemap(aes(area=sumcount, subgroup=v_gene, fill=v_gene), size=0, alpha=0, show.legend=TRUE) +
    labs(fill='V-Gene', 
         title=paste0('V-gene abundance - ', comp, '\n', nclones_title, '\n', subj, '; Patient: ', pat)) +
    scale_fill_manual(values=setNames(cc_df[['color']], cc_df[['gene']])) +
    guides(fill=guide_legend(ncol=3, override.aes=list(alpha=1))) +
    coord_equal()
  
  pdf(tmfile, width=10, height=10)
  print(tm_p1)
  dev.off()
}

