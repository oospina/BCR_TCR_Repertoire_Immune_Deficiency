# ----------------------------------------------------------------------------
# visualize_nd.R
# version: 1.00 (2020-11-11)
# author : Peter Blazso, Kriszitan Csomos
#
# Load and visualize AIRR network diagrams
# -----------------------------------------------------------------------------

library(dplyr)
library(igraph)
library(colorspace)

divergingValueToColor <- function( x, lim=3, gran=10  )
{
  origo = gran+1
  # generate pre-tested Blue-Gray-Red palette
  colcount = gran*2+1
  pal <- diverging_hcl(n=colcount,h1=255,h2=12,c1=130,
                       cmax=80,l1=31,l2=78,p1=1,p2=1.3)

  # set extreme thresholds (values must fall between these two)
  x[ x > lim ] = lim
  x[ x < (-lim) ] = -lim

  # calculate stepping
  step <- lim/gran

  # lookup colors from palette and return this vector
  return( pal[ origo+round(x/step) ] )
}


# isotype_cols <- c( hIGHD="#807F7F",
#                    hIGHM="#CCCCCB",
#                    hIGHG="#FC8008",
#                    hIGHA="#FDCC65",
#                    hIGHE="#7F7F03" )


# MAIN CODE -------------------------------------------------------------------

# input variables (mostly generated from parsed input filename)
args <- commandArgs(TRUE)  # read arguments into "args" vector
grfile <- args[1]          # input GML graph data filename           - 1st arg
#grfile = './data/graph/1index12_IGH.gml'
ndfile <- args[2]          # output SVG1 (compartmental affiliation) - 2nd arg
agfile <- args[3]          # output SVG2 (antigene selection)        - 3rd arg
#itfile <- args[4]          # output SVG3 (isotypes)                  - 4rd arg
itfile = args[4]           # output V-gene clone
mmfile <- args[5]          # output SVG4 (mutations)                 - 5th arg
clfile = args[6]           # V gene colors - 6th arg
#clfile = './data/dest/IGH_v_gene_colorcodes.tsv'

clone_cols = readr::read_delim(clfile, delim='\t')
clone_cols = setNames(clone_cols[['color']], nm=clone_cols[['gene']])

grfile_parts <- strsplit(basename(grfile), "_")[[1]]
subj <- grfile_parts[1]    # subject ID
comp <- grfile_parts[2]    # cellular compartment
rep <- paste0(subj,"_",comp)


# read graph from gml file
graph <- read_graph(grfile, format="gml") # read.graph()` was deprecated

# ---- prepare graph for plot -------------------------------------------------

# disassemble graph - add weight to "destination" vertices - reassemble graph
edf <- igraph::as_data_frame(graph,"edges")
edf$label <- NULL
vdf <- igraph::as_data_frame(graph,"vertices") %>%
       left_join( edf, by=c("name"="to"))      %>%  select( -from )
vdf$weight[ is.na(vdf$weight) ] = 0
V(graph)$weight = vdf$weight

# correct for negative values in case of putative "Inferred" sequences
V(graph)$duplicatecount[ V(graph)$duplicatecount < 0 ] = 1

# calculate node sizes based on duplicate copy number
V(graph)$size = log10( V(graph)$duplicatecount )

# node counts
nodes_real <- sum( !grepl( "Inferred", V(graph)$name ) )
nodes      <- length(V(graph))


# ---- plot graphs with the same layout ---------------------------------------
# recalculate weights to visualize more spread-out graphs
lo <- layout_with_fr( graph, dim=2, weights=E(graph)$weight/4,
                      grid="nogrid", niter = vcount(graph)*2 )

# ---- create SVG file for compartmental view
#svg( filename=ndfile )
graphics.off()
pdf(ndfile, height=8, width=8)

# set colors for plain, (inter-)compartmental view
V(graph)$color = "darkgray"
V(graph)$color[ V(graph)$repertoireid == "MN" ] = "#0000FF"
V(graph)$color[ V(graph)$repertoireid == "AN" ] = "#FF0000"
V(graph)$color[ V(graph)$repertoireid == "MB" ] = "#007D00"

# make the plot
plot(graph, vertex.label=NA, vertex.frame.width=0.2,
     vertex.frame.color=NA,
     edge.width=0.7, edge.label=NA, edge.color = "darkgrey",
     layout=lo, main=paste(rep, "- AIRR network diagram") )

# save plot
dev.off()

# ---- create SVG file for "antigen-selection" view
#svg( filename=agfile )
graphics.off()
pdf(agfile, height=8, width=8)

# recolor vertices based on BASELINe sigma values
V(graph)$color = divergingValueToColor( as.numeric(V(graph)$baselinesigma) )

# make the plot with the same layout
plot(graph, vertex.label=NA, vertex.frame.width=0.2,
     vertex.frame.color=NA,
     edge.width=0.7, edge.label=NA, edge.color = "darkgrey",
     layout=lo, main=paste(rep, "- AIRR network diagram (BASELINe sigma)") )

# save plot
dev.off()

# ---- create SVG file for "isotype" view
#svg( filename=itfile )
graphics.off()
pdf(itfile, height=8, width=8)

# recolor vertices based on isotypes
#V(graph)$color = isotype_cols[ substr( V(graph)$isotype, 0,5 ) ]
V(graph)$color = clone_cols[ substr( V(graph)$v_gene, 0,5 ) ]
V(graph)$color[ is.null(V(graph)$color) ] = "darkgray"

# make the plot with the same layout
plot(graph, vertex.label=NA, vertex.frame.width=0.2,
     vertex.frame.color=NA,
     edge.width=0.7, edge.label=NA, edge.color = "darkgrey",
     layout=lo, main=paste(rep, "- AIRR network diagram (V-gene)") )

# save plot
dev.off()

# ---- create SVG file for "mutation_map" view
#svg( filename=mmfile )
graphics.off()
pdf(mmfile, height=8, width=8)

# recolor vertices based on mutation count
V(graph)$color = divergingValueToColor( V(graph)$weight, lim=10 )

# make the plot with the same layout
plot(graph, vertex.label=NA, vertex.frame.width=0.2,
     vertex.frame.color=NA,
     edge.width=0.7, edge.label=NA, edge.color = "darkgrey",
     layout=lo, main=paste(rep, "- AIRR network diagram (mutations)") )

# save plot
dev.off()
