##
# Calculate "duplicate_count" column from tsv database
#
# By Oscar Ospina
# Created: Feb 02, 2026
# Modified: Feb 23, 2026
#

# Script CMD arguments
args = commandArgs(TRUE) # Read arguments into "args" vector
tab_file = args[1] # tsv file prepared with Change-O from IMGT txz
newtab_file = args[2] # tsv with duplicate_count column
#tab_file = './data/tsv/1index12_IGH_db-pass.tsv'

library('tidyverse')

# Load AIRR TSV
airr = read_tsv(tab_file, show_col_types=FALSE)

# Compute duplicate_count per unique sequence
dup_map = airr %>%
  group_by(sequence) %>%
  summarise(duplicate_count=n(), .groups="drop")

# Attach back to rows
airr_with_dup = airr %>%
  left_join(., dup_map, by='sequence') %>%
  dplyr::select(-1) %>%
  dplyr::distinct() %>%
  mutate(sequence_id=1:nrow(.)) %>%
  dplyr::relocate(sequence_id, .before=1)

# Write to file
write.table(airr_with_dup, file=newtab_file, quote=FALSE, sep='\t', row.names=FALSE, na='')

