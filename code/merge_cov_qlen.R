#!/usr/bin/env Rscript --vanilla
args = commandArgs(trailingOnly=TRUE)

### -- Set environment ---------------------------------------------------------
current <- getwd()
input1 <- args[1]   # mean qlen
input2 <- args[2]   # sum coverage
output <- args[3]

### -- Package -----------------------------------------------------------------
library(tidyverse)

### -- INFILES -----------------------------------------------------------------
infile <- read.table(input1, header=FALSE, sep="\t")
deffile <- read.table(input2, header=TRUE, sep="\t") 
deffile <- deffile %>% select(hit_id, numreads)

### -- Merging files -----------------------------------------------------------
# deffile %>% head()
# infile %>% head()
colnames(infile)[1] <- "HMMs_ID"
mergefile <- left_join(infile, deffile, by="HMMs_ID")


### -- OUTFILE -----------------------------------------------------------------
write.table(mergefile, output, sep="\t", quote = FALSE, col.names = FALSE, row.names = FALSE )


