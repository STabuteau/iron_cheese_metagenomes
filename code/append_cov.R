#!/usr/bin/env Rscript --vanilla
args = commandArgs(trailingOnly=TRUE)

### -- Set environment ---------------------------------------------------------
current <- getwd()
input1 <- args[1]   # hmmscan parsed results
input2 <- args[2]   # mapping coverage
output <- args[3]

### -- Package -----------------------------------------------------------------
library(tidyverse)

### -- INFILES -----------------------------------------------------------------
infile <- read.table(input1, header=FALSE, sep="\t")
deffile <- read.table(input2, header=TRUE, sep="\t") 
deffile <- deffile %>% select(rname, numreads, coverage, meandepth)

### -- Merging files -----------------------------------------------------------
# deffile %>% head()
# infile %>% head()
colnames(infile)[3] <- "rname"
mergefile <- left_join(infile, deffile, by="rname")


### -- OUTFILE -----------------------------------------------------------------
write.table(mergefile, output, sep="\t", quote = FALSE, col.names = FALSE, row.names = FALSE )


