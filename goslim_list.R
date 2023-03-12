#!/usr/bin/env Rscript
## CPCantalapiedra 2022

library(GSEABase)
library(GO.db)

args <- commandArgs(trailingOnly = TRUE)

gos_path = args[1]
fl <- "data/goslim_generic.obo"
slim <- getOBOCollection(fl)

mygos <- read.table(gos_path, sep = "\t", header = FALSE,
      row.names = NULL,
      col.names = c("e6_id", "go_list"),
      colClasses = c("character", "character"))

f_get_goslim <- function(e6_id, go_list) {
         # gos <- as.list(strsplit(go_list, split=",")) #, fixed=TRUE)# [[1]]
         gos <- unlist(strsplit(go_list, split=",")) #, fixed=TRUE)# [[1]]
         gos_col <- GOCollection(gos)
         # BP (biological process), MF (molecular function), or CC (cellular compartment)
         goslims <- goSlim(gos_col, slim, "BP")
         goslims <- goslims[goslims$Count > 0,]

         if (nrow(goslims) == 0) {
           goslims_ids = c("-")
           goslims_terms = c("-")
           goslims_counts = c("-")
         } else {
           goslims_ids = paste(rownames(goslims), collapse = ",")
           goslims_terms = paste(goslims$Term, collapse = "|")
           goslims_counts = paste(goslims$Count, collapse = "|")
         }

         write(sprintf("% s \t % s \t % s \t % s", e6_id, goslims_ids, goslims_terms, goslims_counts), stdout())

         return()
}

goslims <- mapply(f_get_goslim, mygos$e6_id, mygos$go_list)

quit()

