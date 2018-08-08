#requires sp.scRNAseqData package
#run from package root with: source('inst/rawData/testCounts.R')

#counts
packages <- c("sp.scRNAseq", "sp.scRNAseqData", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

m <- c('m.NJB00204.A02', 'm.NJB00204.G04', 'm.NJB00204.D07')
lg1 <- colnames(countsSorted2) %in% m
lg2 <- grepl("^s", colnames(countsSorted2))
lg <- lg1 | lg2

rv <- matrixStats::rowMaxs(countsSorted2[, lg2])
select <- order(rv, decreasing = TRUE)[1:2000]

testCounts <- rbind(countsSortedERCC2[, lg], countsSorted2[select, lg]) %>%
  .[, rev(order(colnames(.)))]

save(testCounts, file = "data/testCounts.rda")

#meta
testMeta <- countsSortedMeta2 %>%
  filter(sample %in% colnames(countsSorted2)[lg]) %>%
  select(sample, cellNumber, cellTypes)

save(testMeta, file = "data/testMeta.rda")
