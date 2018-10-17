#requires sp.scRNAseqData package
#run from package root with: source('inst/rawData/testCounts.R')

#counts
packages <- c("sp.scRNAseq", "sp.scRNAseqData", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

m <- c('m.NJB00204.A02', 'm.NJB00204.G04', 'm.NJB00204.D07')
lg1 <- colnames(SCM.Counts) %in% m
lg2 <- grepl("^s", colnames(SCM.Counts))
lg <- lg1 | lg2

rv <- matrixStats::rowMaxs(SCM.Counts[, lg2])
select <- order(rv, decreasing = TRUE)[1:2000]

testCounts <- rbind(SCM.CountsERCC[, lg], SCM.Counts[select, lg]) %>%
  .[, rev(order(colnames(.)))]

save(testCounts, file = "data/testCounts.rda", compress = "bzip2")

#meta
testMeta <- SCM.Meta %>%
  filter(sample %in% colnames(SCM.Counts)[lg]) %>%
  select(sample, cellNumber, cellTypes)

save(testMeta, file = "data/testMeta.rda", compress = "bzip2")
