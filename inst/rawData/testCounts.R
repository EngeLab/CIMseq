#requires sp.scRNAseqData package
#run from package root with: source('inst/rawData/testCounts.R')

packages <- c("sp.scRNAseq", "sp.scRNAseqData", "tidyverse")
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

s <- grepl("^s", colnames(countsSorted2))
m <- c('m.NJB00204.A02', 'm.NJB00204.G04', 'm.NJB00204.D07')
bool <- colnames(countsSorted2) %in% m

rv <- matrixStats::rowMaxs(countsSorted2[, s])
select <- order(rv, decreasing = TRUE)[1:2000]

c <- countsSorted2[select, s | bool]
e <- countsSortedERCC2[, s | bool]
out <- rbind(e, c)
out <- out[, rev(order(colnames(out)))]
write_tsv(matrix_to_tibble(out, "gene"), path = './inst/rawData/testCounts.txt')

#meta
countsSortedMeta2 %>% 
  filter(sample %in% colnames(out)) %>%
  select(sample, cellNumber, cellTypes) %>%
  write_tsv(path = "./inst/rawData/testMeta.txt")
