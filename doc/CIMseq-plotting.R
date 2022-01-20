## ---- message = FALSE---------------------------------------------------------
packages <- c(
  "CIMseq", "printr", "ggthemes", "dplyr", 
  "tidyr", "ggplot2", "viridis"
)
purrr::walk(packages, library, character.only = TRUE)
rm(packages)

## ---- fig.align="center", fig.width=10, fig.height=8, eval = FALSE------------
#  plotCountsMarkers(
#    CIMseqSinglets_test, CIMseqMultiplets_test, markers = c("CD74", "ANXA3")
#  )

## ---- fig.align="center", fig.width=10, fig.height=8--------------------------
plotCountsERCC(CIMseqSinglets_test, CIMseqMultiplets_test)

## ---- fig.align="center", fig.width=10, fig.height=8--------------------------
plotCountsMarkers(
  CIMseqSinglets_test, CIMseqMultiplets_test, markers = c("CD74", "ANXA3")
) %>%
  plotData() %>%
  head()

plotCountsERCC(CIMseqSinglets_test, CIMseqMultiplets_test) %>%
  plotData() %>%
  head()

## ---- fig.align="center", fig.width=10, fig.height=8--------------------------
plotUnsupervisedClass(CIMseqSinglets_test, CIMseqMultiplets_test)

## ---- fig.align="center", fig.width=10, fig.height=8--------------------------
plotUnsupervisedMarkers(CIMseqSinglets_test, CIMseqMultiplets_test, "CD74")

## ---- fig.align="center", fig.width=10, fig.height=8--------------------------
plotUnsupervisedMarkers(
  CIMseqSinglets_test, CIMseqMultiplets_test,
  markers = c("CD74", "ANXA3", "ACTG2")
)

## ---- fig.align="center", fig.width=10, fig.height=8, message = FALSE---------
#use the data instead
plotUnsupervisedMarkers(
  CIMseqSinglets_test, CIMseqMultiplets_test, 
  markers = c("CD74", "HLA-DRA", "IL13RA2", "MAGEA4")
) %>%
  plotData() %>%
  gather(gene, value, -Sample, -(`Sample type`:Colour)) %>%
  group_by(Sample) %>%
  mutate(`Mean(markers)` = mean(value)) %>%
  ungroup() %>%
  select(`dim.red dim 1`, `dim.red dim 2`, `Mean(markers)`) %>%
  distinct() %>%
  ggplot() +
  geom_point(aes(`dim.red dim 1`, `dim.red dim 2`, colour = `Mean(markers)`)) +
  viridis::scale_colour_viridis(option = "E") +
  theme_few() +
  theme(legend.position = "top", legend.title.align = 0) +
  guides(colour = guide_colourbar(title.position = "top"))

## -----------------------------------------------------------------------------
plotUnsupervisedClass(CIMseqSinglets_test, CIMseqMultiplets_test) %>% 
  plotData() %>%
  head()

plotUnsupervisedMarkers(
  CIMseqSinglets_test, CIMseqMultiplets_test, 
  markers = c("CD74", "ANXA3", "ACTG2")
) %>% 
  plotData %>%
  head()

## ---- fig.align="center", fig.width=10, fig.height=8--------------------------
plotSwarmCircos(swarm=CIMseqSwarm_test, singlets=CIMseqSinglets_test, multiplets=CIMseqMultiplets_test, weightCut=10, maxCellsPerMultiplet=4, alpha=Inf, h.ratio=0.9, depleted=F)

## ---- fig.align="center", fig.width=10, fig.height=8--------------------------
plotSwarmEdgeBar(CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test)

## ---- fig.align="center", fig.width=10, fig.height=8--------------------------
plotSwarmPbar(CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test)

## ---- fig.align="center", fig.width=10, fig.height=8--------------------------
plotSwarmHeat(CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test)

## -----------------------------------------------------------------------------
plotSwarmEdgeBar(CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test) %>%
  plotData() %>%
  head()

## ---- fig.align="center", fig.width=10, fig.height=8--------------------------
plotSwarmGenes(
  CIMseqSwarm_test, CIMseqSinglets_test, CIMseqMultiplets_test,
  c("CD74", "ANXA3", "ACTG2"), rownames(getData(CIMseqSwarm_test, "fractions"))
)

