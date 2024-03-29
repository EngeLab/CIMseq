adj <- adjustFractions(cObjSng, cObjMul, sObj, maxCellsPerMultiplet=4, binary=F)


ncells <- rowSums(adj)
cOrder <- c(
  "SI.Paneth", "SI.Lgr5+", "SI.Lgr5+.Mki67+", "SI.TA.early", "SI.TA.late", "SI.Enterocytes", "SI.Goblet", 
  "Enteroendocrine", "Tufft", "Blood"
)

plotSwarmCircos2(sObj, cObjSng, cObjMul, weightCut=5, maxCellsPerMultiplet=4, h.ratio=0.95, connectionClass='Enteroendocrine', alpha=0.2, classOrder=cOrder)
dev2bitmap("Enteroendocrine_circos_colon_alpha02.pdf", type="pdfwrite")

ps <- calculateEdgeStats(sObj, cObjSng, cObjMul, maxCellsPerMultiplet=4, depleted=F)

swarm2 <- sObj
swarm2@fractions <- swarm2@fractions[ncells == 2,]

plotSwarmCircos2(swarm2, cObjSng, cObjMul, weightCut=2, maxCellsPerMultiplet=4, h.ratio=0.95, alpha=0.05, classOrder=cOrder)
dev2bitmap("swarm_2cells.pdf", type='pdfwrite')

swarm3 <- sObj
swarm3@fractions <- swarm3@fractions[ncells == 3,]
plotSwarmCircos2(swarm3, cObjSng, cObjMul, weightCut=2, maxCellsPerMultiplet=4, h.ratio=0.95, alpha=0.05, classOrder=cOrder)
dev2bitmap("swarm_3cells.pdf", type='pdfwrite')

swarm4 <- sObj
swarm4@fractions <- swarm4@fractions[ncells == 4,]
plotSwarmCircos2(swarm4, cObjSng, cObjMul, weightCut=2, maxCellsPerMultiplet=4, h.ratio=0.95, alpha=0.05, classOrder=cOrder)
dev2bitmap("swarm_4cells.pdf", type='pdfwrite')

swarm5 <- sObj
swarm5@fractions <- swarm5@fractions[ncells == 5,]
plotSwarmCircos2(swarm5, cObjSng, cObjMul, weightCut=2, maxCellsPerMultiplet=4, h.ratio=0.95, alpha=0.05, classOrder=cOrder)
dev2bitmap("swarm_5cells.pdf", type='pdfwrite')


swarm3up <- sObj
swarm3up@fractions <- swarm3up@fractions[ncells > 3,]
my.frac <- sum(ncells < 3) / sum(ncells > 3)
plotSwarmCircos2(swarm3up, cObjSng, cObjMul, weightCut=10/my.frac, maxCellsPerMultiplet=4, h.ratio=0.95, alpha=0.05, classOrder=cOrder)
dev2bitmap("swarm_4upcells_colon.pdf", type='pdfwrite')

swarm3down <- sObj
swarm3down@fractions <- swarm3down@fractions[ncells < 3,]
plotSwarmCircos2(swarm3down, cObjSng, cObjMul, weightCut=5, maxCellsPerMultiplet=4, h.ratio=0.95, alpha=0.05, classOrder=cOrder)
dev2bitmap("swarm_2cells_colon.pdf", type='pdfwrite')
