adj2 <- adjustFractions2(cObjSng, cObjMul, sObj, maxCellsPerMultiplet=4, binary=F)



adjustFractions2 <- function(
  singlets, multiplets, swarm, binary = TRUE, maxCellsPerMultiplet=Inf
){
  medianCellNumber <- sampleType <- estimatedCellNumber <- NULL
  if(!is.matrix(swarm)) {
    fractions <- getData(swarm, "fractions")
  } else {
    fractions <- swarm
  }
  
  #calculate median cell number per singlet class
  cnc <- cellNumberPerClass(singlets, multiplets) %>%
    {setNames(pull(., medianCellNumber), pull(., class))}
  
  cnc <- cnc[match(colnames(fractions), names(cnc))]
  if(!identical(names(cnc), colnames(fractions))) stop("cnc name mismatch")
  
  #calculate cell number per multiplet
  cnm <- estimateCells(singlets, multiplets, maxCellsPerMultiplet=maxCellsPerMultiplet) %>%
    filter(sampleType == "Multiplet") %>%
    {setNames(pull(., estimatedCellNumber), pull(., sample))}
  
  cnm <- cnm[match(rownames(fractions), names(cnm))]
  if(!identical(names(cnm), rownames(fractions))) stop("cnm name mismatch")
  
  #adjust fractions
  frac.renorm <- t(t(fractions) / cnc)
  adjusted <- trunc(frac.renorm * cnm)
  if(binary) adjusted[adjusted > 0] <- 1
  return(adjusted)
}
