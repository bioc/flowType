#Returns cluster membership labels of the given phenotype number.
getLabels <- function(Phenotypes, PhenotypeNumber){
  if (!is.numeric(PhenotypeNumber))
    PhenotypeNumber=which(names(Phenotypes@CellFreqs)==PhenotypeNumber);
  Partitions=Phenotypes@Partitions;
  M=length(Phenotypes@PropMarkers);
  N=length(Partitions[,1]);
  v <- c(as.vector(matrix(0,1, M-length(digitsBase(PhenotypeNumber,3)))),digitsBase(PhenotypeNumber,3));
  v <- v-1;
  index <- as.vector(matrix(1,1,N))
  for (j in 1:M){
    if (v[j]==0)
      next;
    if (v[j]==1)
      index <- index * (Partitions[,j]-1)
    if (v[j]==-1)
      index <- index * (Partitions[,j]-2) * -1
  }
  return(index+1);
}
