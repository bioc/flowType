flowType <- function(Frame, PropMarkers=NA, MFIMarkers=NA, Methods='kmeans', MarkerNames=NULL){
  if (length(Methods)==1)
    Methods=rep(Methods, length(PropMarkers));
  if (length(Methods)!=length(PropMarkers))
    stop('One method per marker is required (unless only one method is provided)')

  ##If list of markers are not supplied, use the all of the available channels.
  if(length(which(is.na(PropMarkers)==TRUE))>0)
    PropMarkers=c(1:length(exprs(Frame)[1,]));
  if(length(which(is.na(MFIMarkers)==TRUE))>0)
    MFIMarkers=PropMarkers;

  ##If list of markers are names, convert them to indexes.
  if (FALSE %in% is.numeric(PropMarkers))
    unlist(lapply(1:length(PropMarkers), function(i){which(PropMarkers[i]==colnames(exprs(Frame)))}))
  if (FALSE %in% is.numeric(MFIMarkers))
    unlist(lapply(1:length(MFIMarkers), function(i){which(MFIMarkers[i]==colnames(exprs(Frame)))}))
  
  
  ##express the flowFrame
  X <- exprs(Frame)[,PropMarkers]
  ##number of markers
  M <- length(X[1,]);
  ##number of events
  N <- length(X[,1]);
  ##the matrix that will hold the partitions
  Partitions<-matrix(0,M,N);

  ##if the method is a number (threshold), use it to partition the cells
  for (i in 1:M){
    if (is.numeric(Methods[i])){
      Partitions[i,which(X[,i]>as.numeric(Methods[i]))] <- 2
      Partitions[i,which(X[,i]<=as.numeric(Methods[i]))] <- 1
    }
  }

  ##Else, use the respective method.
  for (i in 1:M){
    if (Methods[i]=='kmeans'){
      km<-kmeans(X[,i], 2, nstart=50)$cluster;
      ##make sure 2 is the "positive" cluster
      if (mean(X[which(km==1),i])>mean(X[which(km==2),i])){
        km <- replace(km, which(km==1), 3);
        km <- replace(km, which(km==2), 1);
        km <- replace(km, which(km==3), 2);
      }
      Partitions[i,] <- km;
    }
    
    if (Methods[i]=='flowMeans'){
      km<-flowMeans(X[,i], NumC=2, MaxN=10, nstart=10)@Label;
      if (mean(X[which(km==1),i])>mean(X[which(km==2),i])){
        km <- replace(km, which(km==1), 3);
        km <- replace(km, which(km==2), 1);
        km <- replace(km, which(km==3), 2);
      }
      Partitions[i,] <- km;
    }

    if (Methods[i]=='flowClust'){
      res <- flowClust(Frame, varNames=colnames(exprs(Frame))[i], K = 2, level=1);
      ##o<-flowObj(res,Frame);
      ##m<-merge(o);
      ##m<-m[[2]];
      km=map(res@z);
      km=replace(km, which(is.na(km)), 1);
      if (length(which(km==1))<1)
        km[which.min(X[,i])]=1
      if (length(which(km==2))<1)
        km[which.min(X[,i])]=2

      if (mean(X[which(km==1),i])>mean(X[which(km==2),i])){
        km <- replace(km, which(km==1), 3);
        km <- replace(km, which(km==2), 1);
        km <- replace(km, which(km==3), 2);
      }
      Partitions[i,] <- km;
    }

  }

  ##Count the number of cells in every combination of the partitions above.
  Counts <- vector();
  nums=3^M
  Means <- matrix(0, (nums-1), length(MFIMarkers))
  Markers <- colnames(exprs(Frame))
  if (!is.null(MarkerNames))
    Markers <- MarkerNames
  Names <- vector()
  for (i in 1:(nums-1)){
    Names[i] <- getPopName(i, Markers[PropMarkers]);
    v <- c(as.vector(matrix(0,1, M-length(digitsBase(i,3)))),digitsBase(i,3));
    v <- v-1;
    index <- as.vector(matrix(1,1,N))
    for (j in 1:M){
      if (v[j]==0)
        next;
      if (v[j]==1)
        index <- index * (Partitions[j,]-1)
      if (v[j]==-1)
        index <- index * (Partitions[j,]-2) * -1
    }
    Counts[i]=length(which(index==1));

    ##And measure their MFI.
    if (Counts[i]==1)
      Means[i,]=exprs(Frame)[,MFIMarkers][which(index==1),];
    if (Counts[i]==0)
      Means[i,]=0
    if (Counts[i]>1)
      Means[i,]=colMeans(exprs(Frame)[which(index==1),MFIMarkers]);
  }
  names(Counts)=Names;
  rownames(Means)=Names;
  colnames(Means)=Markers[MFIMarkers];
  Partitions=t(Partitions);
  return(new("Phenotypes", CellFreqs=Counts, MFIs=Means, PropMarkers=PropMarkers, MFIMarkers=MFIMarkers, MarkerNames=MarkerNames, Partitions=Partitions));
}
