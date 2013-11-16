flowType <- function(Frame, 
					 PropMarkers=NULL, 
					 MFIMarkers=NULL, 
					 Methods='kmeans', 
					 MarkerNames=NULL, 
					 MaxMarkersPerPop=NULL, 
					 PartitionsPerMarker=2, 
					 Thresholds=NULL, 
					 MemLimit=4,
           verbose=FALSE)
  {
  ##############################################################################################################################
  # Argument processing
  ##############################################################################################################################
  ##If list of markers are not supplied, use all of the available channels.
  if(is.null(PropMarkers))
    PropMarkers=c(1:length(exprs(Frame)[1,]));

  if(is.null(MFIMarkers))
    MFIMarkers=vector();
  
## TODO: typecheck methods parameter.
  
  if (is.null(MaxMarkersPerPop))
    MaxMarkersPerPop <- length(PropMarkers)
  
  ##If list of markers are names, convert them to indexes.
  if (FALSE %in% is.numeric(PropMarkers))
    PropMarkers <- unlist(lapply(1:length(PropMarkers), function(i){which(PropMarkers[i]==colnames(exprs(Frame)))}))
  if (FALSE %in% is.numeric(MFIMarkers)&&length(MFIMarkers) > 0)
    MFIMarkers <- unlist(lapply(1:length(MFIMarkers), function(i){which(MFIMarkers[i]==colnames(exprs(Frame)))}))

  ##Parse method:
  if (length(Methods)==1)
  {
  	Methods=rep(Methods, length(PropMarkers));
  }else
  {
  	stop('Only one method may be specified')
  }
  
  
   if( Methods=='Thresholds' && (!is.list(Thresholds)) )
     stop('Thresholds must be provided as a list of vectors.')
    
  if(length(Thresholds) == 1)
  {
    if(length(unique(PartitionsPerMarker)) > 1)
      stop('When markers have different numbers of partitions, you must specify Thresholds on a per-marker basis.')
    if(length(Thresholds[[1]])!=length(PartitionsPerMarker[1]))
      stop('When a single vector is provided for Thresholds, it must contain exactly PartitionsPerMarker-1 Thresholds.')
    Thresholds <- rep(Thresholds, length(PropMarkers))
  }
  
  ##If PartitionsPerMarker is a single value, replicate for all PropMarkers:
  if (length(PartitionsPerMarker)==1)
    PartitionsPerMarker=rep(PartitionsPerMarker, length(PropMarkers));
  
  if (length(PartitionsPerMarker)!=length(PropMarkers))
    stop('PartitionsPerMarker must either be specified once for all markers, or be of the same length as PropMarkers.')
  
      
  if(length(Thresholds)==0 && 'Thresholds' %in% Methods)
  	stop('When Thresholds is specified as a method, You must provide Thresholds via the "Thresholds" argument.')
   
  #Get marker names if not provided:
  if(is.null(MarkerNames))
    MarkerNames <- as.vector(Frame@parameters@data$name)
  
  #If still NULL, just make it be PropMarkers:
  if(is.null(MarkerNames))
    MarkerNames <- PropMarkers
  
 
  ##############################################################################################################################
  # Memory check
  ##############################################################################################################################
    
  ##Before even starting, perform a sanity check for memory usage
  NumPops <- calcNumPops(PartitionsPerMarker, MaxMarkersPerPop)
  
  MemUse <- calcMemUse(NumPops, length(PropMarkers), length(MFIMarkers), nrow(Frame), MaxMarkersPerPop, max(PartitionsPerMarker)) / 10^9
  if (verbose)
    message(sprintf('Estimated memory required: %f GB of RAM', MemUse))
  if(MemUse > MemLimit)
    stop(paste('Calling flowType with these parameters would require', MemUse, 'GB of RAM, but MemLimit is', MemLimit, 'GB.\n Try reducing MaxMarkersPerPop or the number of MFIMarkers.'))
  
  ##############################################################################################################################
  # Perform partititoning
  ##############################################################################################################################
  
  ##express the flowFrame
  X <- exprs(Frame)[,PropMarkers]
  ##number of markers
  M <- ncol(X)
  ##number of events
  N <- nrow(X)
  ##the matrix that will hold the partitions
  Partitions <- matrix(0,M,N);
 
  ##Use the respective method.
  for (i in 1:M){
    
    #K-means partitioning:
    if (Methods[i]=='kmeans'){
      km<-kmeans(X[,i], PartitionsPerMarker[i], nstart=50)$cluster
      
      ##Sort clusters in increasing order by mean:
      means <- sapply(unique(km), function(x){mean(X[which(km==x),i])})
      names(means) <- unique(km)
      means <- sort(means)
      new.km <- rep(0, length(km))
      for (cluster.ind in 1:PartitionsPerMarker[i])
      {
        to.replace <- which(names(means)==cluster.ind)
        new.km[which(km==cluster.ind)] <- to.replace
      }
      
      Partitions[i,] <- new.km;
    }
    
    #flowMeans partitioning:
    if (Methods[i]=='flowMeans'){
      km<-flowMeans(X[,i], NumC=PartitionsPerMarker[i], MaxN=10, nstart=10)@Label
      ##Sort clusters in increasing order by mean:
      means <- sapply(unique(km), function(x){mean(X[which(km==x),i])})
      names(means) <- unique(km)
      means <- sort(means)
      new.km <- rep(0, length(km))
      for (cluster.ind in 1:PartitionsPerMarker[i])
      {
        to.replace <- which(names(means)==cluster.ind)
        new.km[which(km==cluster.ind)] <- to.replace
      }
      
      Partitions[i,] <- new.km;
    }

    if (Methods[i]=='flowClust'){
      res <- flowClust(Frame, varNames=colnames(exprs(Frame))[i], K = PartitionsPerMarker[i], level=1);
      
      km=map(res@z);
      km=replace(km, which(is.na(km)), 1)
      
            ##Sort clusters in increasing order by mean:
      means <- sapply(unique(km), function(x){mean(X[which(km==x),i])})
      names(means) <- unique(km)
      means <- sort(means)
      new.km <- rep(0, length(km))
      for (cluster.ind in 1:PartitionsPerMarker[i])
      {
        to.replace <- which(names(means)==cluster.ind)
        new.km[which(km==cluster.ind)] <- to.replace
      }
      
      Partitions[i,] <- new.km;
    }

    #If method is Thresholds, still calculate partition membership for later plotting and labelling purposes:
    if (Methods[i]=='Thresholds'){
    	for (Marker in 1:length(Thresholds))
    	{
    		marker.vec <- rep(1, ncol(Partitions))
    		for(partition in 1:(length(Thresholds[[Marker]])+1))
    		{
    			if(partition != length(Thresholds[[Marker]])+1)
    			{
    				marker.vec[which(X[,Marker] < Thresholds[[Marker]][partition])] <- partition
    			}
    			else
    			{
    				marker.vec[which(X[,Marker] >= Thresholds[[Marker]][partition-1])] <- partition
    			}
    		}
    		Partitions[Marker,] <- marker.vec
    	}
    }
  }

  ###Turn partitions into Thresholds to pass down to CPP code:
  partToThresh <- function(ThisChan, Partitions, PropMarkers, ThisExpr)
  {
  	#If Thresholds are pre-specified, return those directly:
  	if(Methods[ThisChan] == 'Thresholds')
  	{
  	  return(Thresholds[[ThisChan]])
  	}
    ThisPart <- t(Partitions)[, ThisChan]
    PartLabels <- unique(ThisPart)
  	
    #Find maxima of all partitions:
    MaxParts <- sapply(PartLabels, function(x){max(ThisExpr[which(ThisPart==PartLabels[x]),ThisChan])})
    Thresholds <- sort(MaxParts)[1:(length(PartLabels)-1)]
    
    Thresholds
  }
  
  
  Thresholds <- lapply(1:length(PropMarkers), partToThresh, Partitions, PropMarkers, X)

  for (i in 1:M){
  ##because flowMeans's 1D Thresholds are not perfectly aligning, we calculate the partitions one more time after the Thresholds are finalized.
    if (Methods[i]=='flowMeans'){
    	for (Marker in 1:length(Thresholds))
    	{
    		marker.vec <- rep(1, ncol(Partitions))
    		for(partition in 1:(length(Thresholds[[Marker]])+1))
    		{
    			if(partition != length(Thresholds[[Marker]])+1)
    			{
    				marker.vec[which(X[,Marker] < Thresholds[[Marker]][partition])] <- partition
    			}
    			else
    			{
    				marker.vec[which(X[,Marker] >= Thresholds[[Marker]][partition-1])] <- partition
    			}
    		}
    		Partitions[Marker,] <- marker.vec
    	}
    }
  }


  
  ##############################################################################################################################
  # Calculate MFIs and counts via CPP code and return
  ##############################################################################################################################
  
  #Call down  to cpp code to calculate:
  if(length(MFIMarkers) > 0)
  {
    MFIData <- matrix(exprs(Frame)[,MFIMarkers], ncol=length(MFIMarkers))
  }else
  {
    MFIData <- matrix()
  }
  res <- .Call('countCells', as.integer(PartitionsPerMarker), Thresholds, as.integer(MaxMarkersPerPop), PropMarkers, MFIData, X, NumPops, verbose)
  Counts <- res$counts;
  if (length(MFIMarkers)>0){
    Means <- res$mfis;
  	Means[which(Counts==0),]=NA
  }else
  {
    Means=matrix();
  }
  Codes <- res$codes;
  
  if (length(MFIMarkers)>0){
   
    colnames(Means)=colnames(exprs(Frame))[MFIMarkers];
  }
  
  Partitions=t(Partitions);
  return(new("Phenotypes", CellFreqs=Counts, PhenoCodes=Codes, MFIs=Means, PropMarkers=PropMarkers, MFIMarkers=MFIMarkers, MarkerNames=MarkerNames, Partitions=Partitions, PartitionsPerMarker=PartitionsPerMarker, Thresholds=Thresholds));
}
