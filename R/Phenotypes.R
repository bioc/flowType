setClass("Phenotypes",
         representation(
           CellFreqs="numeric", 
           MFIs="matrix",
           PhenoCodes="vector",
           PropMarkers="vector", 
           MFIMarkers="vector", 
           MarkerNames="vector", 
           Partitions="matrix",
           MaxPopSize="numeric",
           PartitionsPerMarker="numeric"))

setMethod("summary", signature(object="Phenotypes"), function(object){
  cat(sprintf("Phenotypes object identified by flowType with %d cell-frequency-based on %d MFI-based features.\n", length(object@CellFreqs), dim(object@MFIs)[1]*dim(object@MFIs)[2]));
})

setMethod("plot", signature(x="Phenotypes", y="flowFrame"), function(x, y, ...){
  RawData=y
  object=x
  X=exprs(RawData)[,object@PropMarkers]
  par(mfrow=c(ceiling(sqrt(length(object@PropMarkers))),ceiling(sqrt(length(object@PropMarkers)))))
  for (i in 1:length(X[1,])){
    plot(density(X[,i]), main='', xlab='', ylab='',axes=FALSE, lwd=2, ...);
    axis(1);
    axis(2);
    title(main=object@MarkerNames[object@PropMarkers[i]], ylab='Density', cex=2,cex.lab=1.5);
    for (j in 1:(max(object@Partitions[,i])-1))
        abline(v=max(X[which(object@Partitions[,i]==j),i]),col=2, lwd=2);
  }
})

setMethod("plot", signature(x="Phenotypes", y="numeric"), function(x, y, Frame, ...){
  plot(data.frame(exprs(Frame)[,x@PropMarkers]), col=getLabels(x, y), ...)
})

setMethod("plot", signature(x="Phenotypes", y="character"), function(x, y, Frame, ...){
	phenoCode <- encodePhenotype(pheno.string=y, marker.names=x@MarkerNames[x@PropMarkers])
	plot(data.frame(exprs(Frame)[,x@PropMarkers]), col=getLabels(x, phenoCode), ...)
})
