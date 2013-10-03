#Method to decode a phenotype name into a 0,1,2,etc string
if (!isGeneric("decodePhenotype")) 
{
  setGeneric("decodePhenotype", function(pheno.code, marker.names, partitions.per.marker=rep(2,length(marker.names))) standardGeneric("decodePhenotype"))
}

setMethod("decodePhenotype", signature(pheno.code = "character", marker.names = "character", partitions.per.marker = "numeric"),
		  function(pheno.code, marker.names, partitions.per.marker)
		  {
		  	partitions.per.marker <- as.integer(partitions.per.marker)
		  	if (! is.character(marker.names))
		  		stop("marker.names should be a character vector!")
		  	if (!is.integer(partitions.per.marker))
		  		stop("partitions.per.marker should be an integer vector!")
		  	
		  	if(length(partitions.per.marker)==1)
		  		partitions.per.marker <- rep(partitions.per.marker, length(marker.names))
		  	
		  	makePartition <- function(partitions)
		  	{
		  		stateList = c('-', sapply(1:(partitions-1), function(x){paste(rep('+',x), sep='', collapse='')}))
		  		stateList
		  	}
		  	partitionList <- lapply(2:max(partitions.per.marker), makePartition)
		  	names (partitionList) <- 2:max(partitions.per.marker)
		  	
		  	
		  	decodeOneChannel <- function(phenoPos)
		  	{
		  		phenoChar <- substr(pheno.code, phenoPos, phenoPos)
		  		partitions <- partitions.per.marker[phenoPos]
		  		markerName <- marker.names[phenoPos]
		  		partitionCodes <- partitionList[[as.character(partitions)]]
		  		
		  		if(phenoChar!='0')
		  			return (paste(markerName, partitionCodes[as.integer(phenoChar)], sep=''))
		  		else
		  			return ('')
		  	}
		  	return (paste(sapply(1:length(marker.names), decodeOneChannel), collapse=''))
		  })


## Example of applying to a whole Phenotypes object:
#sapply(res@PhenoCodes, function(x,y){decodePhenotype(y,x)}, res)
