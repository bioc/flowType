\name{flowType}
\alias{flowType}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  flowType: Phenotyping Flow Cytometry Assays
}
\description{
  flowType uses a simple threshold, Kmeans, flowMeans or flowClust to partition every channel to a positive and a negative cell population.
  These partitions are then combined to generate a set of multi-dimensional phenotypes.
}
\usage{
flowType(Frame, PropMarkers=NULL, MFIMarkers=NULL, Methods='kmeans', MarkerNames=NULL, MaxMarkersPerPop=NULL, PartitionsPerMarker=2, Thresholds=NULL, MemLimit=4,verbose=FALSE)
#If upgrading from flowType 1.x to 2.x, please check documentation as some arguments have changed slightly.
}
\arguments{
  \item{Frame}{
    A flowFrame (after transformation) that is going to be phenotyped.
  }
  \item{PropMarkers}{
    A vector of the indexes or names of the markers to partition to specify phenotypes. If \code{NULL}, all markers in the frame will be used.
  }
  \item{MFIMarkers}{
    A vector of the indexes or names of the markers for which MFIs must
    be measured. If \code{NULL}, no markers will be used.
  }
  \item{Methods}{
    A single string specifying the method to use to determine thresholds for partitioning of markers.  Values can be "kmeans", "flowMeans", "flowClust", or "Thresholds". If "Thresholds" is specified, user-specified thresholds must be provided via the Thresholds parameter.
  }
  \item{MarkerNames}{
    A vector of names for the channels. If \code{NULL}, the names in \code{Frame} will be used.
  }
  \item{MaxMarkersPerPop}{
    An integer specifying the maximum number of markers to use to define populations (how "deep" to phenotype). This should be less than or equal to \code{PropMarkers}. If \code{NULL}, will default to the length of \code{PropMarkers}.
  }
  \item{PartitionsPerMarker}{
  	An integer or vector of integers specifing the number of partitions per marker. If a single integer, this number will be used for all markers. If a vector, the numbers will be matched with \code{PropMarkers} in order.
  }
  \item{Thresholds}{
  A list of vectors specifying per-channel thresholds. Each list item corresponds to one marker, and contains the threshold(s) for that marker. If only one vector is provided in the list, then those thresholds will be used for all markers. Otherwise, the list must be of the same length as PropMarkers.
  Note: if \code{Methods == 'thresholds'}, then \code{Thresholds} must be specified. If not, it is ignored.
  }
 \item{MemLimit}{
    Memory limit in GB. flowType will do a sanity check before executing, and if the total size of counts plus MFI values for all populations would exceed \code{MemLimit}, will not run.   
  } 
 \item{verbose}{
    Boolean variable. If TRUE, information about different processing
    tasks will be printed into the standard output.
  } 
}

\value{
  \item{\code{CellFreqs}:}{Object of class \code{"numeric"} containing
    the cell frequencies measured for each phenotype. Phenotype names are assigned as labels.}
  \item{\code{MFIs}:}{Object of class \code{"matrix"} containing the measured MFIs for each phenotype. Phenotype names are assigned as column labels and marker names as row labels.}
  \item{PropMarkers}{A vector of the indexes or names of the markers for which cell proportions must be measured.}
  \item{MFIMarkers}{A vector of the indexes or names of the markers for which MFIs must be measured.}
  \item{MarkerNames}{A vector of names for the channels. If \code{NULL},
    the names provided in \code{Frame} will be used.}
  \item{Partitions}{A matrix where each column shows the partitioning of
    the respective channel. 1 and 2 correspond to negative and positive,
  respectively.}
\item{PhenoCodes}{A vector of strings of length N (the number of markers) for each phenotype measured. For every phenotype, the character corresponding to a given marker can be 0, 1, 2, etc for
    neutral, negative, positive, bright, etc. See
    the provided vigenette for more details
    and examples.}
}

\references{
  Nima Aghaeepour, Pratip K. Chattopadhyay, Anuradha Ganesan, Kieran
  O'Neill, Habil Zare, Adrin Jalali, Holger H. Hoos, Mario Roederer,
  and Ryan R. Brinkman. Early Immunologic Correlates of HIV Protection
  can be Identified from Computational Analysis of Complex Multivariate
  T-cell Flow Cytometry Assays. submitted to Bioinformatics, 2011.
}
\author{
  Nima Aghaeepour
}

\examples{
#Load the library
library(flowType)
data(DLBCLExample)

#These markers will be analyzed
PropMarkers <- 3:5
MFIMarkers <- PropMarkers
MarkerNames <- c('FS', 'SS','CD3','CD5','CD19')

#Run flowType
Res <- flowType(DLBCLExample, PropMarkers, MFIMarkers, 'flowMeans', MarkerNames);

MFIs=Res@MFIs;
Proportions=Res@CellFreqs;
Proportions <- Proportions / max(Proportions)
names(Proportions) <- unlist(lapply(Res@PhenoCodes, 
                      function(x){return(decodePhenotype(
                      x,Res@MarkerNames[PropMarkers],
                      Res@PartitionsPerMarker))}))

#Select the 30 largest phenotypes
index=order(Proportions,decreasing=TRUE)[1:30]
bp=barplot(Proportions[index], axes=FALSE, names.arg=FALSE)
text(bp+0.2, par("usr")[3]+0.02, srt = 90, adj = 0, labels = names(Proportions[index]), xpd = TRUE, cex=0.8)
axis(2);
axis(1, at=bp, labels=FALSE);
title(xlab='Phenotype Names', ylab='Cell Proportion')

#These phenotype can be analyzed using a predictive model (e.g., classification or regression)
}

\keyword{ clustering }
\keyword{ classification }
\keyword{ HIV }
\keyword{ FlowCytData }
