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
flowType(Frame, PropMarkers = NA, MFIMarkers = NA, Methods = "kmeans", MarkerNames="NULL")
}
\arguments{
  \item{Frame}{
    A flowFrame (after transformation) that is going to be phenotyped.
  }
  \item{PropMarkers}{
    A vector of the indexes or names of the markers for which cell proportions must be measured.
  }
  \item{MFIMarkers}{
    A vector of the indexes or names of the markers for which MFIs must be measured.
  }
  \item{Methods}{
    A vector of strings of length 1 or PropMarkers. Values can be "Kmeans", "flowMeans", "flowClust", or a numerical value. If only one method is provided the same method will be used for all of the channels.
  }
  \item{MarkerNames}{
    A vector of names for the channels. If \code{NULL}, the names provided in \code{Frame} will be used.
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
MarkerNames <- c('Time', 'FSC-A','FSC-H','SSC-A','IgG','CD38','CD19','CD3','CD27','CD20', 'NA', 'NA')

#These markers will be analyzed
PropMarkers <- 3:5
MFIMarkers <- PropMarkers
MarkerNames <- c('FS', 'SS','CD3','CD5','CD19')

#Run flowType
Res <- flowType(DLBCLExample, PropMarkers, MFIMarkers, 'flowMeans', MarkerNames);

MFIs=Res@MFIs;
Proportions=Res@CellFreqs;
Proportions <- Proportions / max(Proportions)

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
