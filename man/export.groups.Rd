% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/export.groups.R
\name{export.groups}
\alias{export.groups}
\title{Export the IPCAPS result to a text file}
\usage{
export.groups(result.dir)
}
\arguments{
\item{result.dir}{A result directory as the \code{$output} object returned
from the \code{\link{ipcaps}} function.}
}
\value{
A data frame of exported data containing 4 columns; \code{group},
\code{node}, \code{label}, \code{row.number}, as described below for more
details:
\itemize{
\item \code{group} represents group membership of IPCAPS result.
\item \code{node} represents node numbers of IPCAPS result.
\item \code{label} represents labels of rows in orginal input data.
\item \code{row.number} represents row numbers of orginal input data.
}
}
\description{
Export clustering result of
\code{\link{ipcaps}} to text file called 'groups.txt'.
}
\details{
After running, this function exports the file called 'groups.txt' to the same result
directory. If 'groups.txt' already exists in the result directory, the
exported file is changed to 'groups1.txt', 'groups2.txt', 'groups3.txt', ...,
accordingly.
}
\examples{

# Importantly, bed file, bim file, and fam file are required
# Use the example files embedded in the package

BED.file <- system.file("extdata","ipcaps_example.bed",package="IPCAPS")
LABEL.file <- system.file("extdata","ipcaps_example_individuals.txt.gz",package="IPCAPS")

my.cluster <- ipcaps(bed=BED.file,label.file=LABEL.file,lab.col=2,out=tempdir())

#Here, to export the IPCAPS result to a text file
exported.data <- export.groups(my.cluster$output.dir)
print(dim(exported.data))
head(exported.data)
}
