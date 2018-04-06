\name{plot.3views}
\alias{plot.3views}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
A function to create scatter plots in three views.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function visualizes data in X-Y plane, X-Z plane, and Y-Z plane; an input object (\emph{matrix} or \emph{data.frame}) must contain at least 3 columns.
}
\usage{
plot.3views(X, labels, col.pat.table=NA, only.row=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{x} here~~
A matrix or a data.frame that contains at least 3 columns of numeric data. If there are more than 3 columns in \emph{X}, only the first 3 columns will be used.  
}
  \item{labels}{
%%     ~~Describe \code{x} here~~
A vector containing row labels of \emph{X} for display. All vector elements should be of type "character" (\emph{as.character}). The length of this vector equals the number of rows in \emph{X}.
}
  \item{col.pat.table}{
%%     ~~Describe \code{x} here~~
A data.frame that associates colors and patterns to row labels (see \emph{labels} before). It needs to contain 3 columns which represent unique labels (column 1) and associated patterns (column 2) and colors (column 3). Default = \emph{NA}.
}
  \item{only.row}{
%%     ~~Describe \code{x} here~~
A vector that contains subset of row numbers that are selected to be plotted. Default = \emph{NA}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
\emph{col.pat.table} needs to be set properly. For example, given the vector \emph{labels} consists of 5 elements as in \emph{c('pop1','pop2','pop2','pop2','pop1')}, the data.frame \emph{col.pat.table} must contain 2 rows corresponding to 'pop1' and 'pop2', and 3 columns to describe their characteristics (i.e. patterns and colors). See the section "Examples" for more details.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

data(simSNP) #to load simulated data set
plot.3views(PC[,1:3],label)

#to change colors and patterns
all.labels <- unique(label)
my.colors <- c('pink','yellow','cyan','green')
my.patterns<- c('o','x','&','#')
my.table <- cbind(all.labels,my.patterns,my.colors)
plot.3views(PC[,1:3],label,col.pat.table=my.table)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

