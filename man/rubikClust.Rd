\name{rubikClust}
\alias{rubikClust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
A function for unsupervised clustering to detect rough structures and outliers.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The function operates on a Nx3 matrix, where N is the number of samples and data are collected on 3 variables.
}
\usage{
rubikClust(X, min.space=0.4,rotation=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{x} here~~
A data matrix for which rows represent samples and the 3 columns represent features. Missingness is not allowed.
}
  \item{min.space}{
%%     ~~Describe \code{x} here~~
A value to specify a minimum space between 2 consecutive projected values. Default = 0.4.
}
  \item{rotation}{
%%     ~~Describe \code{x} here~~
To specify if rotation is enabled or not. Default = \emph{TRUE}.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
The function \emph{rubikClust} is able to take up to 3 variables (N x 3 matrix). In case, a matrix contains more than 3 columns, only the first three columns are used;  the other columns are ignored.
}
\value{

  The returned value is a vector of numbers representing cluster memberships.


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
groups <- rubikClust(PC[,1:3])
print(groups)

table(groups)
#groups
#  1   2   3 
#  4 751   5 
}

