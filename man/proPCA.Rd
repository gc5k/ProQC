\name{proPCA}
\alias{proPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PCA analysis for protein matrix
}
\description{
PCA analysis for protein matrix
}
\usage{
proQC(pMat, id, out)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pMat}{
Protein matrix for analysis. A data matrix.
}
  \item{ID}{
ID array for the data matrix. len(ID) = nrow(pMat).
}
  \item{out}{
root of the file name for output. A string.
}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
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
proPCA(pMat, id, "out")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ proQC }% use one of  RShowDoc("KEYWORDS")
