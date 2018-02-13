\name{pEvecCat}
\alias{pEvecCat}
\title{
PCA illustration for categories
}
\description{
PCA illustration for different kind of samples.
}
\usage{
pEvecCat(pEvec, cat, PC=c(1,2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pEvec}{
Eigenvector matrix.
}
  \item{cat}{
Array for different samples. length(cat) = nrow(pEvec).
}
  \item{PC}{
Index for eigenvectors for illustration. length(PC) = 2.
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
Guo-Bo Chen
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
pEvecCat(pEvec, cat, c(1,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ EigenGWAS }% use one of  RShowDoc("KEYWORDS")
