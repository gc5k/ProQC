\name{pcaCatPlot}
\alias{pcaCatPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pVCA(pEg, pc=5, Cat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pEg}{
pEg is the result of eigen analysis.
}
  \item{pc}{
The number of the eigenvectors used for analysis.
}
  \item{Cat}{
The variable for analysis. If it is a continous variable, linear regression is used. If it is a categorical variable, ANOVA is used to determine Rsq.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ proQC }% use one of  RShowDoc("KEYWORDS")
