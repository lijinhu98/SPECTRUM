\name{Gambler_R}
\alias{Gambler_R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
predict
}
\description{
predict
}
\usage{
Gambler_R(test_r, prob_r,gene_list = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test_r}{
  Rows should be cells and columns should be genes.
}
\item{prob_r}{
  trained scibet model  Rows should be genes and columns should be cell types.The genes must be matched with test_r
}
\item{gene_list}{
  a list of genes
}
}

\value{
return the predict
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}

