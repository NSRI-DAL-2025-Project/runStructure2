% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runfaststructure.R
\name{gl.run.faststructure}
\alias{gl.run.faststructure}
\title{Runs a faststructure analysis using a genlight object}
\usage{
gl.run.faststructure(
  x,
  k.range,
  num.k.rep = 1,
  output = getwd(),
  tol = 1e-05,
  prior = "simple",
  cv = 0,
  seed = NULL
)
}
\arguments{
\item{x}{Name of the genlight object containing the SNP data [required].}

\item{k.range}{Range of the number of populations [required].}

\item{num.k.rep}{Number of replicates [default 1].}

\item{output}{Path to output file [default getwd()].}

\item{tol}{Convergence criterion [default 10e-6].}

\item{prior}{Choice of prior: simple or logistic [default "simple"].}

\item{cv}{Number of test sets for cross-validation, 0 implies no CV step
[default 0].}

\item{seed}{Seed for random number generator [default NULL].}

\item{exec}{Full path and name+extension where the fastStructure executable
is located [default working directory "./fastStructure"].}

\item{exec.plink}{path to plink executable [default working directory].}
}
\description{
This function takes a genlight object and runs a faststructure analysis.
}
\details{
Download faststructure binary for your system from here (only runs on Mac or
Linux):
}
\examples{
\dontrun{
# Please note: faststructure needs to be installed
# Please note: faststructure is not available for windows
t1 <- gl.filter.callrate(platypus.gl, threshold = 1)
res <- gl.run.faststructure(t1,
  exec = "./fastStructure", k.range = 2:3,
  num.k.rep = 2, output = paste0(getwd(), "/res_str")
)
qmat <- gl.plot.faststructure(res, k.range = 2:3)
gl.map.structure(qmat, K = 2, t1, scalex = 1, scaley = 0.5)
}
}
\author{
Luis Mijangos (Post to \url{https://groups.google.com/d/forum/dartr})
}
