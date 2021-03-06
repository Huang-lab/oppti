% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze.R
\name{oppti}
\alias{oppti}
\title{Outlier protein and phosphosite target identification}
\usage{
oppti(data, mad.norm = FALSE, cohort.names = NULL, panel = "global",
panel.markers = NULL, tol.nas = 20, ku = 6, miss.pstat = 0.4,
demo.panels = FALSE, save.data = FALSE, draw.sc.plots = FALSE,
draw.vi.plots = FALSE, draw.sc.markers = NULL,
draw.ou.plots = FALSE, draw.ou.markers = NULL, verbose = FALSE)
}
\arguments{
\item{data}{a list object where each element contains a proteomics data for
a different cohort (markers in the rows, samples in the columns) or a
character string defining the path to such data (in .RDS format).}

\item{mad.norm}{logical, to normalize the proteomes to have a unit Median
Absolute Deviation.}

\item{cohort.names}{character array.}

\item{panel}{a character string describing marker panel, e.g., 'kinases'.
Use 'global' to analyze all markers quantified across cohorts (default).
Use 'pancan' to analyze the markers commonly quantified across the cohorts.}

\item{panel.markers}{a character array containing the set of marker names
that user wants to analyze, e.g., panel.markers = c("AAK1", "AATK", "ABL1",
"ABL2", ...).}

\item{tol.nas}{a constant in [0,100], tolerance for the percentage of NAs
in a marker, e.g., tol.nas = 20 will filter out markers containing 20\% or
more NAs across samples.}

\item{ku}{an integer in [1,num.markers], upper bound on the number of
nearest neighbors of a marker.}

\item{miss.pstat}{a constant in [0,1], statistic to estimate potential
outliers. See `artImpute()`.}

\item{demo.panels}{logical, to draw demographics of the panel in each
cohort.}

\item{save.data}{logical, to save intermediate data (background inference
and dysregulation measures).}

\item{draw.sc.plots}{logical, to draw each marker's qqplot of observed vs
inferred (imputed) expressions.}

\item{draw.vi.plots}{logical, to draw each marker's violin plot of observed
vs imputed expressions.}

\item{draw.sc.markers}{character array, marker list to draw scatter plots}

\item{draw.ou.plots}{logical, to draw each marker's outlier prevalence
(by the percentage of outlying samples) across the cohorts.}

\item{draw.ou.markers}{character array, marker list to draw pan-cancer
outlier percentage plots}

\item{verbose}{logical, to show progress of the algorithm.}
}
\value{
dysregulation scores of every marker for each sample.

the imputed data that putatively represents the expressions of the
markers in the (matched) normal states.

the result of Kolmogorov-Smirnov tests that evaluates the
statistical significance of each marker's outlier samples.

a data list containing, for each cohort, the percentage of outlier
samples for every marker.

a data list containing, for each cohort, the outlier significance
threshold.
}
\description{
Find outlying markers and events across cancer types.
}
\examples{
set.seed(1)
dat = setNames(as.data.frame(matrix(runif(10*10),10,10),
row.names = paste('marker',1:10,sep='')), paste('sample',1:10,sep=''))
result = oppti(dat)
}
\seealso{
[artImpute()] for how to set `miss.pstat` and `ku`
}
