% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/analyze.R
\name{markOut}
\alias{markOut}
\title{Display outlying expressions}
\usage{
markOut(marker.proc.list, dat, dat.imp, dat.imp.test, dat.dys,
dys.sig.thr.upp, dataset, num.omit.fit = NULL, draw.sc = TRUE,
draw.vi = TRUE)
}
\arguments{
\item{marker.proc.list}{character array, the rownames of the data to be
processed for outlier analyses and for plotting.}

\item{dat}{an object of log2-normalized protein (or gene) expressions,
containing markers in rows and samples in columns.}

\item{dat.imp}{the imputed data that putatively represents the expressions
of the markers in the (matched) normal states.}

\item{dat.imp.test}{marker's p-value of the statistical significance between
its observed vs imputed values computed by the Kolmogorov-Smirnov test.}

\item{dat.dys}{samples' distances to regression line (i.e., dysregulation)
on the scatter plots.}

\item{dys.sig.thr.upp}{the dysregulation score threshold to elucidate/mark
significantly dysregulated outlier events.}

\item{dataset}{the cohort name to be used in the output files.}

\item{num.omit.fit}{number of outlying events to ignore when fitting a
marker's observed expressions to the imputed ones.}

\item{draw.sc}{logical, to draw a scatter plot for every marker in
marker.proc.list in a separate PDF file.}

\item{draw.vi}{logical, to draw a violin plot for every marker in
marker.proc.list in a separate PDF file.}
}
\value{
the scatter plots of the markers where the outlier dysregulation
events are highlighted by red mark.
}
\description{
Mark outlying expressions on the scatter plot of a given marker
}