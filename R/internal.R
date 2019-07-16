# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#   cbindNA
#
#   Binds arbitrary-length columns given in a list into a dataframe with
#     NA-filled ends.
#
cbindNA = function(X = list(seq_len(2),seq_len(3),seq_len(4))){
    n = 0
    for (x in X){
        n = max(n, length(x))
    }
    df = data.frame(matrix(NA, nrow = n, ncol = length(X)))
    for (c in seq_len(length(X))){
        df[seq_len(length(X[[c]])), c] = X[[c]]
    }
    return(df)
}
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#   outScores
#
#   Calculates potential outlier scores (p-statistic) as a function of distance
#     from 1st/3rd (box) quartiles.
#   If the value lies within the box the score is 1, if it lies further away
#     from the box the score gets closer to 0.
#   Returns, an object of the same class as dat containing the potential
#     outlier scores.
#
outScores = function(dat) {
    n = ncol(dat)
    dat.IQR = matrix(rep(apply(dat, 1, 'IQR', na.rm = TRUE), each = n),
        ncol = n, byrow = TRUE)
    dat.q75 = matrix(rep(apply(dat, 1, function(x){quantile(x, .75,
        na.rm = TRUE)}), each = n), ncol = n, byrow = TRUE)
    dat.q25 = matrix(rep(apply(dat, 1, function(x){quantile(x, .25,
        na.rm = TRUE)}), each = n), ncol = n, byrow = TRUE)
    dat.q75iqr.upp = (dat - dat.q75) / dat.IQR
    dat.q75iqr.low = (dat.q25 - dat) / dat.IQR
    # conflate over/under-expressed outliers
    outlier.pstat.mat = pmax(dat.q75iqr.upp, dat.q75iqr.low) + 1
    outlier.pstat.mat[is.na(outlier.pstat.mat)] = 1
    outlier.pstat.mat[outlier.pstat.mat<1] = 1
    outlier.pstat.mat = 1 / outlier.pstat.mat
    return(outlier.pstat.mat)
}
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#   gqplot
#
#   Draws a scatter plot based on the grammer of graphics (ggplot2), then
#     draws a regression line and displays a confidence interval around it.
#
#   Returns each data point's distance to the regresion line, and the plot.
#
gqplot = function(y, x, ci = 0.95, xlab = NULL, ylab = NULL, dist.sort = FALSE,
    d.thr = 0, na.action = 'omit', samp.names = NULL, marker.name = NULL,
    highlight = NULL, omit.fit = NULL, maxx = NULL, minx = NULL, maxy = NULL,
    miny = NULL, align.xy = TRUE) {
    # marker name
    if (is.null(marker.name)) {marker.name = rownames(y)}
    if (is.null(marker.name)) {marker.name = rownames(x)}
    # align dimentions
    if (nrow(as.matrix(y)) < ncol(as.matrix(y))) {y = t(y)}
    if (nrow(as.matrix(x)) < ncol(as.matrix(x))) {x = t(x)}
    # drop NA values
    if ((any(is.na(y)) | any(is.na(x))) & na.action == 'omit'){
        keep = !is.na(y) & !is.na(x)
        x = x[keep]
        y = y[keep]
        if (!is.null(highlight))
            {highlight[highlight %in% rownames(keep)[keep]]}
    }
    # sample names
    if (is.null(samp.names)) {samp.names = names(y)}
    if (is.null(samp.names)) {samp.names = names(x)}
    # align classes
    if (!methods::is(y)[[2]]=='vector') {y = as.vector(y)}
    if (!methods::is(x)[[2]]=='vector') {x = as.vector(x)}
    # define labels
    if (is.null(xlab)) {xlab = expression('Variable2')}
    if (is.null(ylab)) {ylab = expression('Variable1')}
    # find linear regression fit and compute distances to regression line
    if (is.null(omit.fit)) {
        fit = stats::lm(y ~ x)
        preds = predict(stats::lm(y ~ x, na.action = na.exclude),
                        interval = 'confidence', level = ci)
    } else if (omit.fit < 1) {
        fit = stats::lm(y ~ x)
        preds = predict(stats::lm(y ~ x, na.action = na.exclude),
                        interval = 'confidence', level = ci)
    } else {
        y.s = sort(y, decreasing = TRUE, index.return = TRUE, na.last = TRUE)
        non.out = (y < y[y.s$ix[omit.fit]]) &
            (y > y[y.s$ix[length(y[!is.na(y)])-omit.fit+1]])
        fit = stats::lm(y[non.out] ~ x[non.out], na.action = na.exclude)
        preds = matrix(NA, length(y), 3)
        preds[non.out,] = predict(fit, interval = 'confidence', level = ci)
    }
    # y = slope x + intercept; ax + by + c = 0; by = -ax -c;
    # d_x0_y0 = |ax0 + by0 + c| / sqrt(a^2 + b^2)
    # b = 1; a = -slope; c = -intercept
    d = (y - fit$coefficients[2]*x - fit$coefficients[1]) /
        sqrt(fit$coefficients[2]**2+1)
    if (dist.sort) {
        sd = sort(d, decreasing = TRUE, index.return = TRUE, na.last = TRUE)
        outlier.score = data.frame(sampleID = samp.names[sd$ix[sd$x >=
        quantile(sd$x, d.thr, na.rm = TRUE)]], dist2reg = sd$x[sd$x >=
        quantile(sd$x, d.thr, na.rm = TRUE)])
    } else {
        outlier.score = data.frame(sampleID = samp.names[d >=
        quantile(d, d.thr, na.rm = TRUE)], dist2reg = d[d >=
        quantile(d, d.thr, na.rm = TRUE)])
    }
    # when there is NA, preds is shorten, handle NAs: (NA values dropped)
    df = data.frame(variable1 = y, variable2 = x, clower = preds[,2],
        cupper = preds[,3])
    rownames(df) = samp.names
    # newx = seq(min(df$variable2), max(df$variable2),
        # length.out=length(df$variable2))
    if (is.null(minx)) {minx = min(df$variable2, na.rm=TRUE)}
    if (is.null(maxx)) {maxx = max(df$variable2, na.rm=TRUE)}
    if (is.null(miny)) {miny = min(df$variable1, na.rm=TRUE)}
    if (is.null(maxy)) {maxy = max(df$variable1, na.rm=TRUE)}
    if (align.xy) {
        minx = min(minx, miny); miny = minx
        maxx = max(maxx, maxy); maxy = maxx
    }
    if (is.null(highlight)) {
        gg = ggplot2::ggplot(df, ggplot2::aes(x=variable2, y=variable1)) +
            ggplot2::geom_point(ggplot2::aes(variable2, variable1),
                shape = 1, size = 3) +
            ggplot2::geom_abline(intercept = fit$coefficients[1],
                slope = fit$coefficients[2], alpha = 0.5) +
            ggplot2::geom_line(ggplot2::aes(variable2, cupper), size = .1) +
            ggplot2::geom_line(ggplot2::aes(variable2, clower), size = .1) +
            ggplot2::xlab(xlab) + ggplot2::ylab(ylab) +
            ggplot2::xlim(minx,maxx) + ggplot2::ylim(miny,maxy) +
            ggplot2::ggtitle(marker.name)
    } else {
        gg = ggplot2::ggplot(df, ggplot2::aes(x=variable2, y=variable1)) +
            ggplot2::geom_point(ggplot2::aes(variable2, variable1),
                shape = 1, size = 3) +
            ggplot2::geom_abline(intercept = fit$coefficients[1],
                slope = fit$coefficients[2], alpha = 0.5) +
            ggplot2::geom_line(ggplot2::aes(variable2, cupper), size = .1) +
            ggplot2::geom_line(ggplot2::aes(variable2, clower), size = .1) +
            ggplot2::xlab(xlab) + ggplot2::ylab(ylab) +
            ggplot2::xlim(minx,maxx) + ggplot2::ylim(miny,maxy) +
            ggplot2::ggtitle(marker.name) +
            ggplot2::geom_point(data = df[highlight,],
                ggplot2::aes(x=variable2, y=variable1), colour = 'red')
    }
    return(list(outlier.score, gg))
}
# _|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|_|
#
#   madNorm
#
#   Normalizes the columns of a data frame to have the unit Median Absolute
#     Deviation (MAD), i.e., median(abs(df$i-median(df$i))) = 1, for every 'i'.
#     Optionally, center the medians to the median of all column medians
#     (centering = T).
#
#   Special attention must be paid to columns with zero-MAD, i.e., no
#     variation, static values. By default, they are omitted (unprocessed).
#     Optionally, they can be centered to the median of all column medians
#     (centering.zero.mad = T), however, they can not be scaled to have a
#     unit-MAD due to static values.
#
madNorm = function(df, centering = FALSE, centering.zero.mad = FALSE){
    mads = apply(df,2,function(x)
        {median(abs(x-median(x,na.rm=TRUE)),na.rm=TRUE)})
    omit = mads==0
    mads[omit] = 1 # omit scaling 0-MAD columns
    df = df / matrix(rep(mads,each=dim(df)[1]),nrow=dim(df)[1],byrow=FALSE)
    if (centering) {
        cent = apply(df,2,'median',na.rm=TRUE)-
            median(apply(df,2,'median',na.rm=TRUE),na.rm=TRUE)
        if (!centering.zero.mad) {cent[omit] = 0} #omit centering 0-MAD columns
        df = df - matrix(rep(cent,each=dim(df)[1]),nrow=dim(df)[1],byrow=FALSE)
    }
    return(df)
}
