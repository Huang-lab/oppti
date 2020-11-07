#' @title Filter out markers
#' @description Filters out markers based on the percentage of missing values,
#' low-expression and low-variability rates.
#' @param dat an object of log2-normalized protein (or gene) expressions,
#' containing markers in rows and samples in columns.
#' @param percent_NA a constant in [0,1], the percentage of missing values
#' that will be tolerated in the filtered data.
#' @param low_mean_and_std a constant in [0,inf], the lower-bound of the
#' mean or standard deviation of a marker in the filtered data.
#' @param q_low_var a constant in [0,1], the quantile of marker variances
#' which serves as a lower-bound of the marker variances in the filtered data.
#' @param force_drop character array containing the marker names that user
#' specifically wants to filter out.
#' @return filtered data with the same format as the input data.
#' @return the row names (markers) of the data that are filtered out due to
#' low-expression or low-variability.
#' @examples
#' dat = setNames(as.data.frame(matrix(1:(5*10),5,10),
#' row.names = paste('marker',1:5,sep='')), paste('sample',1:10,sep=''))
#' dat[1,1:2] = NA # marker1 have 20% missing values
#' dropMarkers(dat, percent_NA = .2) # marker1 is filtered out
dropMarkers = function(dat, percent_NA = .2, low_mean_and_std = .05,
    q_low_var = .25, force_drop = NULL){
    if (!is.null(force_drop)) {dat = dat[!(rownames(dat) %in% force_drop),]}
    dat = dat[rowSums(is.na(dat))/ncol(dat) <= percent_NA,]
    if (is.null(low_mean_and_std)) {
        mean_dropped_markers = logical(nrow(dat))
    } else {
        mean_dropped_markers = ((apply(dat, 1, 'mean', na.rm = TRUE) <
        low_mean_and_std) & (apply(dat, 1, 'sd',
        na.rm = TRUE) < low_mean_and_std))
    }
    if (is.null(q_low_var)) {
        vari_dropped_markers = logical(nrow(dat))
    } else {
        vari_dropped_markers = (apply(dat, 1, 'var', na.rm = TRUE) <
        quantile(apply(dat, 1, 'var', na.rm = TRUE), q_low_var))
    }
    dropped_markers = mean_dropped_markers | vari_dropped_markers
    return(list(dat[!dropped_markers,], rownames(dat)[dropped_markers]))
}
#' @title Draw densities
#' @description Draw column densities of an object over multiple plots by using
#' limma::plotDensities() function.
#' @param dat an object of log2-normalized protein (or gene) expressions,
#' containing markers in rows and samples in columns.
#' @param name name tag for the output file.
#' @param per.plot number of densities to be drawn on a single plot. If NULL,
#' ncol(object) will be used.
#' @param main character string, an overall title for the plot.
#' @param group vector or factor classifying the arrays into groups. Should be
#' same length as ncol(object).
#' @param legend character string giving position to place legend. See `legend`
#' for possible values. Can also be logical, with FALSE meaning no legend.
#' @return pdf plot(s).
#' @examples
#' dat = setNames(as.data.frame(matrix(1:(5*10),5,10),
#' row.names = paste('marker',1:5,sep='')), paste('sample',1:10,sep=''))
#' plotDen(dat, name = 'myresults')
plotDen = function(dat, name = '', per.plot = 8, main = NULL, group = NULL,
    legend = TRUE){
    if (is.null(per.plot)) {per.plot = ncol(dat)}
    num.plot = ceiling(ncol(dat)/per.plot)
    for (i in seq_len(num.plot)){
        pdf(paste(name,'_samples_',(per.plot*(i-1)+1),'to',min((per.plot*i),
        ncol(dat)),'.pdf', sep = ''), useDingbats = F)
        limma::plotDensities((dat[,(per.plot*(i-1)+1):min((per.plot*i),
        ncol(dat))]), main = main, group = group, legend = legend)
        dev.off()
    }
}
#' @title Analyze putative outliers
#' @description Calculates a statistical measure of each data entry being a
#' putative outlier
#' @param dat an object of log2-normalized protein (or gene) expressions,
#' containing markers in rows and samples in columns.
#' @return outlier p-statistics
#' @examples
#' dat = setNames(as.data.frame(matrix(1:(5*10),5,10),
#' row.names = paste('marker',1:5,sep='')), paste('sample',1:10,sep=''))
#' result = outScores(dat)
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
#' @title Artificially miss and impute each data entry individually by ignoring
#' outlying values
#' @description Infers the normal-state expression of a marker based on its
#' co-expression network, i.e., the weighted average of the marker's nearest
#' neighbors in the data. The returned imputed data will later be used to
#' elucidate dysregulated (protruding) events.
#' @param dat an object of log2-normalized protein (or gene) expressions,
#' containing markers in rows and samples in columns.
#' @param ku an integer in [1,num.markers], upper bound on the number of
#' nearest neighbors of a marker.
#' @param marker.proc.list character array, the row names of the data to be
#' processed/imputed.
#' @param miss.pstat the score threshold for ignoring potential outliers
#' during imputation. miss.pstat = 1 ignores values outside of the density box
#' (i.e., 1st-3rd quartiles). The algorithm ignores values lying at least
#' (1/miss.pstat)-1 times IQR away from the box; e.g., use miss.pstat=1 to
#' ignore all values lying outside of the box; use miss.pstat=0.4 to ignore
#' values lying at least 1.5 x IQR away from the box; use miss.pstat=0 to
#' employ all data during imputation.
#' @param verbose logical, to show progress of the algorithm.
#' @return the imputed data that putatively represents the expressions of the
#' markers in the (matched) normal states.
#' @examples
#' dat = setNames(as.data.frame(matrix(1:(5*10),5,10),
#' row.names = paste('marker',1:5,sep='')), paste('sample',1:10,sep=''))
#' imputed = artImpute(dat, ku = 2)
artImpute = function(dat, ku = 6, marker.proc.list = NULL, miss.pstat = 4E-1,
    verbose = FALSE) {
    out.pstats = outScores(dat)
    m = nrow(dat); n = ncol(dat)
    dat.imp = data.frame(matrix(nrow = nrow(dat), ncol = ncol(dat)))
    rownames(dat.imp) = rownames(dat); colnames(dat.imp) = colnames(dat)
    if (!methods::is(dat, 'matrix')) {dat = as.matrix(dat)}
    if (is.null(marker.proc.list)){marker.proc.list = seq_len(nrow(dat))}
    k = min(ku, nrow(dat)); dat.mis.ref = dat
    # remove all outlying expressions to not skew imputation
    dat.mis.ref[out.pstats < miss.pstat] = NA
    dat.dis = parallelDist::parDist(x=dat.mis.ref, method='euclidean',
                                    diag=TRUE, upper=TRUE)
    dat.dis = as.matrix(dat.dis)
    dat.cor = as.matrix(cor(t(dat.mis.ref), use = 'pairwise.complete.obs'))
    for (i in marker.proc.list) {
        # choose k euclidean neighbors
        sorted = sort(as.numeric(dat.dis[i,]), decreasing = FALSE,
            index.return = TRUE, na.last = TRUE)
        nei.ix = sorted$ix[2:(k+1)]
        # discard uncorrelated neighbors
        nei.ix = nei.ix[dat.cor[i,nei.ix]>.25]
        dat.knn = dat.mis.ref[nei.ix,]
        # weighted average to impute missing data
        weights = (1/sorted$x[nei.ix])**2
        if (length(weights)==1) {
            dat.imp[i,] = dat.knn #** weights
        } else {
            for (j in seq_len(ncol(dat))) {
                dat.imp[i,j] = stats::weighted.mean(dat.knn[
                    !is.na(dat.knn[,j]),j], weights[!is.na(dat.knn[,j])])
            }
        }
        # replace remaining NAs with the nearest complete values
        for (j in which(is.na(dat.imp[i,]))) {
            dat.imp[i,j] = dat.mis.ref[sorted$ix[which(!is.na(dat.mis.ref[
                sorted$ix[2:n],j]))[1]+1],j]
        }
        if (verbose){
            if (i %in% round(seq(round(length(marker.proc.list)/100),
                length(marker.proc.list), length = 100))){
                message(paste('Background inference: ',
                round(100*i/length(marker.proc.list)),'% done.'))
            }
        }
    }
    return(as.matrix(dat.imp))
}
#' @title Analyze dysregulated (protruding) events
#' @description For each marker processed, draws a scatter plot of matching
#' values of observed vs imputed expressions.
#' @param dat an object of log2-normalized protein (or gene) expressions,
#' containing markers in rows and samples in columns.
#' @param dat.imp the imputed data that putatively represents the expressions
#' of the markers in the (matched) normal states.
#' @param marker.proc.list character array, the row names of the data to be
#' processed for dysregulation.
#' @param verbose logical, to show progress of the algorithm
#' @return samples' distances to regression line (i.e., dysregulation) on the
#' scatter plots.
#' @return the scatter plots.
#' @examples
#' dat = setNames(as.data.frame(matrix(1:(5*10),5,10),
#' row.names = paste('marker',1:5,sep='')), paste('sample',1:10,sep=''))
#' dat.imp = artImpute(dat, ku=2)
#' result = dysReg(dat, dat.imp)
dysReg = function(dat, dat.imp, marker.proc.list = NULL, verbose = FALSE){
    if (is.null(marker.proc.list)) {marker.proc.list = rownames(dat)}
    m = nrow(dat); n = ncol(dat)
    num.omit.fit = round(.1*n)
    dat.dys = data.frame(matrix(NA, m, n))
    rownames(dat.dys) = rownames(dat); colnames(dat.dys) = colnames(dat)
    plot.list = list()
    for (i in which(rownames(dat) %in% marker.proc.list)) {
        out = gqplot(y = dat[i,], x = dat.imp[i,], omit.fit = num.omit.fit,
            ylab = 'Observed', xlab = 'Imputed')
        d=out[[1]]; plot.it=out[[2]]
        plot.list[[i]] = plot.it
        dat.dys[i, as.character(d$sampleID)] = d$dist2reg
        if (verbose){
            if (i %in% round(seq(round(m/100), m, length = 100))){
            message(paste('Dysregulation analysis:', round(100*i/m),'% done.'))
            }
        }
    }
    return(list(dat.dys, plot.list))
}
#' @title Analyze dysregulation significance
#' @description Rank-order markers by the significance of deviation of the
#' observed expressions from the (matched) imputed expressions based on the
#' Kolmogorov-Smirnov (KS) test.
#' @param dat an object of log2-normalized protein (or gene) expressions,
#' containing markers in rows and samples in columns.
#' @param dat.imp the imputed data that putatively represents the expressions
#' of the markers in the (matched) normal states.
#' @param marker.proc.list character array, the row names of the data to be
#' processed for dysregulation significance.
#' @param pval.insig p-value threshold to determine spurious (null)
#' dysregulation events.
#' @return each marker's p-value of the statistical significance between its
#' observed vs imputed values computed by the KS test.
#' @return ranked p-values (KS test) of the significant markers, which are
#' lower than pval.insig.
#' @return ranked significantly dysregulated markers with p-values lower than
#' pval.insig.
#' @return ranked p-values (KS test) of the insignificant markers, which are
#' greater than pval.insig.
#' @return ranked insignificantly dysregulated markers (spurious
#' dysregulations) with p-values greater than pval.insig.
#' @examples
#' set.seed(1)
#' dat = setNames(as.data.frame(matrix(runif(10*10),10,10),
#' row.names = paste('marker',1:10,sep='')), paste('sample',1:10,sep=''))
#' dat.imp = artImpute(dat, ku=6)
#' result = statTest(dat, dat.imp) # the dysregulations on marker4 is
#' # statistically significant with p-value 0.05244755.
statTest = function(dat, dat.imp, marker.proc.list = NULL, pval.insig = 2E-1) {
    if (is.null(marker.proc.list)) {marker.proc.list = rownames(dat)}
    m = nrow(dat); n = ncol(dat)
    dat.imp.test = matrix(NA,m,1); rownames(dat.imp.test) = rownames(dat)
    for (i in marker.proc.list) {dat.imp.test[i,] = stats::ks.test(t(dat[i,]),
        dat.imp[i,], alternative = 'two.sided')$p.value}
    dat.imp.test.sor = sort(dat.imp.test, decreasing = FALSE,
        index.return = TRUE, na.last = TRUE)
    markers.imp.insig.loc =
        dat.imp.test.sor$ix[(dat.imp.test.sor$x>pval.insig) &
        !is.na(dat.imp.test.sor$x)]
    markers.imp.sig.loc = dat.imp.test.sor$ix[(dat.imp.test.sor$x<pval.insig) &
        !is.na(dat.imp.test.sor$x)]
    if (length(markers.imp.insig.loc) < 1)
        {markers.imp.insig.loc=dat.imp.test.sor$ix[!is.na(dat.imp.test.sor$x)]
        markers.imp.insig.loc = markers.imp.insig.loc[
            length(markers.imp.insig.loc)]}
    if (length(markers.imp.sig.loc) < 1)
        {markers.imp.sig.loc = dat.imp.test.sor$ix[1]}
    markers.imp.insig = rownames(dat)[rev(markers.imp.insig.loc)]
    markers.imp.sig = rownames(dat)[markers.imp.sig.loc]
    dat.imp.test.sig = data.frame(pval = dat.imp.test[markers.imp.sig,])
    dat.imp.test.insig = data.frame(pval = dat.imp.test[markers.imp.insig,])
    return(list(dat.imp.test, dat.imp.test.sig, markers.imp.sig,
        dat.imp.test.insig, markers.imp.insig))
}
#' @title Display outlying expressions
#' @description Mark outlying expressions on the scatter plot of a given marker
#' @param marker.proc.list character array, the row names of the data to be
#' processed for outlier analyses and for plotting.
#' @param dat an object of log2-normalized protein (or gene) expressions,
#' containing markers in rows and samples in columns.
#' @param dat.imp the imputed data that putatively represents the expressions
#' of the markers in the (matched) normal states.
#' @param dat.imp.test marker's p-value of the statistical significance between
#' its observed vs imputed values computed by the Kolmogorov-Smirnov test.
#' @param dat.dys samples' distances to regression line (i.e., dysregulation)
#' on the scatter plots.
#' @param dys.sig.thr.upp the dysregulation score threshold to elucidate/mark
#' significantly dysregulated outlier events.
#' @param dataset the cohort name to be used in the output files.
#' @param num.omit.fit number of outlying events to ignore when fitting a
#' marker's observed expressions to the imputed ones.
#' @param draw.sc logical, to draw a scatter plot for every marker in
#' marker.proc.list in a separate PDF file.
#' @param draw.vi logical, to draw a violin plot for every marker in
#' marker.proc.list in a separate PDF file.
#' @param conf.int confidence interval to display around the regression line
#' @param ylab a title for the y axis
#' @param xlab a title for the x axis
#' @return the scatter plots of the markers where the outlier dysregulation
#' events are highlighted by red mark.
#' @examples
#' set.seed(1)
#' dat = setNames(as.data.frame(matrix(runif(10*10),10,10),
#' row.names = paste('marker',1:10,sep='')), paste('sample',1:10,sep=''))
#' dat.imp = artImpute(dat, ku=6)
#' dat.imp.test = statTest(dat, dat.imp)[[1]]
#' dat.dys = dysReg(dat, dat.imp)[[1]]
#' plots = markOut(dat, dat.imp, dat.imp.test, dat.dys, dys.sig.thr.upp = .25)
markOut = function(dat, dat.imp, dat.imp.test, dat.dys, dys.sig.thr.upp,
    marker.proc.list = NULL, dataset = '', num.omit.fit = NULL,
    draw.sc = TRUE, draw.vi = TRUE, conf.int = .95, ylab = 'Observed',
    xlab = 'Inferred'){
    if (is.null(marker.proc.list)) {marker.proc.list = rownames(dat)}
    if (is.null(num.omit.fit)) {num.omit.fit = round(.1*ncol(dat))}
    plot.list.marked = list()
    for (marker in marker.proc.list) {
        marker.loc = which(rownames(dat)==marker)
        # significant outlying events for the given marker
        marker.out.exp.loc = which(dat.dys[marker,]>dys.sig.thr.upp)
        # | dat.dys[marker,]<dys.sig.thr.low)
        marker.out.samples = colnames(dat.dys)[marker.out.exp.loc]
        dat.dys[marker,marker.out.exp.loc]
        out = gqplot(y = dat[marker,], x = dat.imp[marker,], ci = conf.int,
            ylab = ylab, xlab = xlab,
            highlight = marker.out.samples, omit.fit = num.omit.fit)
            #, minl = minl, maxl = maxl)
        d=out[[1]]; plot.it=out[[2]]
        plot.list.marked[[marker.loc]] = plot.it;
        # draw scatter plot
        if (draw.sc) {
            # pdf(paste(dataset,'.',marker,'.sc','.dysreg.pdf',sep=''),
            # width=2.1, height=2.3, useDingbats = F)
            # print(plot.list.marked[[marker.loc]]); dev.off()
            ggsave(filename = paste(dataset,'.',marker,'.sc','.dysreg.pdf'
                ,sep=''), plot = plot.list.marked[[marker.loc]],
                device = 'pdf', width = 2.1, height = 2.3,
                dpi = 'print', useDingbats = F)
        }
        # draw violin plot
        if (draw.vi) {
            df = data.frame(Observed = as.numeric(dat[marker,]),
                Imputed = as.numeric(dat.imp[marker,]))
            df = reshape::melt(df); colnames(df) = c('data','marker')
            pl = ggplot2::ggplot(df, ggplot2::aes(x=data, y=marker)) +
                ggplot2::geom_violin(trim = FALSE) +
                ggplot2::geom_boxplot(width=0.1) +
                ggplot2::xlab('') +
                ggplot2::ylab(marker) +
                ggplot2::ggtitle(as.double(dat.imp.test[marker,]))
            # # ggplot2::ylim(minl,maxl)
            pdf(paste(dataset,'.',marker,'.vi','.dysreg.pdf',sep=''),width=3,
                height=3, useDingbats = F)
            print(pl)
            dev.off()
        }
    }
    return(plot.list.marked)
}
#' @title Rank markers by the percentage of outlying events
#' @description Ranks markers in the order of decreasing percentage of outlying
#' events.
#' @param dat.dys samples' distances to regression line (i.e., dysregulation)
#' on the scatter plots.
#' @param marker.proc.list character array, the row names of the data to be
#' processed for outlier analyses.
#' @param dys.sig.thr.upp the dysregulation score threshold to elucidate/mark
#' significantly dysregulated outlier events.
#' @return markers rank-ordered by the percentage of outliers over the samples.
#' @return the percentages of outliers corresponding to ranked markers.
#' @examples
#' set.seed(1)
#' dat = setNames(as.data.frame(matrix(runif(10*10),10,10),
#' row.names = paste('marker',1:10,sep='')), paste('sample',1:10,sep=''))
#' dat.imp = artImpute(dat, ku=6)
#' dat.dys = dysReg(dat, dat.imp)[[1]]
#' result = rankPerOut(dat.dys, dys.sig.thr.upp = .25)
rankPerOut = function(dat.dys, marker.proc.list = NULL, dys.sig.thr.upp){
    if (is.null(marker.proc.list)) {marker.proc.list = rownames(dat.dys)}
    m = dim(dat.dys)[1]
    marker.out.exp.per = data.frame(percent.dysregulated=array(NA,dim=c(m,1)))
    rownames(marker.out.exp.per) = rownames(dat.dys)
    for (marker in marker.proc.list) {
        marker.out.exp.per[rownames(dat.dys)==marker,] =
            100*rowSums(dat.dys[marker,]>dys.sig.thr.upp, na.rm=TRUE)/
            rowSums(!is.na(dat.dys[marker,]>dys.sig.thr.upp))
    }
    marker.out.exp.per.sor = sort(marker.out.exp.per[,1], decreasing = TRUE,
        index.return = TRUE, na.last = TRUE)
    marker.out.exp.per.sor.names = rownames(marker.out.exp.per)[
        marker.out.exp.per.sor$ix[seq_along(marker.proc.list)]]
    return(list(marker.out.exp.per.sor.names, marker.out.exp.per))
}
#' @title Hierarchical cluster analysis
#' @description Displays the hierarchically clustered data by the "pheatmap"
#' package.
#' The numbers of clusters along the markers/samples can be set by the user,
#' then the cluster structures are estimated by pair-wise analysis.
#' @param data an object of log2-normalized protein (or gene) expressions,
#' containing markers in rows and samples in columns.
#' @param annotation_row data frame that specifies the annotations shown on
#' left side of the heat map. Each row defines the features for a specific
#' row. The rows in the data and in the annotation are matched using
#' corresponding row names. Note that color schemes takes into account if
#' variable is continuous or discrete.
#' @param annotation_col similar to annotation_row, but for columns.
#' @param annotation_colors list for specifying annotation_row and
#' annotation_col track colors manually. It is possible to define the colors
#' for only some of the features.
#' @param main character string, an overall title for the plot.
#' @param legend logical, to determine if legend should be drawn or not.
#' @param clustering_distance_rows distance measure used in clustering rows.
#' Possible values are "correlation" for Pearson correlation and all the
#' distances supported by dist, such as "euclidean", etc. If the value is
#' none of the above it is assumed that a distance matrix is provided.
#' @param clustering_distance_cols distance measure used in clustering
#' columns. Possible values the same as for clustering_distance_rows.
#' @param display_numbers logical, determining if the numeric values are also
#' printed to the cells. If this is a matrix (with same dimensions as original
#' matrix), the contents of the matrix are shown instead of original values.
#' @param number_format format strings (C printf style) of the numbers shown in
#'  cells. For example "%.2f" shows 2 decimal places and "%.1e" shows
#'  exponential notation (see more in sprintf).
#' @param num_clusters_row number of clusters the rows are divided into, based
#' on the hierarchical clustering (using cutree), if rows are not clustered,
#' the argument is ignored.
#' @param num_clusters_col similar to num_clusters_row, but for columns.
#' @param cluster_rows logical, determining if the rows should be clustered;
#' or a hclust object.
#' @param cluster_cols similar to cluster_rows, but for columns.
#' @param border_color color of cell borders on heatmap, use NA if no border
#' should be drawn.
#' @param annotate_new_clusters_col logical, to annotate cluster IDs (column)
#' that will be identified.
#' @param zero_white logical, to display 0 values as white in the colormap.
#' @param color_low, color code for the low intensity values in the colormap.
#' @param color_mid, color code for the medium intensity values in the colormap.
#' @param color_high, color code for the high intensity values in the colormap.
#' @param color_palette vector of colors used in heatmap.
#' @param show_rownames boolean, specifying if row names are be shown.
#' @param show_colnames boolean, specifying if column names are be shown.
#' @param min_data numeric, data value corresponding to minimum intensity in
#' the color_palette
#' @param max_data numeric, data value corresponding to maximum intensity in
#' the color_palette
#' @param treeheight_row the height of a tree for rows, if these are clustered.
#' Default value is 50 points.
#' @param treeheight_col the height of a tree for columns, if these are
#' clustered. Default value is 50 points.
#' @return tree, the hierarchical tree structure.
#' @return cluster_IDs_row, the (row) cluster identities of the markers.
#' @return cluster_IDs_col, the (column) cluster identities of the samples.
#' @examples
#' set.seed(1)
#' dat = setNames(as.data.frame(matrix(runif(10*10),10,10),
#' row.names = paste('marker',1:10,sep='')), paste('sample',1:10,sep=''))
#' result = clusterData(dat)
clusterData = function(data, annotation_row = NULL, annotation_col = NULL,
    annotation_colors = NULL, main = NA, legend = TRUE,
    clustering_distance_rows = 'euclidean',
    clustering_distance_cols = 'euclidean', display_numbers = FALSE,
    number_format = '%.0f', num_clusters_row = NULL, num_clusters_col = NULL,
    cluster_rows = TRUE, cluster_cols = TRUE, border_color = 'gray60',
    annotate_new_clusters_col = FALSE, zero_white = FALSE,
    color_low = '#006699', color_mid = 'white', color_high = 'red',
    color_palette = NULL, show_rownames = FALSE, show_colnames = FALSE,
    min_data = min(data, na.rm=TRUE), max_data = max(data, na.rm=TRUE),
    treeheight_row =
        ifelse(methods::is(cluster_rows, "hclust") || cluster_rows, 50, 0),
    treeheight_col =
        ifelse(methods::is(cluster_cols, "hclust") || cluster_cols, 50, 0)){
    if (is.null(num_clusters_col)) {num_clusters_col = 1}
    if (is.null(num_clusters_row)) {num_clusters_row = 1}
    if (zero_white) {paletteLength = 100
        my.color = grDevices::colorRampPalette(
            c(color_low, color_mid, color_high))(paletteLength)
        my.breaks = c(seq(min_data, 0, length.out = ceiling(paletteLength/2) ),
            seq(max_data/paletteLength, max_data,
            length.out = floor(paletteLength/2)))
        uni.bre = uniq(my.breaks, index.return = TRUE)
        my.breaks = uni.bre$x
        my.color = my.color[uni.bre$ix]
    } else {
        if (is.null(color_palette)){
            my.color=grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(
                n = 7, name = "RdYlBu")))(100)
        } else {
            my.color=grDevices::colorRampPalette(RColorBrewer::brewer.pal(
                n = 7, name = color_palette))(100)
        }
        my.breaks = NA
    }
    tree = pheatmap::pheatmap(data, annotation_col = annotation_col,
        annotation_colors = annotation_colors,
        display_numbers = display_numbers, number_format = number_format,
        clustering_distance_cols = clustering_distance_cols,
        clustering_distance_rows = clustering_distance_rows,
        cluster_cols = cluster_cols, cluster_rows = cluster_rows,
        cutree_cols = num_clusters_col, show_colnames = show_colnames,
        cutree_rows = num_clusters_row, show_rownames = show_rownames,
        main = main, color = my.color, breaks = my.breaks,
        border_color = border_color, treeheight_row = treeheight_row,
        treeheight_col = treeheight_col, legend = legend)
    if (cluster_cols) {cluster_IDs_col = cutree(tree$tree_col,
        k = num_clusters_col)} else {cluster_IDs_col = NULL}
    if (cluster_rows) {cluster_IDs_row = cutree(tree$tree_row,
        k = num_clusters_row)} else {cluster_IDs_row = NULL}
    if (annotate_new_clusters_col) {
        annotation_col_new = cbind(annotation_col,data.frame(cluster_IDs_col)[
            rownames(annotation_col),])
        colnames(annotation_col_new)[length(annotation_col_new)] =
            'New.subtype'
        tree = pheatmap::pheatmap(data, annotation_col = annotation_col_new,
            annotation_colors = annotation_colors,
            display_numbers = display_numbers, number_format = number_format,
            clustering_distance_cols = clustering_distance_cols,
            clustering_distance_rows = clustering_distance_rows,
            cluster_cols = cluster_cols, cluster_rows = cluster_rows,
            cutree_cols = num_clusters_col, show_colnames = show_colnames,
            cutree_rows = num_clusters_row, show_rownames = show_rownames,
            main = main, color = my.color, breaks = my.breaks,
            border_color = border_color, treeheight_row = treeheight_row,
            treeheight_col = treeheight_col, legend = legend)
    }
    return(list(tree, cluster_IDs_row, cluster_IDs_col))
}
#' @title Outlier protein and phosphosite target identification
#' @description Find outlying markers and events across cancer types.
#' @param data a list object where each element contains a proteomics data for
#' a different cohort (markers in the rows, samples in the columns) or a
#' character string defining the path to such data (in .RDS format).
#' @param mad.norm logical, to normalize the proteomes to have a unit Median
#' Absolute Deviation.
#' @param cohort.names character array.
#' @param panel a character string describing marker panel, e.g., 'kinases'.
#' Use 'global' to analyze all markers quantified across cohorts (default).
#' Use 'pancan' to analyze the markers commonly quantified across the cohorts.
#' @param panel.markers a character array containing the set of marker names
#' that user wants to analyze, e.g., panel.markers = c("AAK1", "AATK", "ABL1",
#' "ABL2", ...).
#' @param tol.nas a constant in [0,100], tolerance for the percentage of NAs
#' in a marker, e.g., tol.nas = 20 will filter out markers containing 20\% or
#' more NAs across samples.
#' @param ku an integer in [1,num.markers], upper bound on the number of
#' nearest neighbors of a marker.
#' @param miss.pstat a constant in [0,1], statistic to estimate potential
#' outliers. See `artImpute()`.
#' @param demo.panels logical, to draw demographics of the panel in each
#' cohort.
#' @param save.data logical, to save intermediate data (background inference
#' and dysregulation measures).
#' @param draw.sc.plots logical, to draw each marker's qqplot of observed vs
#' inferred (imputed) expressions.
#' @param draw.vi.plots logical, to draw each marker's violin plot of observed
#' vs imputed expressions.
#' @param draw.sc.markers character array, marker list to draw scatter plots
#' @param draw.ou.plots logical, to draw each marker's outlier prevalence
#' (by the percentage of outlying samples) across the cohorts.
#' @param draw.ou.markers character array, marker list to draw pan-cancer
#' outlier percentage plots
#' @param verbose logical, to show progress of the algorithm.
#' @return dysregulation scores of every marker for each sample.
#' @return the imputed data that putatively represents the expressions of the
#' markers in the (matched) normal states.
#' @return the result of Kolmogorov-Smirnov tests that evaluates the
#' statistical significance of each marker's outlier samples.
#' @return a data list containing, for each cohort, the percentage of outlier
#' samples for every marker.
#' @return a data list containing, for each cohort, the outlier significance
#' threshold.
#' @examples
#' set.seed(1)
#' dat = setNames(as.data.frame(matrix(runif(10*10),10,10),
#' row.names = paste('marker',1:10,sep='')), paste('sample',1:10,sep=''))
#' result = oppti(dat)
#' @seealso [artImpute()] for how to set `miss.pstat` and `ku`
oppti = function(data, mad.norm = FALSE, cohort.names = NULL, panel = 'global',
    panel.markers = NULL, tol.nas = 20, ku = 6, miss.pstat = .4,
    demo.panels = FALSE, save.data = FALSE, draw.sc.plots = FALSE,
    draw.vi.plots = FALSE, draw.sc.markers = NULL,
    draw.ou.plots = FALSE, draw.ou.markers = NULL, verbose = FALSE) {
    # fix flags: let draw.*.markers master draw.*.plots
    if (!is.null(draw.sc.markers)) {draw.sc.plots = TRUE}
    if (!is.null(draw.ou.markers)) {draw.ou.plots = TRUE}
    # Load data
    if (is.character(data)){tryCatch({data = readRDS(data)},
        warning=function(w){message(w)}, finally={})} else if (!is.list(data)){
        warning('Unrecognized format for the "data" object; must be a "list",
        or a "character" string referring to the file path of such object.');
        return()}
    if (!methods::is(data, 'list')) {data = list(data)}
    if (mad.norm) {data = lapply(data, function(x){x=madNorm(x)})}
    if (is.null(cohort.names)) {cohort.names = unlist(lapply(data, function(x){
        x=strsplit(colnames(x)[1],'\\.')[[1]][1]}))} #cohort names
    pan.num = length(cohort.names) #number of cohorts
    tmp.lis = as.list(rep(NA,pan.num)) #template for a cohort-size data object
    # Filter out markers with >= tol.nas percent missing values
    pan.dat = lapply(data, function(x){x=dropMarkers(x,
        percent_NA = tol.nas/100, low_mean_and_std = NULL,
        q_low_var = NULL)[[1]]})
    # Choose samples (tumor or normal) for analyses #***
    pan.dat = lapply(pan.dat, function(x){if (any(regexpr('Tumor',
        colnames(x))>0)){x=x[,regexpr('Tumor', colnames(x)) > 0]} else {x=x}})
    pan.mar = rownames(pan.dat[[1]]); if (pan.num>1) {for (i in 2:pan.num) {
        pan.mar = intersect(pan.mar, rownames(pan.dat[[i]]))}}
        #determine markers quantified pan-cancer
    if (is.null(panel.markers)) {
        switch(panel, 'global' = (panel.markers =
            unique(unlist(lapply(data,function(x) {x=rownames(x)})))),
            'pancan' = (panel.markers = pan.mar))
    }
    # Define which markers to process in each cancer
    pan.proc.markers = lapply(pan.dat, function(x){x=panel.markers[
        panel.markers %in% rownames(x)]})
        # Panel markers detected (and complete) in each cancer
    # Observation
    if (demo.panels) {
        pan.det = lapply(data, function(x)
            {x=sum(panel.markers %in% rownames(x))})
        pan.com = lapply(pan.dat, function(x)
            {x=sum(panel.markers %in% rownames(x))})
        pan.nas = lapply(data, function(x){x=100*rowSums(is.na(x))/ncol(x)})
        for (i in seq_len(pan.num))
            {pdf(paste('',cohort.names[i],'.',panel,'.pdf',sep=''), width=6,
                height=4, useDingbats = F);
            graphics::hist(pan.nas[[i]][rownames(data[[i]]) %in% panel.markers]
            , xlab = '% NAs', ylab = 'Number of markers'
            , main = paste(cohort.names[i],' ', panel, ' > ', pan.det[[i]]
            , ' detected > ', pan.com[[i]], ' complete (<', tol.nas,'% NAs)'
            , sep='')); dev.off()}
        for (i in seq_len(pan.num)) {plotDen(dat = data[[i]],
            name = cohort.names[i], per.plot = NULL, legend = FALSE)}
        plotDen(dat = cbindNA(data), group = cohort.names, name = 'data')
    }
    # Generate inferred expressions by weighted average of nearest neighbors
    tic = Sys.time()
    pan.dat.imp = tmp.lis; for (i in seq_len(pan.num)) {pan.dat.imp[[i]] =
        artImpute(dat = pan.dat[[i]], ku = ku, miss.pstat = miss.pstat,
        marker.proc.list = pan.proc.markers[[i]], verbose = verbose)}
    toc = Sys.time()-tic
    if (save.data)
        {saveRDS(pan.dat.imp, file=paste('pan.dat.imp.',panel,'.RDS',sep=''))}
    # Analyze dysregulation events by a linear fit (qqplot)
    pan.dat.dys = tmp.lis; for (i in seq_len(pan.num)) {pan.dat.dys[[i]] =
        dysReg(pan.dat[[i]], pan.dat.imp[[i]], pan.proc.markers[[i]],
        verbose = verbose)[[1]]}
    if (save.data)
        {saveRDS(pan.dat.dys, file=paste('pan.dat.dys.',panel,'.RDS',sep=''))}
    # Analyze spurious events
    message('Analyzing non-dysregulated markers [statTest]...')
    pan.dat.imp.test = tmp.lis; pan.markers.imp.insig = tmp.lis
    for (i in seq_len(pan.num))
        {out = statTest(pan.dat[[i]], pan.dat.imp[[i]], pan.proc.markers[[i]]);
        pan.dat.imp.test[[i]]=out[[1]]; pan.markers.imp.insig[[i]]=out[[5]]}
    # Define a significance value for outlier samples in contrast to null
        # expressions (x > 95% of null exp. --> pval < 0.05)
    pan.dat.imp.insig.all.dys = tmp.lis; for (i in seq_len(pan.num))
        {pan.dat.imp.insig.all.dys[[i]] =
        unlist(pan.dat.dys[[i]][pan.markers.imp.insig[[i]],])}
    pan.dys.sig.thr.upp = lapply(pan.dat.imp.insig.all.dys, function(x)
        {x=quantile(x, .95, na.rm = TRUE)})
    # Permutation test to associate FDR of the marker overexpressions
    message('Running permutation tests to associate FDR for each marker...')
    pan.sym.tes = tmp.lis;
    for (i in seq_len(pan.num)) {
        message(paste0('Running permutation tests for ', cohort.names[i]))
        dat.ids = colnames(pan.dat[[i]])
        imp.ids = colnames(pan.dat.imp[[i]])
        mar.sym.tes = data.frame(p = array(NA, nrow(pan.dat[[i]])), row.names = rownames(pan.dat[[i]]))
        for (j in which(rownames(pan.dat[[i]]) %in% pan.proc.markers[[i]])) {
            df = data.frame(ids = c(dat.ids, imp.ids),
                exp = c(rep('observed', length(dat.ids)),
                rep('imputed', length(imp.ids))),
                val = as.numeric(c(pan.dat[[i]][j,], pan.dat.imp[[i]][j,])))
            df = df[!df$ids %in% df$ids[is.na(df$val)],]
            mar.sym.tes$p[j] = coin::pvalue(coin::symmetry_test(val ~ exp,
                data = df, alternative = 'two.sided', paired = T))
        }
        mar.sym.tes$FDR = p.adjust(mar.sym.tes$p, method = 'BH')
        pan.sym.tes[[i]] = mar.sym.tes
    }

    if (demo.panels) {
        message('Generating demo for the panel markers...')
        colors = c('red','orange','yellow','green','blue','purple');
        limx=1; limy=1; pdf('pan.null.dys.ecdf.pdf', width=6,height=6,
            useDingbats = F);
        graphics::plot(ecdf(pan.dat.imp.insig.all.dys[[1]]),
        col=colors[1], xlim=c(-limx,limx), ylim=c(1-limy,limy),
        main = 'Empirical CDF')
        if (pan.num>1) {for (i in 2:pan.num)
            {graphics::lines(ecdf(pan.dat.imp.insig.all.dys[[i]]),
            col=colors[i], xlim=c(-limx,limx), ylim=c(1-limy,limy))}}
        graphics::legend(-limx, limy, legend=cohort.names,
            col=colors, lty=rep(1,pan.num), cex=.8);
        dev.off()
    }
    # Mark outlying expressions of a given marker on its qqplot, draw
    # violin plots
    if (draw.sc.plots | draw.vi.plots) {for (i in seq_len(pan.num)){
        if (is.null(draw.sc.markers)) {draw.sc.markers.i=pan.proc.markers[[i]]
        } else {draw.sc.markers.i = draw.sc.markers[draw.sc.markers %in%
            pan.proc.markers[[i]]]}
        if (!is.null(draw.sc.markers.i)) {
            # message(paste(c(pan.num[[i]], '|', draw.sc.markers.i),
            # collapse = ' '))
            message('Drawing scatter plots [markOut]...')
            markOut(pan.dat[[i]], pan.dat.imp[[i]], pan.dat.imp.test[[i]],
                pan.dat.dys[[i]], pan.dys.sig.thr.upp[[i]], draw.sc.markers.i,
                cohort.names[i],draw.sc=draw.sc.plots,draw.vi=draw.vi.plots)}}}
    # Rank markers by the percentage of outlying events
    message('Building heatmaps for percentage of outliers across cancers [rankPerOut] ...')
    pan.marker.out.exp.per = tmp.lis; for (i in seq_len(pan.num))
        {pan.marker.out.exp.per[[i]] = rankPerOut(pan.dat.dys[[i]],
        pan.proc.markers[[i]], pan.dys.sig.thr.upp[[i]])[[2]]}
    if (draw.ou.plots) {
        # Draw markers' percentages of outlying events for each cancer
        if (is.null(draw.ou.markers)) {draw.ou.markers=panel.markers
        } else {draw.ou.markers =
            draw.ou.markers[draw.ou.markers %in% panel.markers]}
        pan.mar.out.exp.per =
            data.frame(matrix(NA,length(panel.markers),pan.num))
        colnames(pan.mar.out.exp.per) = cohort.names
        rownames(pan.mar.out.exp.per) = panel.markers
        for (i in seq_len(pan.num))
            {pan.mar.out.exp.per[rownames(pan.marker.out.exp.per[[i]]),i] =
                pan.marker.out.exp.per[[i]]}
        # Display the predefined marker set
        pan.mar.out.exp.per = pan.mar.out.exp.per[which(rownames(
            pan.mar.out.exp.per) %in% draw.ou.markers),]
        # Pan-cancer highly-outlying markers
        pan.mar.out.exp.per.rat.sor = sort(rowSums(pan.mar.out.exp.per,
            na.rm = TRUE)/ncol(pan.mar.out.exp.per), decreasing = TRUE,
            index.return = TRUE, na.last = TRUE)
        tmp = as.data.frame(pan.mar.out.exp.per[
            pan.mar.out.exp.per.rat.sor$ix[
                pan.mar.out.exp.per.rat.sor$x>0],])
        rownames(tmp) = rownames(pan.mar.out.exp.per)[
            pan.mar.out.exp.per.rat.sor$ix[
                pan.mar.out.exp.per.rat.sor$x>0]]
        colnames(tmp) = colnames(pan.mar.out.exp.per)
        pan.mar.ranked.out.exp.per.tree = clusterData(tmp,
            cluster_cols = FALSE, cluster_rows = FALSE, color_palette = 'Reds',
            display_numbers = FALSE, main = '% of outliers')[[1]]
        pdf(paste('pan.cancer.',panel,'.markers.outlier.scores.pdf', sep = ''),
            width = max(2,ceiling((ncol(tmp))**(.7)-1)),
            height = max(2,ceiling((nrow(tmp))**(.7)-2)), useDingbats = F)
        print(pan.mar.ranked.out.exp.per.tree); dev.off()
        tmp = t(pan.mar.out.exp.per[pan.mar.out.exp.per.rat.sor$ix[
            pan.mar.out.exp.per.rat.sor$x>0][seq_len(min(20,length(
                pan.mar.out.exp.per.rat.sor$ix[
                    pan.mar.out.exp.per.rat.sor$x>0])))],])
        colnames(tmp) = rownames(pan.mar.out.exp.per)[
            pan.mar.out.exp.per.rat.sor$ix[
                pan.mar.out.exp.per.rat.sor$x>0][seq_len(min(20,length(
                    pan.mar.out.exp.per.rat.sor$ix[
                        pan.mar.out.exp.per.rat.sor$x>0])))]]
        rownames(tmp) = colnames(pan.mar.out.exp.per)
        # Display the predefined marker set
        message('Drawing heatmaps for percentage of outliers across cancers [rankPerOut] ...')
        pan.mar.ranked20.t.out.exp.per.tree = clusterData(tmp,
            cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE,
            main = '% of outliers', color_palette = 'Reds')[[1]]
        pdf(paste('pan.cancer.',panel,'.top20.highly.outlying.markers.pdf',
            sep = ''), width = max(2,ceiling((ncol(tmp))**(.7)-1)),
            height = max(2,ceiling((nrow(tmp))**(.7)-2)), useDingbats = F)
        print(pan.mar.ranked20.t.out.exp.per.tree); dev.off()
        # Pan-cancer top-20 variably-highly-outlying markers
        if (pan.num>1) {
            pan.mar.out.exp.per.sd.sor = sort(apply(pan.mar.out.exp.per, 1,
                'sd', na.rm = TRUE), decreasing = TRUE, index.return = TRUE,
                na.last = TRUE)
            tmp = t(pan.mar.out.exp.per[pan.mar.out.exp.per.sd.sor$ix[
                pan.mar.out.exp.per.sd.sor$x>0][seq_len(min(20,length(
                    pan.mar.out.exp.per.sd.sor$ix[
                        pan.mar.out.exp.per.sd.sor$x>0 & !is.na(
                            pan.mar.out.exp.per.sd.sor$x)])))],])
            pan.mar.ranked20.t.sd.out.exp.per.tree = clusterData(tmp,
                cluster_cols = FALSE, cluster_rows = FALSE,
                display_numbers = TRUE, main = '% of outliers',
                color_palette = 'Reds')[[1]]
            pdf(paste('pan.cancer.',panel,
                '.top20.variably.outlying.markers.pdf', sep = ''),
                width = max(2,ceiling((ncol(tmp))**(.7)-1)),
                height = max(2,ceiling((nrow(tmp))**(.7)-2)), useDingbats = F)
            print(pan.mar.ranked20.t.sd.out.exp.per.tree); dev.off()
        }
    }
    message('End of analysis.')
    if (pan.num>1){
        res = list(pan.dat.dys, pan.dat.imp, pan.dat.imp.test,
            pan.marker.out.exp.per, pan.dys.sig.thr.upp, pan.sym.tes)
    } else {
        res = list(pan.dat.dys[[1]], pan.dat.imp[[1]], pan.dat.imp.test[[1]],
            pan.marker.out.exp.per[[1]], pan.dys.sig.thr.upp[[1]], pan.sym.tes[[1]])
        for (i in seq_along(res)) {names(res[[i]]) = cohort.names}
    }
    return(res)
}
