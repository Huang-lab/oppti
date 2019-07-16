#' @title Filter out markers
#' @description Filters out markers based on the percantage of missing values,
#'  low-expression and low-variability rates.
#' @param dat, an object of log2-normalized gene (or protein) expressions,
#'  markers in rows and samples in columns.
#' @param percent_NA, a constant in [0,1], is the percentage of missing values
#'  that will be tolerated in the filtered data.
#' @param low_mean_and_std, a constant in [0,inf], is the lower-bound of the
#'  mean or standard deviation of a marker in the filtered data.
#' @param q_low_var, a constant in [0,1], is the quantile of marker variances
#'  which serves as a lower-bound of the marker variances in the filtered data.
#' @return data_filtered, with same class as the input data;
#'  dropped_marker_names, rownames (markers) of the data that are filtered out
#'  due to low-expression or low-variability.
dropMarkers = function(dat, percent_NA = .2, low_mean_and_std = .05,
    q_low_var = .25, force_drop = NULL){
    if (!is.null(force_drop)) {dat = dat[!(rownames(dat) %in% force_drop),]}
    dat = dat[rowSums(is.na(dat))/ncol(dat) < percent_NA,]
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
#'  limma::plotDensities() function.
#' @param dat, an object of log2-normalized gene (or protein) expressions,
#'  markers in rows and samples in columns.
#' @param name, name tag for the output file.
#' @param per.plot, number of densities to be drawn on a single plot. If NULL,
#'  ncol(object) will be used.
#' @param group, vector or factor classifying the arrays into groups. Should be
#'  same length as ncol(object).
#' @param legend, character string giving position to place legend. See legend
#'  for possible values. Can also be logical, with FALSE meaning no legend.
#' @return a (number) pdf plot(s)
plotDen = function(dat, name = '', per.plot = 8, main = NULL, group = NULL,
    legend = TRUE){
    if (is.null(per.plot)) {per.plot = ncol(dat)}
    num.plot = ceiling(ncol(dat)/per.plot)
    for (i in seq_len(num.plot)){
        pdf(paste(name,'_samples_',(per.plot*(i-1)+1),'to',min((per.plot*i),
        ncol(dat)),'.pdf', sep = ''))
        limma::plotDensities((dat[,(per.plot*(i-1)+1):min((per.plot*i),
        ncol(dat))]), main = main, group = group, legend = legend)
        dev.off()
    }
}
#' @title Artificially miss and impute each data entry individually by ignoring
#'  outlying values
#' @description Infers likely expressions of a marker in the (matched) normal
#'  cohort, based on the weigted average of the marker's nearest neighbors in
#'  the case cohort. The returned imputed data will later be used to elucidate
#'  protruding (dysregulated) events.
#' @param marker.proc.list, the rownames of the data to be processed/imputed.
#' @param miss.pstat, the score threshold for ignoring potential outliers
#'  during imputation. miss.pstat = 1 ignores values outside of the box (i.e.,
#'  1st-3rd quartiles).
#'  The algorithm ignores values lying at least (1/miss.pstat)-1 times IQR away
#'  from the box; e.g., use miss.pstat=1 to ignore all values lying outside of
#'  the box; use miss.pstat=0.4 to ignore values lying at least 1.5 x IQR away
#'  from the box; use miss.pstat=0 to employ all data during imputation.
#' @return dat.imp, the imputed data that likely represents the expressions of
#'  the markers in the (matched) normal cohort.
#' @example
#' m = 5; n = 10
#' dat = matrix(1:(m*n),m,n) # input data
#' out = 1-matrix(.99*round(.6*runif(length(dat))),m,n) # p-statistic of each
#' # data entry for being a potential outlier and ignored during imputation
#' imp = artImpute(dat = dat, ku = 2) # imputed data
artImpute = function(dat, ku = 6, marker.proc.list = NULL, miss.pstat = 4E-1,
    verbose = FALSE) {
    out.pstats = outScores(dat)
    m = nrow(dat); n = ncol(dat)
    dat.imp = data.frame(matrix(nrow = nrow(dat), ncol = ncol(dat)));
    rownames(dat.imp) = rownames(dat); colnames(dat.imp) = colnames(dat)
    if (!methods::is(dat, 'matrix')) {dat = as.matrix(dat)}
    if (is.null(marker.proc.list)){marker.proc.list = seq_len(nrow(dat))}
    k = min(ku, nrow(dat))
    dat.mis.ref = dat
    # remove all outlying expressions to not skew imputation
    dat.mis.ref[out.pstats < miss.pstat] = NA
    dat.dis = as.matrix(dist(dat.mis.ref, method = 'euclidean', diag = TRUE,
        upper = TRUE))
    dat.cor = as.matrix(cor(t(dat.mis.ref), use = 'pairwise.complete.obs'))
    for (i in marker.proc.list) {
        # choose k euclidean neighbors
        sorted = sort(as.numeric(dat.dis[i,]), decreasing = FALSE,
            index.return = TRUE, na.last = TRUE)
        nei.ix = sorted$ix[2:(k+1)]
        # discard uncorrelated neighbors
        nei.ix = nei.ix[dat.cor[i,nei.ix]>.25]
        dat.knn = dat.mis.ref[nei.ix,]
        # weighted average
        weights = (1/sorted$x[nei.ix])**2
        if (length(weights)==1) {
            dat.imp[i,] = dat.knn * weights
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
                print(paste('Background inference: ',
                round(100*i/length(marker.proc.list)),'% done.'))
            }
        }
    }
    return(as.matrix(dat.imp))
}
#' @title Analyze protruding (dysregulated) events
#' @description For each marker processed, draws a 2D scatter plot of maching
#'  values of observed vs imputed expressions.
#' @param dat, an object of log2-normalized gene (or protein) expressions,
#'  markers in rows and samples in columns.
#' @param dat.imp, the imputed data that likely represents the expressions of
#'  the markers in the (matched) normal samples.
#' @param marker.proc.list, the rownames of the data to be processed/imputed.
#' @return dat.dys, samples' distances to regression line (i.e., dysregulation)
#'  in the scatter plots.
#' @return plot.list, the scatter plots.
dysReg = function(dat, dat.imp, marker.proc.list, verbose = FALSE){
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
            print(paste('Dysregulation analysis:', round(100*i/m),'% done.'))
            }
        }
    }
    return(list(dat.dys, plot.list))
}
#' @title Analyze dysregulation significance
#' @description Rank markers by the significance of deviation of the observed
#'  expressions from the (matched) imputed expressions based on the
#'  Kolmogorov-Smirnov (KS) test.
#' @param dat, an object of log2-normalized gene (or protein) expressions,
#'  markers in rows and samples in columns.
#' @param dat.imp, the imputed data that likely represents the expressions of
#'  the markers in the (matched) normal samples.
#' @param marker.proc.list, the rownames of the data to be processed/imputed.
#' @param pval.insig, p-value threshold to determine spurious (null)
#'  dysregulation events.
#' @return dat.imp.test, p-values of the markers significance computed by the
#'  KS test.
#' @return dat.imp.test.sig, ranked p-values (KS test) of the significant
#'  markers that are lower than pval.insig.
#' @return markers.imp.sig, ranked significantly-dysregulated markers with
#'  p-values lower than pval.insig.
#' @return dat.imp.test.insig, ranked p-values (KS test) of the insignificant
#'  markers that are greater than pval.insig.
#' @return markers.imp.insig, ranked markers exhibiting spurious (null)
#'  dysregulation events with p-values greater than pval.insig.
statTest = function(dat, dat.imp, marker.proc.list, pval.insig = 2E-1) {
    m = nrow(dat); n = ncol(dat)
    # Detect weird markers O8: Kolmogorov-Smirnov test, compare empirical CDFs
        # of Observed vs Imputed
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
#' @param marker.proc.list, the given markers to be processed.
#' @param dat, an object of log2-normalized gene (or protein) expressions,
#'  markers in rows and samples in columns.
#' @param dat.imp, the imputed data that likely represents the expressions of
#'  the markers in the (matched) normal samples.
#' @param dat.imp.test, p-values of the markers significance computed by the KS
#'  test.
#' @param dat.dys, samples' distances to regression line (i.e., dysregulation)
#'  in the scatter plots.
#' @param dys.sig.thr.upp, the dysregulation score threshold to elucidate/mark
#'  significantly dysregulated outlier events.
#' @param dataset, the cohort name to be used in the output files.
#' @param num.omit.fit, number of outlying events to ignore when fitting a
#'  marker's observed expressions to the imputed ones.
#' @param draw.sc, if TRUE draws a scatter plot for every marker in
#'  marker.proc.list in a separate PDF file.
#' @param draw.vi, if TRUE draws a violin plot for every marker in
#'  marker.proc.list in a separate PDF file.
#' @return plot.list.marked, the scatter plots of the markers where the outlier
#'  dysregulation events are highlighted by red mark.
markOut = function(marker.proc.list, dat, dat.imp, dat.imp.test, dat.dys,
    dys.sig.thr.upp, dataset, num.omit.fit=NULL, draw.sc=TRUE, draw.vi=TRUE){
    if(is.null(num.omit.fit)) {num.omit.fit = round(.1*ncol(dat))}
    # To align the scatter plot margins in all markers, use:
    # minl = min(min(dat.imp[marker.proc.list,], na.rm=TRUE),
        # min(dat[marker.proc.list,], na.rm=TRUE))
    # maxl = max(max(dat.imp[marker.proc.list,], na.rm=TRUE),
        # max(dat[marker.proc.list,], na.rm=TRUE))
    plot.list.marked = list()
    for (marker in marker.proc.list) {
        marker.loc = which(rownames(dat)==marker)
        # significant outlying events for the given marker
        marker.out.exp.loc = which(dat.dys[marker,]>dys.sig.thr.upp)
        # | dat.dys[marker,]<dys.sig.thr.low)
        marker.out.samples = colnames(dat.dys)[marker.out.exp.loc]
        dat.dys[marker,marker.out.exp.loc]
        out = gqplot(y = dat[marker,], x = dat.imp[marker,], ci = .95,
                    ylab = 'Observed', xlab = 'Imputed',
                    highlight = marker.out.samples, omit.fit = num.omit.fit)
                    #, minl = minl, maxl = maxl)
        d=out[[1]]; plot.it=out[[2]]
        plot.list.marked[[marker.loc]] = plot.it;
        # draw scatter plot
        if (draw.sc) {
            pdf(paste(dataset,'.',marker,'.sc','.dysreg.pdf',sep=''),w=4,h=4)
            print(plot.list.marked[[marker.loc]])
            dev.off()
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
            pdf(paste(dataset,'.',marker,'.vi','.dysreg.pdf',sep=''),w=4,h=4)
            print(pl)
            dev.off()
        }
    }
    return(plot.list.marked)
}
#' @title Rank markers by the percentage of outlying events
#' @description Ranks markers in the order of decreasing percentage of outlying
#'  events.
#' @param dat.dys, samples' distances to regression line (i.e., dysregulation)
#'  in the scatter plots.
#' @param marker.proc.list, the given markers to be processed.
#' @param dys.sig.thr.upp, the dysregulation score threshold to elucidate/mark
#'  significantly dysregulated outlier events.
#' @return ranked markers are returned in marker.out.exp.per.sor.names with
#'  corresponding percentages in marker.out.exp.per.
rankPerOut = function(dat.dys, marker.proc.list, dys.sig.thr.upp){
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
        marker.out.exp.per.sor$ix[seq_len(length(marker.proc.list))]]
    return(list(marker.out.exp.per.sor.names, marker.out.exp.per))
}
#' @title Hierachical cluster analysis
#' @description Displays the hierarchically clustered data by the "pheatmap"
#'  package.
#' The number of cluster along the markers/samples are estimated by a
#'  multivariate algorithm "cafr" package (https://github.com/weiyi-bitw/cafr).
#' Optionally, they can be set by the user then cluster structures are
#'  estimated by pair-wise analysis of the markers/samples.
#' @param num_clusters_row the number of clusters along the markers
#' @param num_clusters_col the number of clusters along the samples
#' @return the hierarchical tree, the cluster identities of the markers and
#'  samples are returned in tree, cluster_IDs_row, and cluster_IDs_col,
#'  respectively.
clusterData = function(data, annotation_col = NULL, annotation_row = NULL,
    annotation_colors = NULL, main = NA,
    stringency_col = 6, stringency_row = 4,
    clustering_distance_cols = 'euclidean',
    clustering_distance_rows = 'euclidean',
    display_numbers = FALSE, num_clusters_col = NULL, num_clusters_row = NULL,
    cluster_cols = TRUE, cluster_rows = TRUE,
    annotate_new_clusters_col = FALSE, zero_white = FALSE){
    if (!is.null(num_clusters_col) | !cluster_cols) {att_col = NULL}
    if (!is.null(num_clusters_row) | !cluster_rows) {att_row = NULL}
    if (is.null(num_clusters_col) & cluster_cols)
        {att_col = cafr::attractorScanning(as.matrix(t(data)), maxIter = 1E2,
        epsilon = 1E-14, a = stringency_col); num_clusters_col=dim(att_col)[1]}
    if (is.null(num_clusters_row) & cluster_rows)
        {att_row = cafr::attractorScanning(as.matrix(data), maxIter = 1E2,
        epsilon = 1E-14, a = stringency_row); num_clusters_row=dim(att_row)[1]}
    if (is.null(num_clusters_col)) {num_clusters_col = 1}
    if (is.null(num_clusters_row)) {num_clusters_row = 1}
    if (zero_white) {
        paletteLength = 100
        my.color = grDevices::colorRampPalette(
            c("#006699", "white", "red"))(paletteLength)
        my.breaks = c(seq(min(data, na.rm = TRUE), 0,
            length.out = ceiling(paletteLength/2) + 1), seq(max(data,
            na.rm = TRUE)/paletteLength, max(data, na.rm = TRUE),
            length.out = floor(paletteLength/2)))
    } else {
        my.color = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(
            n = 7, name = "RdYlBu")))(100)
        my.breaks = NA
    }

    tree = pheatmap::pheatmap(data, annotation_col = annotation_col,
        annotation_colors = annotation_colors,
        display_numbers = display_numbers,
        clustering_distance_cols = clustering_distance_cols,
        clustering_distance_rows = clustering_distance_rows,
        cluster_cols = cluster_cols, cluster_rows = cluster_rows,
        cutree_cols = num_clusters_col, show_colnames = TRUE,
        cutree_rows = num_clusters_row, show_rownames = TRUE,
        main = main, color = my.color, breaks = my.breaks)
    if (cluster_cols) {cluster_IDs_col = cutree(tree$tree_col,
        k = num_clusters_col)} else {cluster_IDs_col = NULL}
    if (cluster_rows) {cluster_IDs_row = cutree(tree$tree_row,
        k = num_clusters_row)} else {cluster_IDs_row = NULL}
    if (annotate_new_clusters_col) {
        annotation_col_new = cbind(annotation_col,data.frame(cluster_IDs_col)[
            rownames(annotation_col),])
        colnames(annotation_col_new)[length(annotation_col_new)] =
            'Putative.subtype'
        tree = pheatmap::pheatmap(data, annotation_col = annotation_col_new,
            annotation_colors = annotation_colors,
            display_numbers = display_numbers,
            clustering_distance_cols = clustering_distance_cols,
            clustering_distance_rows = clustering_distance_rows,
            cluster_cols = cluster_cols, cluster_rows = cluster_rows,
            cutree_cols = num_clusters_col, show_colnames = TRUE,
            cutree_rows = num_clusters_row, show_rownames = TRUE,
            main = main, color = my.color, breaks = my.breaks)
    }
    return(list(tree, cluster_IDs_row, cluster_IDs_col)) #, att_row, att_col
}
#' @title Analyze protein / phosphosite expressions
#' @description Find outlying markers and events across cancer types.
#' @param file.path, character string, path to a data list where each element
#'  contains proteomics data for a different cohort (markers in the rows,
#'   samples in the columns)
#' @param panel, character string in c('global', 'pancan', 'drugge', 'immune',
#'  'kinase'), a panel of markers to analyze accross cohorts
#' @param cohort.names, character array
#' @param tol.nas, a constant in [0,100], tolerance for the percentage of NAs
#'  in a marker, e.g., tol.nas = 20 will filter out markers contatining 20% or
#'   more NAs
#' @param ku, an integer in [1,num.markers], upper bound on the number of
#'  nearest neighbors
#' @param miss.pstat, a constant in [0,1], statistic to estimate potential
#'  outliers, 0.4 ~= q75+1.5*IQR
#' @param demo.panels, logical, to draw demographics of the panel in each
#'  cancer cohort
#' @param save.data, logical, to save intermediate data (background inference
#'  and dysregulation measures)
#' @param draw.sc.plots, logical, to draw each marker's qqplot of observed vs
#'  inferred (imputed) expressions
#' @param draw.vi.plots, logical, to draw each marker's violin plot of observed
#'  vs imputed expressions
#' @param draw.ou.plots, logical, to draw each marker's outlier prevalence
#' (by the percentage of outlying samples) across the cohorts
#' @param verbose, logical, to show progress of the algorithm
#' @return pan.dat.dys, dysregulation scores of every marker per each sample;
#'  pan.dat.imp.test, the result of KS tests that evaluates the statistical
#'  significance of each marker's outlier samples;
#'  pan.marker.out.exp.per, a data list contatining, for each cohort, the
#'  percentage of outlier samples for every marker
oppti = function(data, panel = 'global', cohort.names = NULL, tol.nas = 20,
    ku = 6, miss.pstat = .4, demo.panels = FALSE,
    save.data = FALSE, draw.sc.plots = FALSE,
    draw.vi.plots = FALSE, draw.ou.plots = FALSE,
    verbose = FALSE) {
    # Load data
    if (is.character(data)){tryCatch({data = readRDS(data)},
        warning=function(w){print(w)}, finally={})} else if (!is.list(data)){
        warning('Unrecognized format for the "data" object; must be a "list",
        or a "character" string referring to the file path of such object.');
        return()}
    if (!methods::is(data, 'list')) {data = list(data)}
    if (is.null(cohort.names)) {cohort.names = unlist(lapply(data, function(x){
        x=strsplit(colnames(x)[1],'\\.')[[1]][1]}))} #cohort names
    pan.num = length(cohort.names) #number of cohorts
    tmp.lis = as.list(rep(NA,pan.num)) #template for a cohort-size data object
    # Filter out markers with >= tol.nas percent missing values
    pan.dat = lapply(data, function(x){x=dropMarkers(x,
        percent_NA = tol.nas/100, low_mean_and_std = NULL,
        q_low_var = NULL)[[1]]})
    # Choose samples (tumor or normal) for analyses
    pan.dat = lapply(pan.dat, function(x){if (any(regexpr('Tumor',
        colnames(x))>0)){x=x[,regexpr('Tumor', colnames(x)) > 0]} else {x=x}})
    pan.mar = rownames(pan.dat[[1]]); if (pan.num>1) {for (i in 2:pan.num) {
        pan.mar = intersect(pan.mar, rownames(pan.dat[[i]]))}}
        #determine markers quantified pan-cancer
    # Load panel markers
    switch(panel
        , 'global' = (panel.markers = unique(unlist(lapply(data,function(x)
            {x=rownames(x)}))))
        , 'pancan' = (panel.markers = pan.mar)
        , 'drugge' = (panel.markers = readRDS(file = 'data/drug.genes.RDS'))
        , 'kinase' = (panel.markers = readRDS(file = 'data/kinase.genes.RDS'))
        , 'immune' = (panel.markers =
            readRDS(file = 'data/immune.checkpoint.genes.RDS'))
    )
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
            {pdf(paste('',cohort.names[i],'.',panel,'.pdf',sep=''), w=6,h=4);
            graphics::hist(pan.nas[[i]][rownames(data[[i]]) %in% panel.markers]
            , xlab = '% NAs', ylab = 'Number of markers'
            , main = paste(cohort.names[i],' ', panel, ' > ', pan.det[[i]]
            , ' detected > ', pan.com[[i]], ' complete (<', tol.nas,'% NAs)'
            , sep='')); dev.off()}
        for (i in seq_len(pan.num)) {plotDen(dat = data[[i]],
            name = cohort.names[i], per.plot = NULL, legend = FALSE)}
        plotDen(dat = cbindNA(data), group = cohort.names, name = 'data')
    }
    # Generate inferred expressions by weigted average of nearest neighbors
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
    if (demo.panels) {
        colors = c('red','orange','yellow','green','blue','purple');
        limx=1; limy=1; pdf('pan.null.dys.ecdf.pdf', w=6,h=6);
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
    if (draw.sc.plots | draw.vi.plots) {for (i in seq_len(pan.num))
        {markOut(marker.proc.list = pan.proc.markers[[i]], pan.dat[[i]],
        pan.dat.imp[[i]], pan.dat.imp.test[[i]], pan.dat.dys[[i]],
        pan.dys.sig.thr.upp[[i]], cohort.names[i], draw.sc = draw.sc.plots,
        draw.vi = draw.vi.plots)}}
    # Rank markers by the percentage of outlying events
    pan.marker.out.exp.per = tmp.lis; for (i in seq_len(pan.num))
        {pan.marker.out.exp.per[[i]] = rankPerOut(pan.dat.dys[[i]],
        pan.proc.markers[[i]], pan.dys.sig.thr.upp[[i]])[[2]]}
    if (draw.ou.plots) {
        # Draw markers' percentages of outlying events for each cancer
        pan.mar.out.exp.per =
            data.frame(matrix(NA,length(panel.markers),pan.num))
        colnames(pan.mar.out.exp.per) = cohort.names
        rownames(pan.mar.out.exp.per) = panel.markers
        for (i in seq_len(pan.num))
            {pan.mar.out.exp.per[rownames(pan.marker.out.exp.per[[i]]),i] =
                pan.marker.out.exp.per[[i]]}
        # Pan-cancer highly-outlying markers
        pan.mar.out.exp.per.ratio.sor = sort(rowSums(pan.mar.out.exp.per,
            na.rm = TRUE)/ncol(pan.mar.out.exp.per), decreasing = TRUE,
            index.return = TRUE, na.last = TRUE)
        tmp = as.data.frame(pan.mar.out.exp.per[
            pan.mar.out.exp.per.ratio.sor$ix[
                pan.mar.out.exp.per.ratio.sor$x>0],])
        rownames(tmp) = rownames(pan.mar.out.exp.per)[
            pan.mar.out.exp.per.ratio.sor$ix[
                pan.mar.out.exp.per.ratio.sor$x>0]]
        colnames(tmp) = colnames(pan.mar.out.exp.per)
        pan.mar.ranked.out.exp.per.tree = clusterData(tmp,
            cluster_cols = FALSE, cluster_rows = FALSE,
            display_numbers = FALSE, main = 'In-cohort outlier perce.')[[1]]
        pdf(paste('pan.cancer.',panel,'.markers.outlier.scores.pdf', sep = ''),
            w = max(2,ceiling((ncol(tmp))**(.7)-1)),
            h = max(2,ceiling((nrow(tmp))**(.7)-2)))
        print(pan.mar.ranked.out.exp.per.tree); dev.off()
        # Pan-cancer top-20 highly-outlying markers
        tmp = t(pan.mar.out.exp.per[pan.mar.out.exp.per.ratio.sor$ix[
            pan.mar.out.exp.per.ratio.sor$x>0][seq_len(min(20,length(
                pan.mar.out.exp.per.ratio.sor$ix[
                    pan.mar.out.exp.per.ratio.sor$x>0])))],])
        colnames(tmp) = rownames(pan.mar.out.exp.per)[
            pan.mar.out.exp.per.ratio.sor$ix[
                pan.mar.out.exp.per.ratio.sor$x>0][seq_len(min(20,length(
                    pan.mar.out.exp.per.ratio.sor$ix[
                        pan.mar.out.exp.per.ratio.sor$x>0])))]]
        rownames(tmp) = colnames(pan.mar.out.exp.per)
        pan.mar.ranked20.t.out.exp.per.tree = clusterData(tmp,
            cluster_cols = FALSE, cluster_rows = FALSE, display_numbers = TRUE,
            main = 'In-cohort outlier perce.')[[1]]
        pdf(paste('pan.cancer.',panel,'.top20.highly.outlying.markers.pdf',
            sep = ''), w = max(2,ceiling((ncol(tmp))**(.7)-1)),
            h = max(2,ceiling((nrow(tmp))**(.7)-2)))
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
                display_numbers = TRUE, main = 'In-cohort outlier perce.')[[1]]
            pdf(paste('pan.cancer.',panel,
                '.top20.variably.outlying.markers.pdf', sep = ''),
                w = max(2,ceiling((ncol(tmp))**(.7)-1)),
                h = max(2,ceiling((nrow(tmp))**(.7)-2)))
            print(pan.mar.ranked20.t.sd.out.exp.per.tree); dev.off()
        }
    }
    if (pan.num>1){
        return(list(pan.dat.dys, pan.dat.imp.test, pan.marker.out.exp.per))
    } else {
        return(list(pan.dat.dys[[1]], pan.dat.imp.test[[1]],
                    pan.marker.out.exp.per[[1]]))
    }
}
