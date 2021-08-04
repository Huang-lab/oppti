# Analyze proteomics data of a single cohort

You can easily analyze outlying (dysregulated) markers for each sample
in a cohort. Lets generate a toy proteomics data for a cohort of 30
disease samples, each quantifying 100 proteins.

``` r
set.seed(1)
cohort1.proteomes = as.data.frame(matrix(abs(rnorm(100*30)), 100, 30)) 
rownames(cohort1.proteomes) = paste('marker', 1:100, sep = '')
colnames(cohort1.proteomes) = paste('cohort1.sample', 1:30, sep = '')
```

Outlier analysis is run by the `oppti` function:

``` r
library('oppti')
result = oppti(cohort1.proteomes)
```

The outlier scores of each marker in each sample are then returned in
the first element of the result:

``` r
cohort1.outlier.scores = result[[1]] 
```

|          | cohort1.sample1 | cohort1.sample2 | cohort1.sample3 | cohort1.sample4 |
| :------- | --------------: | --------------: | --------------: | --------------: |
| marker1  |            0.13 |          \-0.28 |          \-0.29 |            0.13 |
| marker2  |          \-0.12 |          \-0.28 |            0.29 |          \-0.10 |
| marker3  |          \-0.12 |            0.16 |            0.48 |            0.53 |
| marker4  |            0.73 |          \-0.60 |          \-0.73 |          \-0.15 |
| marker5  |          \-0.16 |            0.10 |            1.62 |            0.77 |
| marker6  |            0.07 |            0.38 |            0.75 |            0.55 |
| marker7  |          \-0.01 |          \-0.29 |          \-0.02 |          \-0.44 |
| marker8  |            0.17 |            0.29 |            0.00 |          \-0.01 |
| marker9  |          \-0.18 |          \-0.46 |          \-0.77 |          \-0.01 |
| marker10 |          \-0.29 |            0.69 |          \-0.36 |          \-0.55 |

Example matrix of the outlier scores, displayed for the first 10
proteins (rows) and the first 4 samples (columns)

In this toy example, marker5 has a (somewhat) elevated outlier score in
sample3, suggesting a protruding expression in the disease state of
sample3 relative to a normal state (i.e., the consensus co-expression
network inferred for marker5). In contrast, a negative sign in the
outlier score indicates a negative dysregulation event, i.e., relatively
“lower” protein expression is expected in the observed disease state
compared to the normal state. The landscape of these aberrant
expressions analyzed for a cohort of individuals may serve for the
discovery of personalized actionable targets.

The outlier scores correspond to deviations of the observed expressions
from the estimated normal states. The estimated normals are given in the
second element of the result:

``` r
cohort1.normal.states = result[[2]] 
```

|          | cohort1.sample1 | cohort1.sample2 | cohort1.sample3 | cohort1.sample4 |
| :------- | --------------: | --------------: | --------------: | --------------: |
| marker1  |            0.16 |            0.89 |            0.59 |            0.56 |
| marker2  |            0.50 |            0.61 |            0.85 |            0.94 |
| marker3  |            0.85 |            0.52 |            0.77 |            1.11 |
| marker4  |            0.70 |            0.83 |            1.43 |            0.31 |
| marker5  |            0.04 |            0.24 |            0.39 |            1.31 |
| marker6  |            0.31 |            1.23 |            1.69 |            0.51 |
| marker7  |            0.00 |            2.08 |            0.74 |            0.15 |
| marker8  |            0.52 |            0.58 |            0.52 |            0.59 |
| marker9  |            0.36 |            0.62 |            0.54 |            0.92 |
| marker10 |            0.25 |            1.18 |            1.03 |            1.06 |

Example matrix of the normal states

You can evaluate markers by the odds of obtaining these deviations
purely by chance. A Kolmogorov-Smirnov test is performed for each marker
between its observed and estimated states, and the p-values are reported
in the third element of the result:

``` r
cohort1.markers.tests = result[[3]] 
```

|          |      x |
| :------- | -----: |
| marker1  | 0.8080 |
| marker2  | 0.0346 |
| marker3  | 0.3929 |
| marker4  | 0.3929 |
| marker5  | 0.2391 |
| marker6  | 0.1350 |
| marker7  | 0.9988 |
| marker8  | 0.0709 |
| marker9  | 0.3929 |
| marker10 | 0.9578 |

Statistical significance of outlying markers

# Analyze proteomics data of multiple cohorts

For pan-cancer analyses, the normalized proteomics data from different
cohorts can be supplied to `oppti` in a list object. Lets generate
another toy proteomics data for a separate cohort of 20 disease samples,
each quantifying 80 proteins (say, 50 of which are overlapping with
those quantified in the first cohort).

``` r
cohort2.proteomes = as.data.frame(matrix(abs(rnorm(80*20)), 80, 20)) 
rownames(cohort2.proteomes) = paste('marker', 51:130, sep = '')
colnames(cohort2.proteomes) = paste('cohort2.sample', 31:50, sep = '')
```

To run `oppti` for both cohorts, the data are simply fed in a single
list object:

``` r
result = oppti(list(cohort1.proteomes,cohort2.proteomes))
```

Again, the outlier scores of each marker in each sample are returned in
the first element of the result.

``` r
outlier.scores = result[[1]]
```

However, this object is a list of 2 elements per se, corresponding to
two cohorts. To obtain the outlier scores of the first cohort:

``` r
cohort1.outlier.scores = outlier.scores[[1]]
```

|          | cohort1.sample1 | cohort1.sample2 | cohort1.sample3 | cohort1.sample4 |
| :------- | --------------: | --------------: | --------------: | --------------: |
| marker1  |            0.13 |          \-0.28 |          \-0.29 |            0.13 |
| marker2  |          \-0.12 |          \-0.28 |            0.29 |          \-0.10 |
| marker3  |          \-0.12 |            0.16 |            0.48 |            0.53 |
| marker4  |            0.73 |          \-0.60 |          \-0.73 |          \-0.15 |
| marker5  |          \-0.16 |            0.10 |            1.62 |            0.77 |
| marker6  |            0.07 |            0.38 |            0.75 |            0.55 |
| marker7  |          \-0.01 |          \-0.29 |          \-0.02 |          \-0.44 |
| marker8  |            0.17 |            0.29 |            0.00 |          \-0.01 |
| marker9  |          \-0.18 |          \-0.46 |          \-0.77 |          \-0.01 |
| marker10 |          \-0.29 |            0.69 |          \-0.36 |          \-0.55 |

Example outlier scores in cohort1

Similarly, for the second cohort the outlier scores are obtained by:

``` r
cohort2.outlier.scores = outlier.scores[[2]]
```

|          | cohort2.sample31 | cohort2.sample32 | cohort2.sample33 | cohort2.sample34 |
| :------- | ---------------: | ---------------: | ---------------: | ---------------: |
| marker51 |             0.20 |           \-0.41 |             0.15 |           \-0.20 |
| marker52 |           \-0.26 |           \-0.36 |           \-0.40 |             0.11 |
| marker53 |             0.20 |             0.24 |           \-0.01 |           \-0.17 |
| marker54 |             0.10 |           \-0.18 |             0.31 |           \-0.16 |
| marker55 |             0.24 |           \-0.15 |             0.58 |             0.61 |
| marker56 |           \-0.31 |             0.05 |           \-0.17 |             0.23 |
| marker57 |             0.38 |           \-0.10 |             0.26 |           \-0.38 |
| marker58 |           \-0.36 |           \-0.03 |             0.23 |             0.55 |
| marker59 |             0.13 |             0.15 |           \-0.11 |           \-0.02 |
| marker60 |             0.02 |           \-0.35 |             0.26 |           \-0.28 |

Example outlier scores in cohort2

You can evaluate the markers in terms of outlying events they exhibit
across the cohort by using the `draw.sc.plots` flag. The outlier samples
will be marked on a scatter plot displaying disease (observed) vs normal
(estimated) expressions. Note that you can always set `panel.markers`
parameter to restrict your analysis to a specific set of markers.

``` r
result = oppti(list(cohort1.proteomes,cohort2.proteomes), draw.sc.plots = TRUE,
    panel.markers = rownames(cohort1.proteomes)[46:55])
```

To display the summary results of the markers’ outlying events across
cohorts you can use `draw.ou.plots`:

``` r
result = oppti(list(cohort1.proteomes,cohort2.proteomes), draw.ou.plots = TRUE,
    panel.markers = rownames(cohort1.proteomes)[46:55])
```

To narrow down the summary results to a number of markers you can use
`draw.ou.markers`:

``` r
result = oppti(list(cohort1.proteomes,cohort2.proteomes), 
    draw.ou.markers = c('marker50', 'marker55'), 
    panel.markers = rownames(cohort1.proteomes)[46:55])
```
