# TODO for jamma

## ggjammaplot()

* make empty plot panels completely blank, relevant when using `blankPlotPos`
* set `xlab`, `ylab` using summary, difference labels
* add subtitleBox to bottomleft corner
* investigate using `element_text()` instead of `ggtext::element_textbox()`
* allow `geom_text_repel()` to label highlighted points
* consider selectable x- and y-ranges, to highlight a box and points within it

* Needed a workaround to build pkgdown site, see: https://github.com/r-lib/pkgdown/issues/1157

## New functions

* `jammaplotDispEsts()` wrapper for `DESeq2::plotDispEsts()`, although
this function needs the newer colorized `plotSmoothScatterG2()`.


## compatibility with SummarizedExperiment


## Bug fixes

* `jammaplot()` argument `ablineV` is not functioning properly.
* `jammaplot()` highlight legend is not using `highlightPch`.

## Vignettes

1. Guides for MA-plots.

* Basic guide to MA-plots for gene expression data.

  * check data normality, apply appropriate transform, normalize data
  * check within-group variability
  * highlight points
  * common patterns and what they mean:
  
    * shifted up/down
    * skewed up/down
    * low signal-to-noise, no signal
    * the 45-degree lines
    * batch effects
    * technical versus biological controls
    
  * non-parametric (rank-based) MA-plot

* When to do use `centerGroups`, `controlSamples`.
* How to interpret global-, group-, and technical replicate-centered data.
* How to detect sample outliers using a MAD factor threshold.
* How to guide data normalization by MA-plot review.


## useful functions

* `log2fold_to_fold()` and `fold_to_log2fold()` - to interconvert:

   * normal space fold changes, which could be represented as positive and
   negative values (e.g. 2-fold and -2-fold) or as ratios
   (e.g. 2-fold and 0.5-fold); and
   * log2 fold change values
   * Consider `fold_sign()` function that takes either as input and
   returns either `1` or `-1`, and by default never returns `0` since
   1-fold change cannot easily be multiplied by its `fold_sign()` without
   causing it to become zero. More thought to be had as to the workflow.


## Interface with other R packages

### DESeq2

* `plotMA()` is used for the `DESeq2::DESeqResults` object,
highlighting points that meet `alpha < 0.1` by default.
* Add new method to run `jammaplot()` on the count data,
optionally transformed using their recommended approach,
or use `useRank=TRUE`.


## Bug fixes

* Version 0.0.10.900 fixed a bug where rowGroupMeans() was used to
center values, but used default `na.rm=FALSE`, which caused groups
with missing values not to display a centered value for other non-NA
samples. The new version should provide two options:

    1. Hide or display values where groups with n>1 members have
    only one non-NA value. The reasoning is that the MA-plot is
    intended to show "difference from mean" and has far less utility
    when a fraction of rows has only one non-NA value in the group, 
    thus emphasizing the display of "zero difference" in those
    samples. Hiding these values therefore helps display the
    variability where data is available to judge the difference
    from mean.
    2. Hide groups with any NA values.



## Refactoring ideas

* Consider using `jamba::rowGroupMeans()` within `centerGeneData()`
which would allow optional outlier detection, which could substantially
improve the quality of MA-plot panels.
* `jammagg()` or something similarly named, to produce a ggplot2 object.
* One goal for ggplot2 output is to produce
interactive visualizations with plotly; however, if that is not
sufficient for the desired properties, e.g. useful hover text over
the smooth scatter density, then a custom interactive plot function
may be needed, for example direct calls to plotly instead of using
the ggplotly wrapper function.


### Update to handle SummarizedExperiment objects

* centerGroups should be able to take colnames from colData(x)
to define groupings.
* `jammaplot.SE()` could be specific to SummarizedExperiment objects.
It would require the names(assays(SE)) to define the data matrix to use.
It could even recognize sample groups from `colData(SE)`, or
from internal design matrix used for statistical testing.


### Highlighted points

* Consider allowing highlighted points per panel, for example,
showing gene hits only in the panels relevant to that comparison.
However, there is a workaround, using `whichSamples` which will
create the MA-plot data, but only display the samples of interest.


### Color legend for highlighted points

* Add optional color legend at bottom of the figure labeling the
colors highlighted on the plot.

### Enhance the returned data

* Currently returns a list of 2-column matrices
* Return data.frame with properly labeled colnames
* Add columns:
    * indicating the centerGroups value if defined
    * whether the sample is a control sample within its centerGroup
    * highlight label, highlight point color
    * title box color

User should be able to take the returned data.frame list, and make
a tall data.frame sufficient for ggplot2, or sufficient to answer
the question "What is that outlier datapoint?"


### Interactive plots using plotly

Future idea:

* Convert plotSmoothScatter to generate a heatmap-ready density
plot for each panel.
* For each cell in the panel, optionally update the label so it
lists the point label(s) represented in that cell; alternatively
simply label with the number of points represented.

### Accessory function to select points

* A common workflow request stems from the question "What are those
points?" A wrapper function was written to allow clicking twice inside
an active plot device to define a rectangle, then returning the
rownames of the corresponding points.
* The steps above could be used to supply coordinates from something
like R-shiny, to request the rownames of the points specified.
