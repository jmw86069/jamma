## TODO for jamma

### Update to handle SummarizedExperiment objects

* centerGroups should be able to take colnames from colData(x)
to define groupings.
* 


### Outlier detection enhancements

* consider calculating per-sample MAD values within centerGroups,
so the variability of another centerGroup does not define the MAD
used across all centerGroups.

### Highlighted points

* Consider allowing highlighted points per panel, for example,
showing gene hits only in the panels relevant to that comparison.
However, there is a workaround, using `whichSamples` which will
create the MA-plot data, but only display the samples of interest.

### Base plotting enhancements

* Add argument for title and subtitle, then use outer margin to
display one title acros the top of the plot.
* change `colrampOutlier` to use more bold yellow color, for example
c("palegoldenrod", "lightblue", "blue", "navy", "orange", "orangered2").
However, it may be best to change the first color in `colramp` to
`"palegoldenrod"`, maintaining the rest of the color ramp as-is.
That mechanism would ensure that custom `colramp` would be consistent
with `colrampOutlier`.

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


