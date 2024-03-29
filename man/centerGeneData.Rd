% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jamma-center.R
\name{centerGeneData}
\alias{centerGeneData}
\title{Center gene data}
\usage{
centerGeneData(
  x,
  centerGroups = NULL,
  na.rm = TRUE,
  controlSamples = NULL,
  useMedian = TRUE,
  rmOutliers = FALSE,
  madFactor = 5,
  controlFloor = NA,
  naControlAction = c("na", "row", "floor", "min"),
  naControlFloor = 0,
  rowStatsFunc = NULL,
  returnGroupedValues = FALSE,
  returnGroups = FALSE,
  mean = NULL,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{x}{\code{numeric} matrix of input data. See assumptions,
that data is assumed to be log2-transformed, or otherwise
appropriately transformed.}

\item{centerGroups}{\code{character} vector of group names, or
\code{NULL} if there are no groups.}

\item{na.rm}{\code{logical} indicating whether NA values should be
ignored for summary statistics. This argument is passed
to the corresponding row stats function. Frankly, this
value should be \code{na.rm=TRUE} for all stat functions by default,
for example \code{mean(..., na.rm=TRUE)} should be default.}

\item{controlSamples}{\code{character} vector of values in \code{colnames(x)}
which defines the columns to use when calculating group
summary values.}

\item{useMedian}{\code{logical} indicating whether to use group median
values when calculating summary statistics \code{TRUE}, or
group means \code{FALSE}. In either case, when \code{rowStatsFunc}
is provided, it is used instead.}

\item{rmOutliers}{\code{logical} indicating whether to perform outlier
detection and removal prior to row group stats. This
argument is passed to \code{jamba::rowGroupMeans()}. Note that
outliers are only removed during the row group summary step,
and not in the centered data.}

\item{madFactor}{\code{numeric} value passed to \code{jamba::rowGroupMeans()},
indicating the MAD factor threshold to use when \code{rmOutliers=TRUE}.
The MAD of each row group is computed, the overall group median
MAD is used to define 1x MAD factor, and any MAD more than
\code{madFactor} times the group median MAD is considered an outlier
and is removed. The remaining data is used to compute row
group values.}

\item{controlFloor}{\code{numeric} value used as a minimum for any control
summary value during centering.
Use \code{NA} to skip this behavior.
When defined, all control group summary values are calculated,
then any values below \code{controlFloor} are set to the \code{controlFloor}
for the purpose of data centering.
By default \code{controlFloor=NA} which imposes no such floor value.
However, \code{controlFloor=0} would be appropriate when zero is defined
as effective noise floor after something like background subtraction
during the upstream processing or upstream normalization.
Using a value above zero would be appropriate when the effective
noise floor of a platform is above zero, so that values are not
centered relative to noise. For example, if the effective noise
floor is 5, then centering should not "amplify" differences from
any value less than 5, since in this scenario a value of 5 or less
is effectively the same as a value of 5. It has the effect of returning
fold changes relative to the effective platform minimum detectable
signal.}

\item{naControlAction}{\code{character} string indicating how to handle the specific
scenario when the control group summary value is \code{NA} for a particular
centering operation.
\itemize{
\item \code{"na"}: default is to return \code{NA} since 15 - NA = NA.
\item \code{"row"}: use the summary value across all relevant samples,
so the centering is against all non-NA values within the center group.
\item \code{"floor"}: use the numeric value defined by \code{naControlFloor},
to indicate a practical noise floor for the centering operation.
When \code{naControlFloor=0} (default) this option effectively keeps
non-NA values without centering these values.
\item \code{"min"}: use the minimum control value as the floor, which effectively
defines the floor by the lowest observed summary value across all
rows. It assumes rows are generally on the same range of detection,
even if not all rows have the same observed range. For example,
microarray probes have reasonably similar theoretical range of
detection, even if some probes to highly-expressed genes are
commonly observed with higher signal. The lowest observed signal
effectively sets the minimum detected value.
}}

\item{rowStatsFunc}{\code{optional} function used to calculate row group
summary values. This function should take a numeric matrix as
input, and return a one-column numeric matrix as output, or
a numeric vector with length \code{nrow(x)}. The function should
also accept \code{na.rm} as an argument.}

\item{returnGroupedValues}{\code{logical} indicating whether to include
the numeric matrix of row group values used during centering,
returned in the attributes with name \code{"x_group"}.}

\item{returnGroups}{\code{logical} indicating whether to return the
centering summary data.frame in attributes with name "center_df".}

\item{verbose}{\code{logical} indicating whether to print verbose output.}

\item{...}{additional arguments are passed to \code{jamba::rowGroupMeans()}.}
}
\description{
Performs per-row centering on a numeric matrix
}
\details{
This function centers data by subtracting the median or
mean for each row.

Columns can be grouped using argument \code{centerGroups}.
Each group group of columns defined by \code{centerGroups}
is centered independently.

Data can be centered relative to specific control columns
using argument \code{controlSamples}.
When \code{controlSamples} is not supplied, the default behavior
is to use all columns. This process is consistent with
typical MA-plots.

It may be preferred to define \code{controlSamples} in cases where
there are known reference samples, against which other samples
should be compared.

The \code{controlSamples} logic is applied independently to each
group defined in \code{centerGroups}.

You can confirm the \code{centerGroups} and \code{controlSamples} are
correct in the result data, by accessing the attribute
\code{"center_df"}, see examples below.

Note: This function assumes input data is suitable for
centering by subtraction.
This data requirement is true for:
\itemize{
\item most log-transformed gene expression data
\item quantitative PCR (QPCR) cycle threshold (CT) values
\item other numeric data that has been suitably transformed
to meet reasonable parametric assumption of normality,
\item rank-transformed data which results in difference in rank
\item generally speaking, any data where the difference between 5 and 7 (2)
is reasonably similar to the difference between 15 and 17 (2).
\item it may be feasible to perform background subtraction on straight
count data, for example sequence coverage at a particular location
in a genome.
}

The data requirement is not true for:
\itemize{
\item most gene expression data in normal space
(hint: if any value is above 100, it is generally not log-transformed)
\item numeric data that is strongly skewed
\item generally speaking, any data where the difference between 5 and 7
is not reasonably similar to the difference between 15 and 17. If
the percent difference is more likely to be the interesting measure,
data may be log-transformed for analysis.
}

For special cases, \code{rowStatsFunc} can be supplied to perform
specific group summary calculations per row.
\subsection{Control groups with NA values (since version 0.0.28.900)}{

When \code{controlSamples} is supplied, and contains all \code{NA} values
for a given row of data, within relevant \code{centerGroups} subsets,
the default behavior is defined by \code{naControlAction="NA"} below:
\enumerate{
\item \code{naControlAction="na"}: values are centered versus \code{NA} which
results in all values \code{NA} (current behavior, default).
\item \code{naControlAction="row"}: values are centered versus the row,
using all samples in the same center group. This action effectively
"centers to what we have".
\item \code{naControlAction="floor"}: values are centered versus a \code{numeric}
floor defined by argument \code{naControlFloor}. When \code{naControlFloor=0}
then values are effectively not centered. However, \code{naControlFloor=10}
could for example be used to center values versus a practical noise
floor, if the range of detection for a particular experiment starts
at 10 as a low value.
\item \code{naControlAction="min"}: values are centered versus the minimum
observed summary value in the data, which effectively uses the data
to define a value for \code{naControlFloor}.
}

The motivation to center versus something other than \code{controlSamples}
when all measurements for \code{controlSamples} are \code{NA} is to have
a \code{numeric} value to indicate that a measurement was detected in
non-control columns. This situation occurs in technologies when
control samples have very low signal, and in some cases report
\code{NA} when no measurement is detected within the instrument range
of detection.
}
}
\examples{
x <- matrix(1:100, ncol=10);
colnames(x) <- letters[1:10];
# basic centering
centerGeneData(x);

# grouped centering
centerGeneData(x,
   centerGroups=rep(c("A","B"), c(5,5)));

# centering versus specific control columns
centerGeneData(x,
   controlSamples=letters[c(1:3)]);

# grouped centering versus specific control columns
centerGeneData(x,
   centerGroups=rep(c("A","B"), c(5,5)),
   controlSamples=letters[c(1:3, 6:8)]);

# confirm the centerGroups and controlSamples
x_ctr <- centerGeneData(x,
   centerGroups=rep(c("A","B"), c(5,5)),
   controlSamples=letters[c(1:3, 6:8)],
   returnGroups=TRUE);
attr(x_ctr, "center_df");

}
\seealso{
Other jam matrix functions: 
\code{\link{jammacalc}()},
\code{\link{jammanorm}()},
\code{\link{matrix_to_column_rank}()}
}
\concept{jam matrix functions}
