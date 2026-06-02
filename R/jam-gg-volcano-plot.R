## Todo
# - include caption text outside the plot (somehow)
#   for total points (bottom-left), and cutoffs (bottom-right)

#' Draw a volcano plot using ggplot2 with reasonable default arguments.
#'
#' Draw a volcano plot using ggplot2 with reasonable default
#' arguments, and with a large number of customization options.
#' The default plot uses smooth scatter plot for much improved
#' display of point density.
#' 
#' See `volcano_plot()` for the detailed customization options.
#' 
#' Note that ggplot2 does not currently display block arrows with
#' the number of hits that met statistical thresholds.
#' 
#' ## To label genes
#' 
#' * Define `hi_points` for the points to be highlighted.
#' 
#'    * It can be a `list` of `character` vectors, which are
#'    matched with gene values in `gene_colname`.
#'    The points are colored by `hi_colors` which should be
#'    a `character` vector of colors the same length as
#'    `hi_points` list length, optionally named using
#'    `names(hi_points)`,
#'    * It can be a `character` vector of gene values.
#'    This input is treated as a `list` containing one
#'    `character` vector.
#'    * **The color legend is not yet implemented**
#'    for `ggvolcano_plot()`.
#' 
#' * Define `hi_do_label` with `logical` vector, `character`
#' vector, or `integer` index.
#' 
#'    * Labels use the column matched by `label_colname`, which
#'    by default is the same column as `gene_colname`.
#'    * The label color is inherited from `hi_colors` except it
#'    it made much lighter, and desaturated, using
#'    `jamba::makeColorDarker(..., darkFactor=-1.6, sFactor=-1.5)`
#' 
#' @family jam plot functions
#' 
#' @returns `ggplot2` graphics object, where the `data.frame`
#'    is obtained with `gg@data`.
#' 
#' @param x `data.frame` with data suitable for volcano plot.
#' @param base_size `numeric` base font point size used for the
#'    ggplot2 theme `colorjam::theme_jam()`, default 20 intended
#'    for presentation or publication figures.
#' @param panel.grid.major.colour,panel.grid.major.colour `character`
#'    colors used for the major and minor grid lines.
#' @param axis.text.x.angle `numeric` angle used for the x-axis
#'    label text, default 0.
#' @param titleCex `numeric` title font adjustment
#' @param detail_factor `numeric` adjustment to increase overall
#'    detail in the density plot.
#' @param nbin `integer` number of bins to use for density plot.
#' @param bw_factor `numeric` adjustment to increase detail used
#'    in calculating the detail in the point density.
#' @param transFactor `numeric` exponent value used to adjust the
#'    point density, similar to that used by `smoothScatter()`
#'    in the form `function(x)x^transFactor`.
#' @param smooth_colors `character` vector of colors, or single color
#'    or single color ramp name recognized by `jamba::getColorRamp()`.
#' @param ... additional arguments are passed to `volcano_plot()`.
#' 
#' @examples
#' n <- 15000;
#' set.seed(12);
#' x_lfc <- (rnorm(n) * 1);
#' x_lfc <- x_lfc^2 * sign(x_lfc);
#' x_lfc <- x_lfc[order(-abs(x_lfc) + rnorm(n) / 2)];
#' x_pv <- sort(10^-(rnorm(n)*1.5)^2);
#' x <- data.frame(
#'    Gene=paste("gene", seq_len(n)),
#'    `log2fold Group-Control`=x_lfc,
#'    `P.Value Group-Control`=x_pv[order(-abs(x_lfc))],
#'    `mgm Group-Contol`=((rnorm(1500)+5)^2)/5,
#'    check.names=FALSE);
#'
#' ggvdf <- ggvolcano_plot(x);
#' plot(ggvdf)
#' 
#' @export
ggvolcano_plot <- function
(x,
 base_size=20,
 fold_cutoff=1.5,
 sig_cutoff=0.05,
 cutoff_color="#00000077",
 cutoff_linetype="dashed",
 cutoff_linewidth=1,
 fold_style=c('fold',
    'log2fold',
    'log2',
    'lfc'),
 baseColor="white",
 panel.grid.major.colour="grey90",
 panel.grid.minor.colour="transparent",
 axis.text.x.angle=60,
 titleCex=1,
 detail_factor=1,
 nbin=300,
 bw_factor=1,
 transFactor=0.24,
 smooth_colors=c("transparent",
    "lightblue",
    "lightskyblue3",
    "royalblue",
    "darkblue",
    "orange",
    "darkorange1",
    "orangered2"),
 aspect_adjust=1,
 verbose=FALSE,
 ...)
{
   #
   fold_style <- match.arg(fold_style);
   #
   if (verbose) {
      jamba::printDebug("ggvolcano_plot(): ",
         "Calling volcano_plot().");
   }
   vdf <- volcano_plot(x,
      do_plot=FALSE,
      transFactor=transFactor,
      ...)
   
   # h values for MASS::kde2()
   #
   if (verbose) {
      jamba::printDebug("ggvolcano_plot(): ",
         "Preparing ggplot2 object.");
   }

   xlim <- range(vdf$x)
   ylim <- range(vdf$y)
   ylim <- round(digits=2,
      ylim);
   bw_factor <- rep(bw_factor,
      length.out=2);
   hx <- diff(range(xlim, na.rm=TRUE)) / (80 * bw_factor[1] * (1/1.33) * aspect_adjust);
   # hx <- diff(range(xlim, na.rm=TRUE)) / (80 * bw_factor[1] * 6/8);
   hy <- diff(range(ylim, na.rm=TRUE)) / (80 * bw_factor[2] * 1);
   
   #
   x_exp <- xlim + diff(range(xlim)) * c(-0.6, 0.6);
   y_exp <- ylim + diff(range(ylim)) * c(-0.6, 0.6);
   x_exp <- xlim + diff(range(xlim)) * c(-0.3, 0.3);
   y_exp <- ylim + diff(range(ylim)) * c(-0.3, 0.3);
      
   nrow_x <- nrow(vdf)

   argsList <- attr(vdf, "argsList")
   use_xlab <- argsList$xlab;
   use_ylab <- argsList$ylab;

   p <- ggplot2::ggplot(vdf, ggplot2::aes(x=x, y=y)) +
      ggplot2::stat_density_2d(
         geom="raster",
         # data=subset(jp2tall, !outlier),
         ggplot2::aes(
            fill=(ggplot2::after_stat(count) / nrow_x)^transFactor),
         n=round(
            c(nbin * 2 * detail_factor, 
               nbin * detail_factor)),
         adjust=c(1, 1)/2,
         show.legend=FALSE,
         h=c(hx, hy * 2) * 1/2,   # slightly smaller point density
         # h=c(hx, hy * 2) * 2/3, # slightly larger point density
         contour=FALSE) +
      ggplot2::lims(
         x=xlim) +
         # y=ylim) +
      ggplot2::xlab(use_xlab) +
      ggplot2::ylab(use_ylab)

   # adjust axis expansion
   # x-axis labels
   if ("fold" %in% fold_style) {
      # convert to normal space
      p <- p +
         ggplot2::scale_x_continuous(
            # breaks=scales::breaks_width(1),
            breaks=function(y) {
               sort(unique(c(
                  scales::breaks_width(1)(y),
                  log2(fold_cutoff) * c(-1, 1)
               )))
            },
            labels=function(y){
               ifelse(y < 0,
                  (-1 * 2^(abs(y))),
                  (2^(abs(y))))
            },
            expand=ggplot2::expansion(mult=0.01))
   } else {
      p <- p +
      ggplot2::scale_x_continuous(
         # breaks=scales::breaks_width(1),
         breaks=function(y) {
            sort(unique(c(
               scales::breaks_width(1)(y),
               log2(fold_cutoff) * c(-1, 1)
            )))
         },
         labels=function(y){
            sapply(y, format)
            # ifelse(y < 0,
            #    (-1 * 2^(abs(y))),
            #    (2^(abs(y))))
         },
      expand=ggplot2::expansion(mult=0.01))
   }
   ## y-axis labels
   p <- p +
      ggplot2::scale_y_continuous(
         limits=ylim,
         breaks=function(y) {
            sort(unique(c(
               scales::breaks_width(1)(y),
               -log10(sig_cutoff)
            )))
         },
         # breaks=scales::breaks_width(1),
         labels=function(y){
            ifelse(y <= 2 | y == -log10(sig_cutoff),
               10^-y,
               scales::label_math(10^{- .x})(y))
         },
         expand=ggplot2::expansion(mult=0.01))

   ## fold cutoff ablines
   if (length(fold_cutoff) == 1 && fold_cutoff > 1) {
      p <- p +
         ggplot2::geom_vline(
            xintercept=log2(fold_cutoff) * c(-1, 1),
            linewidth=cutoff_linewidth,
            linetype=cutoff_linetype,
            color=cutoff_color)
   }
   ## significance cutoff abline
   if (length(sig_cutoff) == 1 && sig_cutoff < 1) {
      p <- p +
         ggplot2::geom_hline(
            yintercept=-log10(sig_cutoff),
            linetype=cutoff_linetype,
            linewidth=cutoff_linewidth,
            color=cutoff_color)
   }

   ## Set stable aspect ratio to accomodate xlim,ylim
   exp_aspect <- diff(x_exp) / diff(y_exp);
   p <- p +
      ggplot2::coord_fixed(exp_aspect / 1.25 * aspect_adjust)


   # FIll background solid?
   # if (fillBackground) {
   #    p <- p +
   #       ggplot2::geom_rect(
   #          fill=baseColor,
   #          # data=subset(jp2tall, !outlier & !duplicated(name)),
   #          alpha=1,
   #          xmin=x_exp[1],
   #          xmax=x_exp[2],
   #          ymin=y_exp[1],
   #          ymax=y_exp[2])
   # }

   # Plot theme
   p <- p + colorjam::theme_jam(
      # strip.background.fill=strip_bg,
      strip.text.size=ggplot2::rel(0.6 * titleCex),
      panel.grid.major.colour=panel.grid.major.colour,
      panel.grid.minor.colour=panel.grid.minor.colour,
      axis.text.x.angle=axis.text.x.angle,
      base_size=base_size);

   # Density colors
   if (length(smooth_colors)) {
      smooth_colors <- jamba::getColorRamp(smooth_colors,
         ...)
   }
   p <- p +
      ggplot2::scale_fill_gradientn(colours=smooth_colors);

   # optional hi points
   if (any(nchar(vdf$hi_points) > 0)) {
      # add non-overlapping labels
      pt_vdf <- subset(vdf, hi_points > 0);
      color_set <- attr(vdf, "color_set")
      border_set <- attr(vdf, "border_set")
      # Todo: Add second color scale,
      # look up proper gg method to do that
      p <- p + 
         ggplot2::geom_point(
            ggplot2::aes(
               # bg=point_type,
               # color=point_type,
               label=gene_values),
            bg=color_set[pt_vdf$point_type],
            color=border_set[pt_vdf$point_type],
            shape=21,
            size=2.5,
            # show.legend=FALSE,
            data=pt_vdf) +
         ggplot2::scale_color_identity()
   }

   # optional hi labels
   if (any(nchar(vdf$hi_points_label_shown) > 0)) {
      # add non-overlapping labels
      ptl_vdf <- subset(vdf, hi_points_label_shown > 0);
      color_set_lt <- jamba::makeColorDarker(color_set,
         darkFactor=-1.6,
         sFactor=-1.5)
      p <- p + ggrepel::geom_label_repel(
         ggplot2::aes(label=gene_values),
         bg=color_set_lt[ptl_vdf$point_type],
         color=jamba::setTextContrastColor(
            color_set_lt[ptl_vdf$point_type]),
         # border=border_set[ptl_vdf$point_type],
         data=ptl_vdf)
      p
   }
   p

}

# vdf <- volcano_plot(x,
#    hi_points=sample(subset(x, abs(x[,2]) > 1 & x[,3] < 0.01)[,1], 200),
#    hi_do_labels=TRUE, sig_max_range=1e-50, fold_max_range=200, max_labels=10, blockarrow=FALSE)

