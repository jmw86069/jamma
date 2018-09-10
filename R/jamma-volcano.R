#' Create Volcano plot of a statistical comparison
#'
#' Create Volcano plot of a statistical comparison
#'
#' This function takes a data.frame of statistical results in the form
#' of contrasts with fold changes and P-values, and produces a volcano
#' plot with -log10 P-value on the y-axis, versus log2 fold change on
#' the x-axis. The primary purpose of the plot is to observe the distribution
#' of up- and down- changes, the distribution of P-values, relative to
#' optionally supplied statistical cutoffs used to define statistical hits.
#'
#' To do:
#'
#' \itemize{
#'    \item{Adjust y-axis margin so the y-axis label is displayed, since it
#'       indicates on the plot which P-value column is being displayed.}
#'    \item{Aesthetics, Increase default bandwidthN. Increase default
#'       ylab and xlab font sizes. Consider higher contrast color for the
#'       y-axis statistical hits block arrow.}
#'    \item{Include hits with returned data, it should be easy to retrieve
#'       the exact entries marked as hits on the plot.}
#' }
#'
#' Note: This function is not exported, as it is still being ported
#' to the Jam way of doing things.
#'
volcanoPlot <- function
(x,
 n=NULL, coef=NULL,

 ## Data columns
 geneColumn=head(provigrep(geneColGrep, colnames(x)),1),
 lfcCol=head(provigrep(lfcColGrep, colnames(x)),1),
 pvCol=head(provigrep(pvColGrep, colnames(x)),1),
 intensityCols=provigrep(intensityColGrep, colnames(x)),

 ## grep expressions
 geneColGrep=c("^gene","symbol"),
 lfcColGrep=c("log.{0,2}fc", "lfc", "log.{0,4}fold", "fold.*ch",
   "fc", "fold", "log.*ratio", "ratio"),
 pvColGrep=c("^padj|adj.*p","fdr","q.*val","p.*val"),
 intensityColGrep=c("maxgroupmean|maxmean","groupMean","aveexpr","^fpkm"),

 ## Stats cutoffs
 pvalueCutoff=0.01,
 lfcCutoff=0,
 intensityCutoff=NULL,
 minPvalue=1e-20,
 pvalueFloor=minPvalue,
 maxLog2FC=10,
 xRange=NULL,
 yRange=NULL,
 lfcCeiling=NULL,

 ## Plot features, title, axis labels, margins
 margins=c(5,6,4,2),
 main="Volcano Plot",
 lineMain=3,
 subMain=NULL,
 font.main=1,
 symmetricAxes=TRUE,
 xlab=paste("Fold change (", lfcCol, ")", sep=""),
 ylab=paste("Significance (", pvCol, ")", sep=""),
 nXlabels=13,
 nYlabels=7,
 ablineColor="#000000AA",
 orientation="up",
 addPlot=FALSE,

 ## Plot subtitle with stats cutoffs applied
 doStatsSubtitles=TRUE,
 statsLine=4,
 statsCex=0.8,

 ## Point sizes and colors
 ptCex=0.9,
 ptPch=21,
 hitType="hits",
 hitColorSet=c("#77777777", "#99000055"),
 pointColors=NULL,
 pointBgColors=NULL,
 pointColorSet=c("base"="#77777733",
   "up"="#99000055",
   "down"="#00009955",
   "base.highlight"="#FFDD55DD",
   "up.highlight"="#FFDD55DD",
   "down.highlight"="#FFDD55DD"),
 pointBorderSet=c("base"="#33333333",
   "up"="#77000088",
   "down"="#00007788",
   "base.highlight"="#773322DD",
   "up.highlight"="#773322DD",
   "down.highlight"="#773322DD"),
 contrastFactor=1.5,

 ## Highlighted points
 highlightPoints=NULL,
 highlightHits=FALSE,
 highlightItemType=paste("highlighted", hitType, sep=" "),
 highlightCex=1,

 ## Non-overlapping label
 pointsToLabel=NULL,
 noObscureMode=TRUE,
 initialAngle=0,
 initialRadius=0,
 labelCex=1,
 labelMar=labelCex*1.1,
 boxColor="#EEBB7788",
 boxBorderColor="#000000AA",
 labelFixedCoords=NULL,
 doHighlightLabels=TRUE,
 sectionLabelSpacing=12,
 labelN=30,

 ## Jitter points to make overlaps more visibly clear
 doJitter=10,
 jitterAmount=NULL,

 ## Logic to group probes by gene before counting the number of hits
 hitsByGene=FALSE,
 printMultiGeneHits=FALSE,

 ## smooth scatter plotting
 doSmoothScatter=TRUE,
 smoothScatterFunc=plotSmoothScatter,
 colramp=c("white", "lightblue", "lightskyblue3", "royalblue",
   "darkblue", "orange", "darkorange1", "orangered2"),
 useRaster=TRUE,
 transFactor=0.24,
 transformation=function(x){x^transFactor}, # use 0.14 for very large datasets

 doBoth=FALSE,
 bothColorOnlyHits=TRUE,
 pointSubset=NULL,
 plotOnlySubset=FALSE,
 nbin=256,
 nrpoints=0,
 ## Block arrow parameters
 doBlockArrows=TRUE,
 blockArrowFont=1,
 blockArrowHitColor="#E67739FF", blockArrowUpColor="#990000FF", blockArrowDownColor="#000099FF",
 blockArrowLabelHitColor="#FFFFAA", blockArrowLabelUpColor="#FFFFAA", blockArrowLabelDownColor="#FFFFAA",
 blockArrowCex=c(1,1), blockArrowLabelCex=1,
 doShadowTextArrows=TRUE,
 ## Top histogram
 doTopHist=FALSE, doTopHistCutoffs=c("pvalue", "foldchange"),
 topHistBreaks=100, topHistColor="#000099FF",
 topHistCexSub=1, topHistSubLine=0,
 topHistPlotFraction=1/5, topHistBy=0.20,
 ## Label font sizing
 cexSub=0.8, cexMain=1.5, cex.axis=1.2,
 lineXlab=2.5, lineYlab=4,
 labelHits=!doBlockArrows, verbose=FALSE,
 ...)
{
   ## Purpose is to wrapper a simple volcano plot with log-friendly axis labels
   ##
   ## orientation can be "up" for upright, "right" to tilt 90 degree on its side
   ##
   ## minPvalue is designed to set NA or '0' values to a fixed value
   ## pvalueFloor is designed to set extremely low P-values (1e-150) to a
   ##    minimum value to control the plot boundaries
   ## sectionLabelSpacing is the fraction of the y-axis range away from each
   ##    border used to position the highlight hit counts labels
   ##    (if highlightPoints is not NULL)
   ##
   ## The graph plots fold change on the x-axis, but with log2 scale,
   ##    the numbers are technically being reported in normal space.
   ##    However, to give it a log2 axis label:
   ##    xlab=expression(log[2]*"(Fold change)"
   ##

   ## Handle input object type
   if (igrepHas("MArrayLM", class(x))) {
      ## MArrayLM is direct output from lmFit() or eBayes()
      if (verbose) {
         printDebug("volcanoPlot() applying topTable(x)");
      }
      if (length(coef) == 0) {
         coef <- 1;
      }
      x <- topTable(x, coef=coef, number=Inf);
      if ((length(geneColumn) == 0 || !geneColumn %in% colnames(x)) &&
         !is.null(rownames(x))) {
         x[,geneColumn] <- rownames(x);
      }
   }

   ## Validate parameters
   blockArrowCex <- rep(blockArrowCex, length.out=2);

   ## Define point colors
   ## Only alter parameters which were provided in the function call,
   ## but keep default values which were not defined
   pointColorSet <- updateFunctionParamList(functionName="volcanoPlot",
      paramName="pointColorSet", newValues=pointColorSet);
   pointBorderSet <- updateFunctionParamList(functionName="volcanoPlot",
      paramName="pointBorderSet", newValues=pointBorderSet);
   if (verbose) {
      printDebug(" pointColorSet:");
      printDebug((pointColorSet), fgText=list(c(setTextContrastColor(pointColorSet))),
         bgText=list(pointColorSet));
      printDebug("pointBorderSet:");
      printDebug((pointBorderSet), fgText=list(c(setTextContrastColor(pointBorderSet))),
         bgText=list(pointBorderSet));
   }

   ##
   if (doBoth && bothColorOnlyHits) {
      hitColorSet[1] <- "#FFFFFF00";
      pointColorSet["base"] <- hitColorSet[1];
   }

   if (verbose) {
      printDebug("lfcCol: ", lfcCol, c("orange", "lightblue"));
      printDebug("pvCol: ", pvCol, c("orange", "lightblue"));
   }
   origPar <- par();

   ## Allow processing only a subset 'n' rows of data
   if (is.null(n)) {
      n <- nrow(x);
   } else if (n < nrow(x)) {
      x <- x[1:n,,drop=FALSE];
   }

   if (is.null(rownames(x)) || length(tcount(rownames(x), minCount=2)) > 0) {
      rownames(x) <- makeNames(rep("row", n));
   }

   ## Optionally use intensity/expression level as a filter for hits
   if (!is.null(intensityCutoff)) {
      intensityCols <- intersect(colnames(x), intensityCols);
      if (length(intensityCols) > 0) {
         if (verbose) {
            printDebug("Applying intensityCutoff.");
         }
         if (!suppressPackageStartupMessages(require(matrixStats, quietly=TRUE))) {
            stop("The matrixStats package is required.");
         }
         intensityValues <- rowMaxs(as.matrix(x[,intensityCols,drop=FALSE]));
         metIntensity <- intensityValues > intensityCutoff;
         if (verbose) {
            print(table(metIntensity));
         }
      } else {
         metIntensity <- rep(TRUE, nrow(x));
      }
   } else {
      metIntensity <- rep(TRUE, nrow(x));
   }
   pvHitsUpWhich <- which(x[,pvCol] <= pvalueCutoff &
      (x[,lfcCol]) >= lfcCutoff & metIntensity);
   pvHitsDownWhich <- which(x[,pvCol] <= pvalueCutoff &
      -(x[,lfcCol]) >= lfcCutoff & metIntensity);
   pvHitsBothWhich <- which(x[,pvCol] <= pvalueCutoff &
      abs(x[,lfcCol]) >= lfcCutoff & metIntensity);
   pvHitsUp <- length(pvHitsUpWhich);
   pvHitsDown <- length(pvHitsDownWhich);
   pvHits <- length(pvHitsBothWhich);
   upHits <- rownames(x)[pvHitsUpWhich];
   downHits <- rownames(x)[pvHitsDownWhich];

   ## Define point colors
   if (is.null(pointColors)) {
      pvColors1 <- nameVector(rep(pointColorSet["base"], n), rownames(x));
      pvColors1[upHits] <- pointColorSet["up"];
      pvColors1[downHits] <- pointColorSet["down"];
   } else {
      pvColors1 <- pointColors;
   }
   if (is.null(pointBgColors)) {
      #pvBgColors1 <- makeColorDarker(pvColors1, darkFactor=1.5);
      pvBgColors1 <- nameVector(rep(pointBorderSet["base"], n), rownames(x));
      pvBgColors1[upHits] <- pointBorderSet["up"];
      pvBgColors1[downHits] <- pointBorderSet["down"];
   } else {
      pvBgColors1 <- pointBgColors;
   }
   if (verbose) {
      printDebug("head(pvColors1): ", c("orange", "lightblue"));
      print(head(pvColors1));
   }
   if (class(x[,pvCol]) %in% c("numeric", "integer")) {
      pvValues <- nameVector(x[,pvCol], rownames(x)[1:n]);
   } else {
      pvValues <- nameVector(as.numeric(x[,pvCol]), rownames(x)[1:n]);
   }
   pvValues[is.na(pvValues)] <- 1;
   if (verbose) {
      printDebug("minPvalue: ", minPvalue, c("orange", "lightblue"));
      printDebug("pvValues: ", head(pvValues), c("orange", "lightblue"));
      printDebug("any(pvValues < minPvalue): ", any(pvValues < minPvalue), c("orange", "lightblue"));
   }
   if (!is.null(minPvalue) && any(pvValues < minPvalue)) {
      pvValues[pvValues < minPvalue] <- minPvalue;
   }
   if (!is.null(pvalueFloor) && any(pvValues < pvalueFloor)) {
      pvValues[pvValues < pvalueFloor] <- pvalueFloor;
   }
   yValues1 <- -log10(pvValues);

   fcValues <- nameVector(x[,lfcCol], rownames(x));

   fcValues[is.na(fcValues)] <- 0;
   fcValues[abs(fcValues) > maxLog2FC] <- sign(fcValues[abs(fcValues) > maxLog2FC]) * maxLog2FC;
   xValues1 <- fcValues;

   ## Allow highlighting a subset of points, which overrides the pointSubset argument
   if (is.null(highlightPoints) && highlightHits) {
      if (verbose) {
         printDebug("   volcanoPlot defining highlightPoints as up/down hits.");
      }
      highlightPoints <- rownames(x)[unique(c(pvHitsUpWhich, pvHitsDownWhich))];
      highlightItemType <- hitType;
   }
   if (!is.null(pointsToLabel)) {
      row2name <- c(nameVector(rownames(x), x[,geneColumn]), nameVector(rownames(x)));
      pointsToLabel <- flipVector(row2name[pointsToLabel], makeNamesFunc=c);
   }
   if (verbose) {
      printDebug("   volcanoPlot length(highlightPoints):",
         format(length(highlightPoints), big.mark=","));
   }
   if (verbose) {
      printDebug("   doBoth=", doBoth);
   }
   if (!is.null(highlightPoints)) {
      if (verbose) {
         printDebug("Creating subset of points for highlighting.");
      }
      if (is.null(doSmoothScatter) || doSmoothScatter) {
         if (verbose) {
            printDebug("Setting doBoth=", "TRUE");
         }
         doBoth <- TRUE;
      }
      pointSubset <- which(rownames(x) %in% highlightPoints | x[,geneColumn] %in% highlightPoints);
      if (verbose) {
         printDebug("   volcanoPlot head(pointSubset):", head(pointSubset));
         printDebug("   volcanoPlot length(pointSubset):",
            format(length(pointSubset), big.mark=","));
      }
      ## If we're plotting all individual points, then we can color-code the hits, and the highlighted
      ## items all independently
      highlightUpHits <- intersect(rownames(x)[pointSubset], upHits);
      highlightDownHits <- intersect(rownames(x)[pointSubset], downHits);
      highlightBaseHits <- setdiff(rownames(x)[pointSubset], c(upHits, downHits));
      pvColors1[highlightBaseHits] <- pointColorSet["base.highlight"];
      pvColors1[highlightUpHits] <- pointColorSet["up.highlight"];
      pvColors1[highlightDownHits] <- pointColorSet["down.highlight"];
      pvBgColors1[highlightBaseHits] <- pointBorderSet["base.highlight"];
      pvBgColors1[highlightUpHits] <- pointBorderSet["up.highlight"];
      pvBgColors1[highlightDownHits] <- pointBorderSet["down.highlight"];

      pointSubsetLabels <- x[pointSubset,geneColumn];
      if (verbose) {
         printDebug("highlightUpHits: ", head(highlightUpHits), c("orange", "lightblue"));
         printDebug("highlightDownHits: ", head(highlightDownHits), c("orange", "lightblue"));
         printDebug("highlightBaseHits: ", head(highlightBaseHits), c("orange", "lightblue"));
         printDebug("pointSubsetLabels: ", head(pointSubsetLabels), c("orange", "lightblue"));
         printDebug("pvColors1[highlightDownHits]:");
         printDebug(names(pvColors1[highlightDownHits]), list(pvColors1[highlightDownHits]));
         printDebug("pvBgColors1[highlightDownHits]:");
         printDebug(names(pvBgColors1[highlightDownHits]), list(pvBgColors1[highlightDownHits]));
      }
      names(pointSubset) <- names(xValues1)[pointSubset];
      if (length(highlightCex) == 1) {
         pointSubsetCex <- highlightCex;
      } else {
         ## Order the given highlightCex alongside the matching rows in the data
         names(highlightCex) <- highlightPoints;
         pointSubsetCex <- highlightCex[names(pointSubset)];
      }
      if (length(pointSubset) == 0) {
         pointSubset <- NULL;
         pointSubsetCex <- 0;
         if (verbose) {
            printDebug("No points found to highlight.", fgText="yellow");
         }
      } else {
         ## Generate hit counts among the highlighted points
         pvHighlightHitsUp <- length(which(x[1:n,pvCol][pointSubset] <= pvalueCutoff & (x[1:n,lfcCol][pointSubset]) >= lfcCutoff));
         pvHighlightHitsDown <- length(which(x[1:n,pvCol][pointSubset] <= pvalueCutoff & -(x[1:n,lfcCol][pointSubset]) >= lfcCutoff));
         pvHighlightHitsOther <- length(pointSubset) - (pvHighlightHitsUp + pvHighlightHitsDown);
         if (verbose) {
            printDebug("   pvHighlightHitsUp: ", pvHighlightHitsUp, c("orange", "lightblue"));
            printDebug("   pvHighlightHitsDown: ", pvHighlightHitsDown, c("orange", "lightblue"));
            printDebug("   pvHighlightHitsOther: ", pvHighlightHitsOther, c("orange", "lightblue"));
         }
      }
   }
   if (!is.null(pointSubset) && plotOnlySubset) {
      xValues1 <- xValues1[pointSubset];
      yValues1 <- yValues1[pointSubset];
      pvColors1 <- pvColors1[pointSubset];
      pvBgColors1 <- pvBgColors1[pointSubset];
      pointSubset <- NULL;
      pointSubsetCex <- 0;
   }

   ## Switch pvBgColors1 and pvColors1
   pvColors2 <- pvColors1;
   pvColors1 <- pvBgColors1;
   pvBgColors1 <- pvColors2;

   if (!is.null(lfcCeiling)) {
      xValues1[abs(xValues1) > lfcCeiling] <- sign(xValues1[abs(xValues1) > lfcCeiling]) * lfcCeiling;
   }

   if (doJitter) {
      if (verbose) {
         printDebug("Jittering points by:", format(digits=2, jitterAmount));
      }
      xValues2 <- xValues1;
      yValues2 <- yValues1;
      xValues1 <- jitter2(xValues1, factor=as.numeric(doJitter), amount=jitterAmount);
      yValues1 <- jitter2(yValues1, factor=as.numeric(doJitter), amount=jitterAmount);
   }

   if (is.null(xRange)) {
      xRange <- range(xValues1);
   }
   if (is.null(yRange)) {
      yRange <- range(yValues1);
   }
   if (symmetricAxes) {
      xRange <- c(-max(abs(xRange)), max(abs(xRange)));
   }

   if (!is.null(pointSubset)) {
      ###############################################
      ## Position labels in quadrants of the plot
      namesX <- c(xRange[1] - xRange[1]/sectionLabelSpacing, xRange[2] - xRange[2]/sectionLabelSpacing);
      namesY <- rep(yRange[2] - yRange[2]/sectionLabelSpacing, 2);
      namesLabel <- paste(c(format(big.mark=",", pvHighlightHitsDown), format(big.mark=",", pvHighlightHitsUp)), highlightItemType, sep=" ");
      namesLabel <- gsub("^(1 .+)s$", "\\1", namesLabel);
      if (verbose) {
         printDebug("namesLabel: ", namesLabel, c("orange", "seagreen"));
      }
   }

   parList <- list();
   parList[["prePlots"]] <- par();

   ## labelCoords will have the return data from addNonOverlappingLabels() but only
   ## if we end up calling that method
   labelCoords <- NULL;
   ##########################################################
   ## Histogram along the top border, set up the details here
   if (doTopHist) {
      if (length(grep("pval", doTopHistCutoffs)) > 0) {
         ## Apply P-value filtering
         pvHitsWhich <- which(x[1:n,pvCol] <= pvalueCutoff);
      } else {
         pvHitsWhich <- 1:n;
      }
      if (length(grep("fc|fold", doTopHistCutoffs)) > 0) {
         ## Apply fold change filtering
         fcHitsWhich <- which(abs(x[1:n,lfcCol]) >= lfcCutoff);
      } else {
         fcHitsWhich <- 1:n;
      }
      histWhich <- intersect(pvHitsWhich, fcHitsWhich);
      #return(histWhich);

      ## Color the bars consistent with the block arrows
      topBarColors <- c("#000099FF", "#990000FF");
      parMar <- par("mar");
      topHistBreaks <- seq(from=xRange[1], to=xRange[2], by=topHistBy);
      if (max(topHistBreaks) < xRange[2]) {
         if (verbose) {
            printDebug("Fixing topHistBreaks: ", cPaste(topHistBreaks), c("orange", "lightblue"));
         }
         topHistBreaks <- c(topHistBreaks, xRange[2]);
      }
      if (orientation %in% c("right")) {
         topHist <- hist(xValues2[histWhich], breaks=topHistBreaks, plot=FALSE);
         plotZones <- matrix(c(2,1), nrow=1);
         layout(plotZones, widths=c(1-topHistPlotFraction, topHistPlotFraction), heights=c(1));
         ## Change margins, then plot the top histogram
         ## default is par(mar=c(5.1, 4.1, 4.1, 2.1));
         #par(mar=c(5.1, 0, 4.1, 2.1));
         par("mar"=c(parMar[1], 0, parMar[3], parMar[4]));
         parList[["preTopHist"]] <- par();
         barplot(topHist$counts, axes=FALSE, xlim=c(0, max(topHist$counts)), ylim=yRange, space=0, col=topHistColor, horiz=TRUE);
         parList[["postTopHist"]] <- par();
         #par(mar=c(5.1, 4.1, 4.1, 2.1));
         par("mar"=c(parMar[1], parMar[2], parMar[3], parMar[4]));
      } else {
         topHist <- hist(xValues2[histWhich], breaks=topHistBreaks, plot=FALSE);
         plotZones <- matrix(c(1,2), nrow=2);
         layout(plotZones, widths=c(1), heights=c(topHistPlotFraction, 1-topHistPlotFraction));
         ## Change margins, then plot the top histogram
         ## default is par(mar=c(5.1, 4.1, 4.1, 2.1));
         par("mar"=c(1, parMar[2], parMar[3], parMar[4]));
         parList[["preTopHist"]] <- par();
         if (verbose) {
            printDebug("par(mgp): ", par("mgp"));
            printDebug("length(topHistBreaks): ", length(topHistBreaks));
         }
         r1 <- as.integer(length(topHistBreaks)/2);
         topHistCol <- rep(topBarColors, c(r1, r1+1));
         barplot(topHist$counts, axes=FALSE, ylim=c(0, max(topHist$counts)),
            space=0, col=topHistCol, border=makeColorDarker(topHistCol),
            horiz=FALSE, las=2, cex.axis=1);
         prettyAt1 <- pretty(c(0, max(topHist$counts)), n=4);
         prettyAt1 <- prettyAt1[prettyAt1 <= max(topHist$counts)];
         axis(2, las=2, cex.axis=1.3, at=prettyAt1, ...);
         title(xlab=paste("Distribution of ", hitType, sep=""),
            cex.lab=topHistCexSub, xpd=TRUE, outer=FALSE, line=topHistSubLine);
         parList[["postTopHist"]] <- par();
         par("mar"=c(parMar[1], parMar[2], parMar[3]-1.5, parMar[4]));
      }
   }

   ######################
   ## Smooth scatter plot
   if (doSmoothScatter || doBoth) {
      if (orientation %in% c("right")) {
         ## Volcano plot on its side, tiled 90 degrees to the right
         smoothScatterFunc(y=xValues1, x=yValues1,
            xaxt="n", yaxt="n", xlab="", ylab="", transformation=transformation,
            useRaster=useRaster, xlim=yRange, ylim=xRange, nbin=nbin,
            colramp=colramp, add=addPlot, nrpoints=nrpoints, ...);
         title(xlab=ylab, line=lineXlab, cex.axis=cex.axis, ...);
         title(ylab=xlab, line=lineYlab, cex.axis=cex.axis, ...);
         parUsr <- par("usr");
         if (doBoth) {
            if (!is.null(pointSubset)) {
               points(y=xValues1[pointSubset], x=yValues1[pointSubset], pch=ptPch, cex=pointSubsetCex,#ptCex,
                  col=pvColors1[pointSubset], bg=pvBgColors1[pointSubset], xaxt="n", yaxt="n",
                  xlim=yRange, ylim=xRange, ...);
               ## Optionally labels the highlighted points, using the addNonOverlappingLabels() function
               if (doHighlightLabels) {
                  labelCoords <- addNonOverlappingLabels(y=xValues1[pointSubset], x=yValues1[pointSubset],
                     labelCex=labelCex, labelMar=labelMar, boxColor=boxColor,
                     boxBorderColor=boxBorderColor, initialRadius=initialRadius,
                     txt=pointSubsetLabels, n=labelN, useNewCex=TRUE, initialAngle=initialAngle,
                     fixedCoords=labelFixedCoords,
                     ...);
               }
               ## Add text labels indicating the number of highlighted points
               textAdjX <- c(1, 1);
               if (labelHits) {
                  for(i in seq_along(namesLabel)) {
                     ## Note the flipped x and y values
                     printDebug("   namesLabel[i]: ", namesLabel[i], c("orange", "seagreen"));
                     text(y=namesX[i], x=namesY[i], label=namesLabel[i], adj=c(textAdjX[i], 0.5));
                  }
               }
            } else {
               points(y=xValues1, x=yValues1, pch=ptPch, cex=ptCex,
                      col=pvColors1, bg=pvBgColors1, xaxt="n", yaxt="n", xlab=ylab, ylab=xlab,
                      xlim=yRange, ylim=xRange, ...);
            }
         }
      } else {
         ## Upright volcano plot (standard orientation)
         smoothScatterFunc(x=xValues1, y=yValues1, #col=pvColors1, bg=pvBgColors1,
              xaxt="n", yaxt="n", xlab="", ylab="", transformation=transformation, useRaster=useRaster,
              xlim=xRange, ylim=yRange, nbin=nbin, colramp=colramp, add=addPlot, nrpoints=nrpoints, ...);
         title(xlab=xlab, line=lineXlab, cex.axis=cex.axis, ...);
         title(ylab=ylab, line=lineYlab, cex.axis=cex.axis, ...);
         parUsr <- par("usr");
         if (verbose) {
            printDebug("xlab: ", xlab, c("orange", "lightblue"));
         }
         if (doBoth) {
            if (!is.null(pointSubset)) {
               if (verbose) {
                  printDebug("Following doBoth to label pointSubset.");
               }
               points(x=xValues1[pointSubset], y=yValues1[pointSubset], pch=ptPch, cex=pointSubsetCex,#ptCex,
                      col=pvColors1[pointSubset], bg=pvBgColors1[pointSubset], xaxt="n", yaxt="n",
                      xlim=xRange, ylim=yRange, ...);
               ## Optionally labels the highlighted points, using the addNonOverlappingLabels() function
               if (doHighlightLabels) {
                  printDebug("length(pointSubset):", length(pointSubset), ", labelN:", labelN);
                  labelCoords <- addNonOverlappingLabels(x=xValues1[pointSubset], y=yValues1[pointSubset], initialAngle=initialAngle,
                                          txt=pointSubsetLabels, n=labelN,
                                          labelCex=labelCex,  labelMar=labelMar, boxColor=boxColor,
                                          boxBorderColor=boxBorderColor, initialRadius=initialRadius,
                                          fixedCoords=labelFixedCoords,
                                          ...);
               }
               ## Add text labels indicating the number of highlighted points
               textAdjX <- c(0, 1);
               if (labelHits) {
                  for(i in seq_along(namesLabel)) {
                     printDebug("   namesLabel[i]: ", namesLabel[i], c("orange", "seagreen"));
                     text(x=namesX[i], y=namesY[i], label=namesLabel[i], adj=c(textAdjX[i], 0.5));
                  }
               }
            } else {
               if (verbose) {
                  printDebug("Following doBoth but pointSubset is NULL.");
               }
               points(x=xValues1, y=yValues1, pch=ptPch, cex=ptCex,
                      col=pvColors1, bg=pvBgColors1, xaxt="n", yaxt="n",
                      xlim=xRange, ylim=yRange, ...);
            }
         }
      }
      parList[["postSmoothScatter"]] <- par();
   }
   if (!doSmoothScatter || doBoth) {
      ## Non-smoothScatter
      ## Add points to an existing plot
      if (addPlot || doBoth) {
         if (orientation %in% c("right")) {
            if (!is.null(pointSubset)) {
               points(y=xValues1[pointSubset], x=yValues1[pointSubset], pch=ptPch, cex=pointSubsetCex,#ptCex,
                    col=pvColors1[pointSubset], bg=pvBgColors1[pointSubset], xaxt="n", yaxt="n",
                    xlim=yRange, ylim=xRange, ...);
            } else {
               points(y=xValues1, x=yValues1, pch=21, cex=ptCex,
                    col=pvColors1, bg=pvBgColors1, xaxt="n", yaxt="n",
                    xlim=yRange, ylim=xRange, ...);
            }
         } else {
            if (!is.null(pointSubset) && !plotOnlySubset) {
               points(y=yValues1[pointSubset], x=xValues1[pointSubset], pch=ptPch, cex=pointSubsetCex,#ptCex,
                  col=pvColors1[pointSubset], bg=pvBgColors1[pointSubset], xaxt="n", yaxt="n",
                  xlim=xRange, ylim=yRange, ...);
               if (verbose) {
                  printDebug("pvColors1[pointSubset]:");
                  printDebug((pvColors1[pointSubset]), list(pvBgColors1[pointSubset]));
                  printDebug("pvBgColors1[pointSubset]:");
                  printDebug((pvBgColors1[pointSubset]), list(pvBgColors1[pointSubset]));
               }
            } else {
               points(x=xValues1, y=yValues1, pch=ptPch, cex=ptCex,
                    col=pvColors1, bg=pvBgColors1, xaxt="n", yaxt="n",
                    xlim=xRange, ylim=yRange, ...);
            }
         }
      } else {
         ## Create new all-point scatter plot
         if (orientation %in% c("right")) {
            doPoints <- which(!names(xValues1) %in% pointSubset);
            plot(y=xValues1[doPoints], x=yValues1[doPoints], pch=ptPch, cex=ptCex,
                 col=pvColors1[doPoints], bg=pvBgColors1[doPoints], xaxt="n", yaxt="n",
                 xlab="", ylab="",
                 xlim=yRange, ylim=xRange, ...);
            title(xlab=ylab, line=lineXlab, cex.axis=cex.axis, ...);
            title(ylab=xlab, line=lineYlab, cex.axis=cex.axis, ...);
            ######################
            ## Borrowed from above
            if (!is.null(pointSubset)) {
               points(y=xValues1[pointSubset], x=yValues1[pointSubset], pch=ptPch, cex=pointSubsetCex,#ptCex,
                      col=pvColors1[pointSubset], bg=pvBgColors1[pointSubset], xaxt="n", yaxt="n",
                      xlim=yRange, ylim=xRange, ...);
               ## Optionally labels the highlighted points, using the addNonOverlappingLabels() function
               if (doHighlightLabels) {
                  labelCoords <- addNonOverlappingLabels(y=xValues1[pointSubset], x=yValues1[pointSubset],
                                          txt=pointSubsetLabels, n=labelN, useNewCex=TRUE,
                                          labelCex=labelCex, labelMar=labelMar, boxColor=boxColor,
                                          boxBorderColor=boxBorderColor, initialRadius=initialRadius,
                                          fixedCoords=labelFixedCoords,
                                          initialAngle=initialAngle, ...);
               }
               ## Add text labels indicating the number of highlighted points
               textAdjX <- c(1, 1);
               if (labelHits) {
                  for(i in seq_along(namesLabel)) {
                     ## Note the flipped x and y values
                     text(y=namesX[i], x=namesY[i], label=namesLabel[i], adj=c(textAdjX[i], 0.5));
                  }
               }
            } else {
               points(y=xValues1, x=yValues1, pch=ptPch, cex=ptCex,
                      col=pvColors1, bg=pvBgColors1, xaxt="n", yaxt="n", xlab=ylab, ylab=xlab,
                      xlim=yRange, ylim=xRange, ...);
            }
            #######
         } else {
            doPoints <- which(!names(xValues1) %in% pointSubset);
            if (verbose) {
               printDebug("length(doPoints): ", length(doPoints), c("orange", "lightblue"));
            }
            plot(x=xValues1[doPoints], y=yValues1[doPoints], pch=ptPch, cex=ptCex,
                 col=pvColors1[doPoints], bg=pvBgColors1[doPoints], xaxt="n", yaxt="n",
                 xlab="", ylab="",
                 xlim=xRange, ylim=yRange, ...);
            title(xlab=xlab, line=lineXlab, cex.axis=cex.axis, ...);
            title(ylab=ylab, line=lineYlab, cex.axis=cex.axis, ...);
            if (verbose) {
               printDebug("ping", "orange");
            }
            ######################
            ## Optionally label a specific subset of points
            if (!is.null(pointsToLabel) && length(pointsToLabel) > 0) {
               if (verbose) {
                  printDebug("Labeling ", length(pointsToLabel), " pointsToLabel.", c("orange", "lightblue"));
               }
               ## If pointsToLabel has names, they're used to match rownames(x), and the values
               ## are used as text labels on the plot.  This allows for adding a label which is
               ## different from the data matrix rownames.
               if (is.null(names(pointsToLabel))) {
                  pointsToLabel <- nameVector(pointsToLabel);
               }
               labelCoords <- addNonOverlappingLabels(x=xValues1[names(pointsToLabel)], y=yValues1[names(pointsToLabel)],
                                       labelCex=labelCex, labelMar=labelMar, boxColor=boxColor, boxBorderColor=boxBorderColor,
                                       initialRadius=initialRadius,
                                          txt=pointsToLabel, noObscureMode=noObscureMode, initialAngle=initialAngle,
                                          fixedCoords=labelFixedCoords);
            }
            ######################
            ## Borrowed from above
            if (!is.null(pointSubset)) {
               if (verbose) {
                  printDebug("length(pointSubset): ", length(pointSubset), c("orange", "lightblue"));
               }
               points(x=xValues1[pointSubset], y=yValues1[pointSubset], pch=ptPch, cex=pointSubsetCex,#ptCex,
                      col=pvColors1[pointSubset], bg=pvBgColors1[pointSubset], xaxt="n", yaxt="n",
                      xlim=xRange, ylim=yRange, cex.axis=cex.axis, ...);
               ## Optionally labels the highlighted points, using the addNonOverlappingLabels() function
               if (doHighlightLabels) {
                  labelCoords <- addNonOverlappingLabels(x=xValues1[pointSubset], y=yValues1[pointSubset],
                                          labelCex=labelCex, labelMar=labelMar, boxColor=boxColor, boxBorderColor=boxBorderColor,
                                          initialRadius=initialRadius,
                                          txt=pointSubsetLabels, n=labelN, initialAngle=initialAngle,
                                          fixedCoords=labelFixedCoords);
               }
               ## Add text labels indicating the number of highlighted points
               textAdjX <- c(0, 1);
               if (labelHits) {
                  for(i in seq_along(namesLabel)) {
                     text(x=namesX[i], y=namesY[i], label=namesLabel[i], adj=c(textAdjX[i], 0.5));
                  }
               }
            } else {
               points(x=xValues1, y=yValues1, pch=ptPch, cex=ptCex,
                      col=pvColors1, bg=pvBgColors1, xaxt="n", yaxt="n",
                      xlim=xRange, ylim=yRange, ...);
            }
            #######
         }
      }
      parList[["postScatter"]] <- par();
   }

   multiGenesUp <- character(0);
   multiGenesDown <- character(0);
   if (!addPlot) {
      #labelUp <- paste(format(big.mark=",", pvHitsUp), "up-regulated", hitType);
      #labelUp <- gsub("((^|[ ])1 .+)s$", "\\1", labelUp);
      #labelDown <- paste(format(big.mark=",", pvHitsDown), "down-regulated", hitType);
      #labelDown <- gsub("((^|[ ])1 .+)s$", "\\1", labelDown);

      ## Define labels
      labelUp <- paste(format(big.mark=",", pvHitsUp), hitType, "up");
      labelDown <- paste(format(big.mark=",", pvHitsDown), hitType, "down");

      ## Get rid of "1 genes", instead a more civilized "1 gene"
      labelUp <- gsub("((^|[ ])1 .+)s$", "\\1", labelUp);
      labelDown <- gsub("((^|[ ])1 .+)s$", "\\1", labelDown);
      if (!is.null(geneColumn) && hitsByGene && geneColumn %in% colnames(x)) {
         ## Aggregate hit counts by gene instead of by assay/probe/whatever
         if (verbose) {
            printDebug("Aggregating hit totals by gene column '", geneColumn, "'", c("orange", "lightblue"));
         }
         pvHitsUp1 <- x[(which(x[1:n,pvCol] <= pvalueCutoff & (x[1:n,lfcCol]) >= lfcCutoff)),geneColumn];
         multiGenesUp <- tcountPaste(pvHitsUp1[pvHitsUp1 %in% names(tcount(pvHitsUp1)[tcount(pvHitsUp1) > 1])]);
         pvHitsDown1 <- x[(which(x[1:n,pvCol] <= pvalueCutoff & -(x[1:n,lfcCol]) >= lfcCutoff)),geneColumn];
         multiGenesDown <- tcountPaste(pvHitsDown1[pvHitsDown1 %in% names(tcount(pvHitsDown1)[tcount(pvHitsDown1) > 1])]);

         if (printMultiGeneHits) {
            title(adj=0, sub=paste("Genes Multiply-Up:", wordWrap(multiGenesUp, 70, cleanText=FALSE, justify="left"), "\n",
                                   "Genes Multiply-Down:", wordWrap(multiGenesDown, 70, cleanText=FALSE, justify="left"), sep=""),
                  cex.sub=0.6, line=-2, outer=TRUE);
         }
         pvHitsUp <- length(unique(pvHitsUp1));
         pvHitsDown <- length(unique(pvHitsDown1));
         hitType <- "genes";
         labelUp1 <- labelUp;
         labelDown1 <- labelDown;
         labelUp <- paste("(", format(big.mark=",", pvHitsUp), hitType, ")");
         labelUp <- gsub("((^|[^0-9])1 .+)s$", "\\1", labelUp);
         labelUp <- paste(labelUp1, labelUp, sep="\n");
         labelDown <- paste("(", format(big.mark=",", pvHitsDown), hitType, ")");
         labelDown <- gsub("((^|[^0-9])1 .+)s$", "\\1", labelDown);
         labelDown <- paste(labelDown1, labelDown, sep="\n");
      }
      labelX <- c(range(xValues1)[1] - range(xValues1)[1]/5,
                  range(xValues1)[2] - range(xValues1)[2]/5);
      labelY <- rep(range(yValues1)[2] - range(xValues1)[2]/5, 2);
      if (!addPlot && labelHits && !doBlockArrows) {
         if (orientation %in% c("right")) {
            text(y=labelX, x=labelY, label=c(labelDown, labelUp), srt=270);
         } else {
            text(x=labelX, y=labelY, label=c(labelDown, labelUp));
         }
      }

      ## Block Arrows
      parList[["preBlockArrows"]] <- par();
      par("xpd"=FALSE);
      if (orientation %in% c("right")) {
         logAxis(1, at=unique(as.integer(pretty(yRange, n=nYlabels))), value=FALSE, base=10, makeNegative=TRUE, cex.axis=cex.axis*0.8, ...);
         logAxis(2, at=unique(as.integer(pretty(xRange, n=nXlabels))), value=TRUE, base=2, cex.axis=cex.axis*0.8, ...);
         abline(v=-log10(pvalueCutoff), h=unique(c(-lfcCutoff, lfcCutoff)), lty="dashed", col=ablineColor);
         if (doBlockArrows) {
            hitCol <- hsv(h=0.06, s=0.75, v=0.9, alpha=1);
            #if (doShadowTextArrows) {
            #   assign("text", get("shadowText"), envir=.GlobalEnv);
            #}
            # blockArrowLabelHitColor="#FFFFAA", blockArrowLabelUpColor="#FFFFAA", blockArrowLabelDownColor="#FFFFAA"
            blockArrow(axisPosition="topAxis", xleft=-log10(pvalueCutoff), arrowPosition="right",
                        blockWidthPercent=4*blockArrowCex[2], arrowLengthPercent=4*blockArrowCex[2]/2,
                        arrowLabel=paste(format(big.mark=",", pvHits), " significant ", hitType, sep=""), col=hitCol,
                        labelCex=blockArrowLabelCex, labelFont=blockArrowFont, arrowLabelColor=blockArrowLabelHitColor,
                        arrowLabelBorder=alpha2col(setTextContrastColor(blockArrowLabelHitColor), 0.3));
            blockArrow(axisPosition="rightAxis", ybottom=lfcCutoff, arrowPosition="top",
                        blockWidthPercent=4*blockArrowCex[1], arrowLengthPercent=4*blockArrowCex[1],
                        arrowLabel=labelUp, col="#990000FF", labelCex=blockArrowLabelCex,
                        labelFont=blockArrowFont, arrowLabelColor=blockArrowLabelUpColor,
                        arrowLabelBorder=alpha2col(setTextContrastColor(blockArrowLabelUpColor), 0.3));
            blockArrow(axisPosition="rightAxis", ytop=-lfcCutoff, arrowPosition="bottom",
                        blockWidthPercent=4*blockArrowCex[1], arrowLengthPercent=4*blockArrowCex[1],
                        arrowLabel=labelDown, col="#000099FF", labelCex=blockArrowLabelCex,
                        labelFont=blockArrowFont, arrowLabelColor=blockArrowLabelDownColor,
                        arrowLabelBorder=alpha2col(setTextContrastColor(blockArrowLabelDownColor), 0.3));
            if (doShadowTextArrows) {
               assign("text", graphics:::text, envir=.GlobalEnv);
            }
         }
      } else {
         logAxis(2, at=unique(as.integer(pretty(yRange, n=nYlabels))), value=FALSE, base=10, makeNegative=TRUE, cex.axis=cex.axis*0.8, ...);
         logAxis(1, at=unique(as.integer(pretty(xRange, n=nXlabels))), value=TRUE, base=2, cex.axis=cex.axis*0.8, ...);
         abline(h=-log10(pvalueCutoff), v=unique(c(-lfcCutoff, lfcCutoff)), lty="dashed", col=ablineColor);
         if (doBlockArrows) {
            hitCol <- hsv(h=0.06, s=0.75, v=0.9, alpha=1);
            #if (doShadowTextArrows) {
            #   assign("text", get("shadowText", envir=CBioRUtils), envir=.GlobalEnv);
            #}
            if (doTopHist) {
               rightAdj <- 1.2;
            } else {
               rightAdj <- 1;
            }
            blockArrow(
               axisPosition="rightAxis",
               ybottom=-log10(pvalueCutoff),
               arrowPosition="top",
               blockWidthPercent=4*blockArrowCex[2]*rightAdj,
               arrowLengthPercent=4*blockArrowCex[2]*rightAdj/2,
               arrowLabel=paste(format(big.mark=",", pvHits), " significant ", hitType, sep=""),
               col=blockArrowHitColor,
               labelCex=blockArrowLabelCex,
               labelFont=blockArrowFont,
               arrowLabelColor=blockArrowLabelHitColor,
               arrowLabelBorder=alpha2col(setTextContrastColor(blockArrowLabelHitColor), 0.3));
            blockArrow(
               axisPosition="topAxis",
               xleft=lfcCutoff,
               arrowPosition="right",
               blockWidthPercent=4*blockArrowCex[1],
               arrowLengthPercent=4*blockArrowCex[1],
               arrowLabel=labelUp,
               col=blockArrowUpColor,
               labelCex=blockArrowLabelCex,
               labelFont=blockArrowFont,
               arrowLabelColor=blockArrowLabelUpColor,
               arrowLabelBorder=alpha2col(setTextContrastColor(blockArrowLabelUpColor), 0.3));
            blockArrow(
               axisPosition="topAxis",
               xright=-lfcCutoff,
               arrowPosition="left",
               blockWidthPercent=4*blockArrowCex[1],
               arrowLengthPercent=4*blockArrowCex[1],
               arrowLabel=labelDown,
               col=blockArrowDownColor,
               labelCex=blockArrowLabelCex,
               labelFont=blockArrowFont,
               arrowLabelColor=blockArrowLabelDownColor,
               arrowLabelBorder=alpha2col(setTextContrastColor(blockArrowLabelDownColor), 0.3));
            if (doShadowTextArrows) {
               assign("text", graphics:::text, envir=.GlobalEnv);
            }
         }
      }
      parList[["postBlockArrows"]] <- par();

      ## Display the significance cutoff used
      lfcCutoffLabel <- signif(digits=2, lfcCutoff);

      if (doStatsSubtitles) {
         if (lfcCutoff > 0) {
            subTitle <- paste(sep="", "Significance cutoff ≤ ", pvalueCutoff,
               ", and ", "log2 ",
               "fold change cutoff ≥ ", lfcCutoffLabel);
         } else {
            subTitle <- paste("Significance cutoff ≤ ", pvalueCutoff);
         }
         title(sub=subTitle, cex.sub=statsCex, line=statsLine, ...);

         ## Display the total points
         totalPointsSub <- paste("Total points: ", format(big.mark=",", nrow(x)));
         parMar <- par("mar");
         par("mar"=c(parMar[1], 0.5, parMar[3:4]));
         title(sub=totalPointsSub, adj=0.01, cex.sub=statsCex, line=statsLine, ...);
         par("mar"=parMar);
      }

      ## Overall Title
      if (doTopHist && !orientation %in% c("right")) {
         parMfg(mfg=c(1,1), usePar=parList[["postTopHist"]]);
         ## For some reason, the title gets cropped unless I run it with outer=TRUE first,
         ## thereafter it isn't cropped, although no par() values change at all! Bug.
         title(main="", outer=TRUE, line=0);
      }
      par("xpd"=TRUE);
      if (!is.null(subMain)) {
         title(main=main, line=lineMain, font.main=font.main, cex.main=cexMain);
         title(main=subMain, line=lineMain-1, cex.main=cexSub, font.main=font.main);
      } else {
         title(main=main, line=lineMain, font.main=font.main, cex.main=cexMain);
      }
      par("xpd"=FALSE);
   }

   if (doTopHist) {
      par("mar"=origPar$mar);
      par("plt"=origPar$plt);
      par("usr"=origPar$usr);
   }
   if (length(pointSubset) > 0) {
      vp1DF <- data.frame(check.names=FALSE,
         stringsAsFactors=FALSE,
         "pointSubset"=pointSubset,
         "pointSubsetX"=xValues1[pointSubset],
         "pointSubsetY"=yValues1[pointSubset],
         "pointSubsetCex"=pointSubsetCex,
         "pointSubsetCol"=pvColors1[pointSubset]);
      rownames(vp1DF) <- names(pointSubset);
      invisible(list(
         x=xValues1,
         y=yValues1,
         col=pvColors1,
         pointSubsetDF=vp1DF,
         parList=parList,
         labelCoords=labelCoords));
   } else if (length(multiGenesUp) > 0 || length(multiGenesDown) > 0) {
      invisible(list(
         x=xValues1,
         y=yValues1,
         col=pvColors1,
         multiGenesUp=multiGenesUp,
         multiGenesDown=multiGenesDown,
         parList=parList,
         labelCoords=labelCoords));
   } else {
      invisible(list(
         x=xValues1,
         y=yValues1,
         col=pvColors1,
         parList=parList,
         labelCoords=labelCoords));
   }
}

blockArrow <- function
(axisPosition="rightAxis",
 arrowPosition="top",
 arrowLabel="",
 arrowDirection="updown",
 xleft=NULL, xright=NULL, ybottom=NULL, ytop=NULL,
 col="#660000FF", labelFont=1, labelCex=1,
 border="#000000FF", xpd=TRUE, parUsr=par("usr"),
 arrowWidthPercent=10, arrowLengthPercent=4,
 blockWidthPercent=4,
 blankFirst=FALSE,
 blankColor="white",
 arrowLabelColor="#FFFFFFFF",
 arrowLabelBorder="#000000FF",
 doExample=FALSE,
 doBlockGradient=TRUE,
 gradientDarkFactor=1.5,
 gradientSFactor=1.5,
 verbose=FALSE,
 ...)
{
   ## Purpose is to draw a block arrow, like ones you see in PowerPoint, but using
   ## rect() syntax as if drawing a rectangle.
   ##
   ## Currently the function draws block arrows outside the plot, as if to label
   ## an axis.
   ##
   ## Trim some of the width away to give room for the arrow to be drawn.
   ##
   ## Some examples:
   if (doExample) {
      nullPlot();
      blockArrow(axisPosition="rightAxis",
         ybottom=1.52,
         arrowPosition="top",
         col="#660000FF",
         arrowLabel="Up-regulated",
         arrowWidthPercent=arrowWidthPercent,
         blockWidthPercent=blockWidthPercent,
         labelCex=labelCex,
         ...);
      blockArrow(axisPosition="rightAxis",
         ytop=1.48,
         arrowPosition="bottom",
         col="#000066FF",
         arrowLabel="Down-regulated",
         arrowWidthPercent=arrowWidthPercent,
         blockWidthPercent=blockWidthPercent,
         labelCex=labelCex,
         ...);
      blockArrow(axisPosition="topAxis",
         xright=1.48,
         arrowPosition="left",
         col="#660000FF",
         arrowLabel="Down-regulated",
         arrowWidthPercent=arrowWidthPercent,
         blockWidthPercent=blockWidthPercent,
         labelCex=labelCex,
         ...);
      blockArrow(axisPosition="topAxis",
         xleft=1.52,
         arrowPosition="right",
         col="#000066FF",
         arrowLabel="Up-regulated",
         arrowWidthPercent=arrowWidthPercent,
         blockWidthPercent=blockWidthPercent,
         labelCex=labelCex,
         ...);
      return(NULL);
   }
   if (length(axisPosition) > 1) {
      retVals <- lapply(axisPosition, function(iAxisPosition) {
         blockArrow(axisPosition=iAxisPosition,
            arrowPosition=arrowPosition,
            arrowLabel=arrowLabel,
            arrowDirection=arrowDirection,
            xleft=xleft,
            xright=xright,
            ybottom=ybottom,
            ytop=ytop,
            col=col,
            labelFont=1,
            labelCex=1,
            border=border,
            xpd=xpd,
            parUsr=parUsr,
            arrowWidthPercent=arrowWidthPercent,
            arrowLengthPercent=arrowLengthPercent,
            blockWidthPercent=blockWidthPercent,
            blankFirst=blankFirst,
            arrowLabelColor=arrowLabelColor,
            arrowLabelBorder=arrowLabelBorder,
            doExample=doExample,
            doBlockGradient=doBlockGradient,
            gradientDarkFactor=gradientDarkFactor,
            gradientSFactor=gradientSFactor,
            ...);
      });
      invisible(retVals);
   }

   if (doBlockGradient) {
      if (!suppressPackageStartupMessages(require(plotrix))) {
         if (verbose) {
            printDebug("Cannot use color gradient, it requires the 'plotrix' package.",
               fgText="red");
         }
         doBlockGradient <- FALSE;
      } else {
         col1 <- col[1];
         if (length(col) <= 2) {
            col2 <- makeColorDarker(col[length(col)],
               darkFactor=gradientDarkFactor,
               sFactor=gradientSFactor);
            colGradient <- colorRampPalette(c(col, col2),
               alpha=TRUE)(25);
         } else {
            colGradient <- colorRampPalette(col, alpha=TRUE)(25);
            col2 <- col[length(col)];
         }
      }
   }
   #usrBox(fill="#FFFF9933", label="before");

   arrowSets <- list("empty"=list(x=NULL, y=NULL));
   if (igrepHas("rightAxis|leftAxis", axisPosition)) {
      if (is.null(ybottom)) {
         ybottom <- parUsr[3];
      }
      if (is.null(ytop)) {
         ytop <- parUsr[4];
      }
      if (igrepHas("rightAxis", axisPosition)) {
         if (is.null(xleft)) {
            xleft <- parUsr[2];
         }
         if (is.null(xright)) {
            xright <- parUsr[2] +
               (parUsr[2] - parUsr[1]) * blockWidthPercent/100;
         }
      } else {
         if (is.null(xright)) {
            xright <- parUsr[1];
         }
         if (is.null(xleft)) {
            xleft <- parUsr[1] -
               (parUsr[2] - parUsr[1]) * blockWidthPercent/100;
         }
      }
      arrowDirection <- "updown";
      srtLabel <- 90;
      gradientXY <- "y";
   } else if (igrepHas("topAxis|bottomAxis", axisPosition)) {
      if (igrepHas("topAxis", axisPosition)) {
         if (is.null(ybottom)) {
            ybottom <- parUsr[4];
         }
         if (is.null(ytop)) {
            ytop <- parUsr[4] +
               (parUsr[4] - parUsr[3]) * blockWidthPercent/100;
         }
      } else {
         if (is.null(ytop)) {
            ytop <- parUsr[3];
         }
         if (is.null(ybottom)) {
            ybottom <- parUsr[3] -
               (parUsr[4] - parUsr[3]) * blockWidthPercent/100;
         }
      }
      if (is.null(xleft)) {
         xleft <- parUsr[1];
      }
      if (is.null(xright)) {
         xright <- parUsr[2];
      }
      arrowDirection <- "leftright";
      srtLabel <- 0;
      gradientXY <- "x";
   }
   if (blankFirst) {
      #printDebug("blanking the area first");
      polygon(
         x=c(xleft, xright, xright, xleft),
         y=c(ybottom, ybottom, ytop, ytop),
         col=blankColor,
         border=blankColor,
         xpd=TRUE);
   }

   if (length(igrep("up|down", arrowDirection)) > 0) {
      ## Trim away the y-coordinates
      arrowWidth <- abs(xright - xleft);
      arrowWidthDiff <- arrowWidth * (arrowWidthPercent / 100);
      if (xright > xleft) {
         xright1 <- xright - arrowWidthDiff;
         xleft1 <- xleft + arrowWidthDiff;
      } else {
         xright1 <- xright + arrowWidthDiff;
         xleft1 <- xleft - arrowWidthDiff;
      }
      #printDebug(c("xright1: ", xright1));
      #printDebug(c("xleft1: ", xleft1));
      arrowLength <- abs(ytop - ybottom);
      arrowLengthDiff <- arrowLength * (arrowLengthPercent / 100);
      if (ytop > ybottom) {
         if (length(igrep("top", arrowPosition)) > 0) {
            ytop1 <- ytop - arrowLengthDiff;
            yArrowPoints <- c(ytop1, ytop1, ytop, ytop1, ytop1);
            xArrowPoints <- c(xleft1, xleft, mean(c(xleft, xright)), xright, xright1);
            arrowSets <- c(arrowSets, list("top"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop1, ybottom));
            xLabel <- mean(c(xleft, xright));
         } else {
            yArrowPoints <- c(ytop, ytop);
            xArrowPoints <- c(xleft1, xright1);
            arrowSets <- c(arrowSets, list("top"=list(x=xArrowPoints, y=yArrowPoints)));
            #print(arrowSets);
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright));
         }
         if (length(igrep("bottom", arrowPosition)) > 0) {
            ybottom1 <- ybottom + arrowLengthDiff;
            yArrowPoints <- c(ybottom1, ybottom1, ybottom, ybottom1, ybottom1);
            xArrowPoints <- rev(c(xleft1, xleft, mean(c(xleft, xright)), xright, xright1));
            arrowSets <- c(arrowSets, list("bottom"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom1));
            xLabel <- mean(c(xleft, xright));
         } else {
            yArrowPoints <- c(ybottom, ybottom);
            xArrowPoints <- c(xright1, xleft1);
            arrowSets <- c(arrowSets, list("bottom"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright));
         }
      } else {
         if (length(igrep("bottom", arrowPosition)) > 0) {
            ytop1 <- ytop + arrowLengthDiff;
            yArrowPoints <- c(ytop1, ytop1, ytop, ytop1, ytop1);
            xArrowPoints <- rev(c(xleft1, xleft, mean(c(xleft, xright)), xright, xright1));
            arrowSets <- c(arrowSets, list("bottom"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop1, ybottom));
            xLabel <- mean(c(xleft, xright));
         } else {
            yArrowPoints <- c(ytop, ytop);
            xArrowPoints <- c(xright1, xleft1);
            arrowSets <- c(arrowSets, list("bottom"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright));
         }
         if (length(igrep("top", arrowPosition)) > 0) {
            ybottom1 <- ybottom - arrowLengthDiff;
            yArrowPoints <- c(ybottom1, ybottom1, ybottom, ybottom1, ybottom1);
            xArrowPoints <- c(xleft1, xleft, mean(c(xleft, xright)), xright, xright1);
            arrowSets <- c(arrowSets, list("top"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom1));
            xLabel <- mean(c(xleft, xright));
         } else {
            yArrowPoints <- c(ybottom, ybottom);
            xArrowPoints <- c(xleft1, xright1);
            arrowSets <- c(arrowSets, list("top"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright));
         }
      }
   } else {
      ## Arrow position is left-to-right
      ## Trim away the x-coordinates
      arrowWidth <- abs(ytop - ybottom);
      arrowWidthDiff <- arrowWidth * (arrowWidthPercent / 100);
      if (ytop > ybottom) {
         ytop1 <- ytop - arrowWidthDiff;
         ybottom1 <- ybottom + arrowWidthDiff;
      } else {
         ytop1 <- ytop + arrowWidthDiff;
         ybottom1 <- ybottom - arrowWidthDiff;
      }
      #printDebug(c("ytop1: ", ytop1));
      #printDebug(c("ybottom1: ", ybottom1));
      arrowLength <- abs(xright - xleft);
      arrowLengthDiff <- arrowLength * (arrowLengthPercent / 100);

      if (xright > xleft) {
         if (length(igrep("right", arrowPosition)) > 0) {
            xright1 <- xright - arrowLengthDiff;
            xArrowPoints <- c(xright1, xright1, xright, xright1, xright1);
            yArrowPoints <- c(ybottom1, ybottom, mean(c(ybottom, ytop)), ytop, ytop1);
            arrowSets <- c(arrowSets, list("right"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright1));
         } else {
            xright1 <- xright;
            xArrowPoints <- c(xright, xright);
            yArrowPoints <- c(ybottom1, ytop1);
            arrowSets <- c(arrowSets, list("right"=list(x=xArrowPoints, y=yArrowPoints)));
            #print(arrowSets);
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft, xright));
         }
         if (length(igrep("left", arrowPosition)) > 0) {
            xleft1 <- xleft + arrowLengthDiff;
            xArrowPoints <- c(xleft1, xleft1, xleft, xleft1, xleft1);
            yArrowPoints <- rev(c(ybottom1, ybottom, mean(c(ybottom, ytop)), ytop, ytop1));
            arrowSets <- c(arrowSets, list("left"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft1, xright1));
         } else {
            xleft1 <- xleft;
            xArrowPoints <- c(xleft, xleft);
            yArrowPoints <- rev(c(ybottom1, ytop1));
            arrowSets <- c(arrowSets, list("bottom"=list(x=xArrowPoints, y=yArrowPoints)));
            yLabel <- mean(c(ytop, ybottom));
            xLabel <- mean(c(xleft1, xright1));
         }
      } else {
      }

   }
   yLabel <- mean(c(ytop, ybottom));
   xLabel <- mean(c(xleft, xright));

   allX <- unlist(lapply(unvigrep("^empty$", names(arrowSets)), function(i){
      j <- arrowSets[[i]];
      j["x"];
   }))
   allY <- unlist(lapply(unvigrep("^empty$", names(arrowSets)), function(i){
      j <- arrowSets[[i]];
      j["y"];
   }));

   if (doBlockGradient) {
      arrowSides <- sapply(unvigrep("^empty$", names(arrowSets)), function(i){
         length(arrowSets[[i]][["x"]]) > 2;
      });
      arrowSides <- paste(names(arrowSides)[arrowSides], collapse="");
      if (arrowSides %in% c("bottom", "left")) {
         colGradient <- rev(colGradient);
         col21 <- col1;
         col1 <- col2;
         col2 <- col1;
      }
      arrowBox <- lapply(nameVector(unvigrep("^empty$", names(arrowSets))), function(i){
         xi <- arrowSets[[i]][["x"]];
         yi <- arrowSets[[i]][["y"]];
         xI <- unique(c(xi[1], xi[length(xi)]));
         yI <- unique(c(yi[1], yi[length(yi)]));
         list(x=xI, y=yI);#, xi=xi, yi=yi);
      });
      arrowBoxX <- unique(sort(c(arrowBox[[1]][["x"]], arrowBox[[2]][["x"]])));
      arrowBoxY <- unique(sort(c(arrowBox[[1]][["y"]], arrowBox[[2]][["y"]])));

      xpdPar <- par("xpd");
      par("xpd"=TRUE);
      polygon(
         x=allX,
         y=allY,
         col=tail(colGradient,1),
         border=border,
         xpd=TRUE);
      gr1 <- gradient.rect(
         col=colGradient,
         gradient=gradientXY,
         xleft=arrowBoxX[1],
         xright=arrowBoxX[2],
         ybottom=arrowBoxY[1],
         ytop=arrowBoxY[2],
         border="#00000000");
      par("xpd"=xpdPar);
      arrowsDrawn1 <- lapply(nameVector(vgrep("top|right", names(arrowSets))), function(i){
         as1 <- arrowSets[[i]];
         polygon(
            x=as1[["x"]],
            y=as1[["y"]],
            col=col2,
            border=NA,
            xpd=TRUE);
         as1;
      });
      arrowsDrawn2 <- lapply(nameVector(vgrep("bottom|left", names(arrowSets))), function(i){
         as1 <- arrowSets[[i]];
         polygon(
            x=as1[["x"]],
            y=as1[["y"]],
            col=col1,
            border=col1,
            xpd=TRUE);
         as1;
      });
      polygon(
         x=allX,
         y=allY,
         col=NA,
         border=border,
         xpd=TRUE);
   } else {
      polygon(
         x=allX,
         y=allY,
         col=col,
         border=border,
         xpd=TRUE);
      arrowsDrawn1 <- NULL;
      arrowsDrawn2 <- NULL;
   }

   #usrBox(fill="#FFFF9933", label="after");

   ## Optionally label the block arrows
   if (!is.null(arrowLabel) && !arrowLabel %in% c(NA, "")) {
      if (verbose) {
         printDebug("blockArrow(): ",
            "xLabel: ",
            round(digits=2, xLabel),
            ",  yLabel: ",
            round(digits=2, yLabel));
      }
      text(
         label=arrowLabel,
         x=xLabel,
         y=yLabel,
         srt=srtLabel,
         xpd=TRUE,
         col=arrowLabelColor,
         font=labelFont,
         cex=labelCex,
         bg=arrowLabelBorder,
         adj=c(0.5,0.5),
         ...);
      #mtext(text=arrowLabel, x=xLabel, y=yLabel, srt=srtLabel);
   }
   retVals <- list(
      arrowSets=arrowSets,
      arrowsDrawn1=arrowsDrawn1,
      arrowsDrawn2=arrowsDrawn2,
      allX=allX,
      allY=allY)
   invisible(retVals);
}

updateFunctionParamList <- function
(functionName=NULL,
 paramName=NULL,
 newValues=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to facilitate updating default parameters which are present in a list,
   ## but where defaults are defined in a function name.  So if someone wants to change
   ## one of the default values, but keep the rest of the list of defaults, this function
   ## does it.
   ##
   ## The default values are taken from the function formals, using eval(formals(functionName)).
   defaultList <- eval(formals(functionName)[[paramName]]);
   defaultList <- updateListElements(defaultList, newValues, verbose=verbose);
   return(defaultList);
}

updateListElements <- function
(sourceList,
 updateList,
 listLayerNum=1,
 verbose=TRUE,
 ...)
{
   ## Purpose is to update elements in a list, allowing for multi-layered lists
   ## of lists. In case of a list-of-list, it will call this function again with
   ## each successive layer of listedness.
   ##
   ## Handy for updating lattice graphics settings, which are impossibly nested
   ## tangled ball of textual yarn.  An example:
   ## tp2 <- updateListElements(trellis.par.get(), list(fontsize=list(points=6)));
   ## trellis.par.set(tp2);
   ##
   ## Or in one line:
   ## trellis.par.set(updateListElements(trellis.par.get(), list(fontsize=list(points=6))));
   if (class(updateList) %in% c("list") && class(updateList[[1]]) %in% c("list")) {
      for (updateListName in names(updateList)) {
         if (verbose) {
#            printDebug(paste(rep("   ", listLayerNum-1), sep=""), "listLayerNum: ", listLayerNum, c("lightblue", "orange"));
#            printDebug(paste(rep("   ", listLayerNum-1), sep=""), "updateListName: ", updateListName, c("lightblue", "orange"));
         }
         if (updateListName %in% names(sourceList)) {
            ## If the name already exists, we must update items within the list
            sourceList[[updateListName]] <- updateListElements(
               sourceList=sourceList[[updateListName]],
               updateList=updateList[[updateListName]],
               listLayerNum=listLayerNum+1);
         } else {
            ## If the name does not already exist, we can simply add it.
            sourceList[[updateListName]] <- updateList[[updateListName]];
         }
      }
   } else {
      if (!is.null(names(updateList))) {
         if (verbose) {
#            printDebug(paste(rep("   ", listLayerNum-1), sep=""), "Updating sourceList[", names(updateList), "]: ", updateList, c("lightblue", "orange"));
         }
         sourceList[names(updateList)] <- updateList;
      } else {
         if (verbose) {
#            printDebug(paste(rep("   ", listLayerNum-1), sep=""), "Updating sourceList: ", updateList, c("lightblue", "orange"));
         }
         sourceList <- updateList;
      }
   }
   return(sourceList);
}

jitter2 <- function
(x,
 factor=1,
 amount=NULL,
 verbose=FALSE,
 ...)
{
   ## Purpose is to replicate the jitter() function, except removing duplicated
   ## points beforehand, which otherwise cause problems when calculating z.
   if (length(x) == 0L) {
      return(x);
   }
   if (!is.numeric(x)) {
      stop("'x' must be numeric");
   }
   z <- diff(r <- range(x[is.finite(x)]));
   if (z == 0) {
      z <- abs(r[1L]);
   }
   if (z == 0) {
      z <- 1;
   }
   if (verbose) {
      printDebug("jitter2() output:", c("yellow"));
      printDebug("================", c("yellow"));
      printDebug("        z: ", z, c("orange", "lightblue"));
      printDebug("   amount: ", amount, c("orange", "lightblue"));
      printDebug("   factor: ", factor, c("orange", "lightblue"));
   }
   if (is.null(amount)) {
      ## Here is the subtle change:
      zFloorOld <- floor(log10(z));
      zFloor <- floor(round(log10(z)));
      dOld <- diff(xx <- unique(sort.int(       round(x, 3 - zFloorOld))));
      d    <- diff(xx <- unique(sort.int(unique(round(x, 3 - zFloor)))));
      if (verbose) {
         printDebug("zFloorOld: ", zFloorOld, c("orange", "lightgreen"));
         printDebug("   zFloor: ", zFloor, c("orange", "lightgreen"));
         printDebug("        d: ", head(d), c("orange", "lightgreen"));
      }
      if (length(d)) {
         d <- min(d);
      } else if (xx != 0) {
         d <- xx/10;
      } else {
         d <- z/10;
      }
      if (verbose) {
         printDebug("     dOld: ", head(dOld), c("orange", "lightgreen"));
         printDebug("        d: ", d, c("orange", "lightgreen"));
      }
      amount <- factor/5 * abs(d);
      if (verbose) {
         printDebug("   amount: ", amount, c("orange", "lightgreen"));
      }
   } else if (amount == 0) {
      amount <- factor * (z/50);
      if (verbose) {
         printDebug("   amount: ", amount, c("orange", "lightblue"));
      }
   }
   x + stats::runif(length(x), -amount, amount)
}


#' Draw block arrows onto a figure
#'
#' Draw block arrows onto a figure
#'
#' This function draws block arrows onto an existing figure, such that
#' the arrow is scaled proportionally relative to plot coordinates. It
#' by default fills the block arrow with a color gradient, and optionally
#' label text which is drawn with \code{\link{shadowText}}.
#'
#' To see an example, run \code{doTest=TRUE}.
#'
#' If \code{location} is supplied, with 'top', 'bottom', 'left', or 'right'
#' then the arrows are drawn outside plot space in the axis space adjacent
#' to the plot. This function is used to label volcano plot gene hits, for
#' example.
#'
#' The shape of the arrow can be configured by adjusting
#' \code{arrowLengthPercent}, \code{arrowWidthPercent}, and
#' \code{blockWidthPercent}, which collectively determine the width
#' of the block arrow stem, and width of the arrow head at the end.
#'
#' Arrow labels are drawn so text is upright, regardless the angle of the
#' arrow, unless \code{flipLabelSrt=FALSE}. Labels by default use
#' \code{\link{shadowText}} which mostly just draws a contrasting border
#' around the text, hopefully making it visually easier to read.
#'
#' The color gradient uses \code{\link{color2gradient}} and notably uses
#' nGradientSteps to define how many subsections are used for the gradient.
#' If multiple colors are supplied, then \code{\link{colorRampPalette} is
#' used to expand those colors to the number from \code{nGradientSteps}.
#'
#'
#' @param fromX,toX,fromY,toY coordinates defining the start and end
#'    coordinates for the blockArrow. Currently these values are not
#'    vectorized, which means only one block arrow is drawn per function
#'    call.
#' @param location NULL or character string indicating whether the block
#'    arrow should be placed outside plot space in the axis margins.
#'    Allowed values are 'top', 'bottom', 'left', and 'right'. When one of
#'    those substrings is detected, the corresponding fromX,toX, or fromY,toY
#'    are fixed to the appropriate margin coordinates.
#' @param col vector of one or more R colors, used to color the block arrow.
#'    If one color is supplied, and \code{doBlockGradient=TRUE} then it is
#'    expanded into a color gradient using \code{\link{color2gradient}}.
#'    If multiple values are supplied, then \code{\link{colorRampPalette}} is
#'    used to expand colors to \code{nGradientSteps}.
#' @param border single R color used as a border around the block arrow.
#'    Supply "transparent" to suppress the border.
#' @param arrowLabel character string with a label to print inside the block
#'    arrow. Note: text is not scaled to fit the block arrow, it is simply
#'    placed at the center.
#' @param arrowLabelCex,arrowLabelFont,arrowLabelColor,arrowLabelSrt parameters
#'    controlling the label text.
#' @param flipLabelSrt logical whether to ensure the arrow label is upright,
#'    regardless the angle of the arrow. Otherwise text will flow from the stem
#'    to the arrow head, which may cause text to appear upside down.
#' @param doShadowTextArrows logical whether to render the arrow label using
#'    \code{\link{shadowText}} which provides an outline to each letter.
#' @param arrowLengthPercent,arrowWidthPercent,blockWidthPercent parameters
#'    used to define the shape of the arrow. The parameter
#'    \code{arrowSizeRelativeTo} determines how these values are calculated.
#' @param arrowSizeRelativeTo character value, either 'plot' or 'length', which
#'    helps define relative scaling of the block arrow. The 'plot' scaling
#'    ensures that multiple block arrows share common reference, being the
#'    size of the plot itself. However 'length' will allow custom sizing of
#'    the arrow head for short arrows, for example.
#' @param doBlockGradient logical whether to create a color gradient to
#'    fill the block arrow.
#' @param nGradientSteps integer number of colors used when
#'    \code{doBlockGradient=TRUE}.
#' @param gradientWtFactor when supplying one color for \code{col}, and
#'    \code{doBlockGradient=TRUE}, the \code{\link{color2gradient}} function
#'    is used to create a color ramp. The gradient width can be controlled with
#'    gradientWtFactor, where a value 1 or higher is a dramatic gradient, but
#'    values below 1 are more subtle.
#' @param lwd line width used for block arrow border
#' @param parUsr,parPin values derived from \code{par("usr")} and
#'    \code{par("pin")}, supplied here in the rare case that a custom
#'    range should be supplied.
#' @param xpd,add parameters controlling plot features. The \code{par("xpd")}
#'    value controls clipping block arrows when they are outside the plot area,
#'    and \code{add} controls whether arrows are added to an existing plot,
#'    are used to create a new plot.
#' @param verbose logical whether to print verbose output.
#' @param doTest logical whether to perform a short test to demonstrate
#'    the capabilities of this function.
#'
#' @examples
#' drawBlockArrow(doTest=TRUE);
#' drawBlockArrow(fromX=20, toX=80, fromY=1500, toY=1500,
#'    col="blue4", arrowLabel="Blue block arrow")
#'
#' @export
drawBlockArrow <- function
(fromX=NULL,
 toX=NULL,
 fromY=NULL,
 toY=NULL,
 location=NULL,
 col="#660000FF",
 border="#000000FF",
 arrowLabel="",
 arrowLabelCex=1,
 arrowLabelFont=1,
 arrowLabelColor="white",
 arrowLabelSrt=0,
 flipLabelSrt=TRUE,
 arrowLengthPercent=4,
 arrowWidthPercent=10,
 doShadowTextArrows=TRUE,
 blockWidthPercent=4,
 arrowSizeRelativeTo=c("plot", "length"),
 doBlockGradient=TRUE,
 nGradientSteps=20,
 gradientWtFactor=1/2,
 lwd=1,
 parUsr=NULL,
 parPin=NULL,
 xpd=TRUE,
 add=TRUE,
 verbose=FALSE,
 doTest=FALSE,
 ...)
{
   ## An extension of the "blockArrow()" function, which only
   ## currently draws axis block arrows. This function allows
   ## any arbitrary arrow to be drawn, shaded with gradient
   ## colors, then rotated within the plot space
   ##
   ## The hardest part may be making sure the proper aspect ratio
   ## is maintained, when the x- and y-axes are on different scales
   ##
   ## flipLabelSrt=TRUE will try to keep labels upright when
   ## they may otherwise appear upside down.
   ##
   ## TODO: enable arcSegments() type logic, in order to draw curved
   ## block arrows. It would require drawing arced rectangles, so the
   ## width of the rectangle is maintained along the right angle to the
   ## path of the arrow. See addArrows3d()
   ##
   ## Save graphical parameters, since we may need to adjust line
   ## parameters for lmitre, lend, ljoin
   if (!suppressPackageStartupMessages(require(sp))) {
      stop("drawBlockArrow() requires the sp package for Polygons objects.");
   }
   if (!suppressPackageStartupMessages(require(maptools))) {
      stop("drawBlockArrow() requires the maptools package for the elide() function.");
   }


   ## test example:
   if (doTest) {
      nullPlot(xlim=c(1,100), ylim=c(1,2000));
      par("las"=2);
      axis(1); axis(2);
      fromX <- 20;
      toX  <- 60;
      fromY <- 100;
      toY <- 1600;
      db1 <- drawBlockArrow(fromX=fromX, toX=toX,
         fromY=fromY, toY=toY,
         arrowLabel="Example Block Arrow",
         col="red4");
      db2 <- drawBlockArrow(fromX=60, toX=100,
         location="top", arrowLabel="location='top'",
         col="green4");
      db3 <- drawBlockArrow(fromX=40, toX=90,
         fromY=1500, toY=500,
         arrowLabel="Another arrow",
         col="blue4");
      return(NULL);
   }

   oPar <- par();
   par("lend"="butt");    # 1=butt, 2=square
   par("ljoin"="mitre");   # 0=round, 1=mitre, 2=bevel
   par("lmitre"=50); # num lines when mitre is converted to bevel
   par("xpd"=xpd);

   ## par("usr") is of the form c(x1, x2, y1, y2) with plot dimensions
   ## par("pin") is of the form c(x, y) with plot dimensions in inches
   ## getPlotLayout() literally combines the effects of the plotted x range, plotted y range, and visual x/y aspect ratio
   if (is.null(parUsr)) {
      parUsr <- par("usr");
   }
   if (is.null(parPin)) {
      parPin <- par("pin");
   }
   plotAspect <- getPlotAspect(parUsr=parUsr, parPin=parPin);

   if (length(location) > 0) {
      diffY <- diff(parUsr[3:4])/34;
      diffX <- diff(parUsr[1:2])/34;
      if (igrepHas("top", location)) {
         ## Top axis
         fromY <- parUsr[4] + diffY;
         toY <- parUsr[4] + diffY;
      } else if (igrepHas("right", location)) {
         ## Right axis
         fromX <- parUsr[2] + diffX;
         toX <- parUsr[2] + diffX;
      } else  if (igrepHas("bottom", location)) {
         ## Bottom axis
         fromY <- parUsr[3] - diffY;
         toY <- parUsr[3] - diffY;
      } else if (igrepHas("left", location)) {
         ## Left axis
         fromX <- parUsr[1] - diffX;
         toX <- parUsr[1] - diffX;
      }
   }

   if (verbose) {
      printDebug("parUsr: ", parUsr);
      printDebug("parPin: ", paste(collapse=" x ",
         paste0(format(trim=TRUE, digits=2, parPin), '"')));
      printDebug("plotAspect: ", format(digits=2, plotAspect));
   }
   plotWidth <- diff(parUsr[1:2]);
   plotHeight <- diff(parUsr[3:4]);

   ## Define dimensions of the plot using the visible x-axis,
   ## then use plot aspect to scale after we rotate the polygons
   arrowSizeRelativeTo <- match.arg(arrowSizeRelativeTo);
   arrowWidth <- (toX - fromX);
   arrowHeight <- (toY - fromY);
   arrowLength <- sqrt(arrowWidth^2 + arrowHeight^2);
   arrowLengthScaled <- sqrt((toX - fromX)^2 + ((toY - fromY)/plotAspect)^2);
   arrowAngleScaled <- rad2deg(atan2(arrowHeight/plotAspect, arrowWidth));
   arrowLabelSrtScaled <- ifelse(arrowAngleScaled > 90,
      arrowAngleScaled-180,
      arrowAngleScaled+0);

   arrowLength <- arrowLengthScaled;
   if (arrowSizeRelativeTo %in% c("plot")) {
      arrowLengthDiff <- plotWidth * (arrowLengthPercent / 100);
   } else {
      arrowLengthDiff <- arrowLength * (arrowLengthPercent / 100);
   }
   arrowWidth <- plotWidth * blockWidthPercent/100;
   arrowWidthDiff <- arrowWidth * (arrowWidthPercent / 100);

   ## Define a block arrow body with "unit coordinates",
   ## that we will rotate and scale into proper user coordinates.
   ## First the wide stem of the arrow
   stemXleft <- 0;
   stemXright <- arrowLength-arrowLengthDiff;
   stemYbottom <- -arrowWidth/2;
   stemYtop <- arrowWidth/2;
   if (doBlockGradient && nGradientSteps > 1) {
      gradientX <- seq(from=stemXleft, to=stemXright, length.out=nGradientSteps);
      stemPolySet <- lapply(seq_along(gradientX)[-1], function(i) {
         xStep1 <- gradientX[i-1];
         xStep2 <- gradientX[i];
         polyStepX <- c(xStep1, xStep2, xStep2, xStep1, xStep1);
         polyStepY <- c(stemYtop, stemYtop, stemYbottom, stemYbottom, stemYtop);
         stemPoly <- Polygon(do.call(cbind, list(x=c(polyStepX), y=c(polyStepY))));
         stemPolys <- Polygons(list(stemPoly), ID=paste("arrowStem", i, sep="_"));
      });
      #stemSP <- SpatialPolygons(stemPolySet);
      if (length(col) == 1) {
         colSet <- color2gradient(col, nGradientSteps,
            gradientWtFactor=gradientWtFactor);
      } else {
         colSet <- colorRampPalette(col, space="Lab")(nGradientSteps);
         colsetAlpha <- alpha2hex(approx(col2alpha(col), n=nGradientSteps)$y);
         colSet <- paste(colSet, colsetAlpha, sep="");
         printDebug("colSet: ", colSet, c("orange", "lightblue"));
      }
      wholeArrowColor <- "transparent";
   } else {
      polyStepX <- c(stemXleft, stemXright, stemXright, stemXleft, stemXleft);
      polyStepY <- c(stemYtop, stemYtop, stemYbottom, stemYbottom, stemYtop);
      stemPoly <- Polygon(do.call(cbind, list(x=c(polyStepX), y=c(polyStepY))));
      stemPolys <- Polygons(list(stemPoly), ID="arrowStem");
      #stemSP <- SpatialPolygons(list(stemPolys));
      stemPolySet <- list(stemPolys);
      nGradientSteps <- 2;
      wholeArrowColor <- col[1];
      colSet <- rep("transparent", length.out=2);
   }

   ## Draw arrow but include points from the stem so we allow
   ## line joining to happen
   arrowXset <- c(stemXright, arrowLength, stemXright,
      stemXright, stemXright, stemXright);
   arrowYset <- c(stemYtop+arrowWidthDiff, 0, stemYbottom-arrowWidthDiff,
      stemYbottom, stemYtop, stemYtop+arrowWidthDiff);
   arrowPoly <- Polygon(do.call(cbind, list(x=c(arrowXset), y=c(arrowYset))));
   arrowPolys <- Polygons(list(arrowPoly), ID="arrowHead");

   ## Create a polygon representing the border around the whole block arrow
   wholeArrowXset <- c(stemXleft, stemXright, stemXright,
      arrowLength, stemXright, stemXright, stemXleft, stemXleft);
   wholeArrowYset <- c(stemYtop, stemYtop, stemYtop+arrowWidthDiff,
      0, stemYbottom-arrowWidthDiff, stemYbottom, stemYbottom, stemYtop);
   wholeArrowPoly <- Polygon(do.call(cbind,
      list(x=c(wholeArrowXset), y=c(wholeArrowYset))));
   wholeArrowPolys <- Polygons(list(wholeArrowPoly), ID="wholeArrow");

   ## Put together all the pieces
   blockArrowSP <- SpatialPolygons(c(stemPolySet, arrowPolys, wholeArrowPolys),
      pO=1:(nGradientSteps+1));
   blockArrowSPsub <- SpatialPolygons(c(stemPolySet, arrowPolys),
      pO=1:(nGradientSteps));
   ## We need to draw borders around each rectangle, otherwise there
   ## are tiny gaps between them where the rendering is not precise.
   borderSet <- c(colSet, border);
   colorSet <- c(colSet, wholeArrowColor);

   ## transform the arrow into user coordinates
   blockArrowSP2 <- elide(blockArrowSP,
      rotate=-arrowAngleScaled,
      center=c(0,0));
   blockArrowSP2@plotOrder <- 1:(nGradientSteps+1);
   ## Scale the polygons
   blockArrowSP3 <- resizeSpatialPolygons(blockArrowSP2,
      size=c(1, plotAspect),
      offsetXY=c(fromX, fromY));
   blockArrowSP3bbox <- bbox(blockArrowSP3);
   arrowXlabel <- mean(blockArrowSP3bbox["x",]);
   arrowYlabel <- mean(blockArrowSP3bbox["y",]);

   ## Plot using sp::plot generic function
   plot(blockArrowSP3, col=colorSet, border=borderSet, lwd=lwd, add=add);

   ## Optionally print a label inside the block arrow
   if (!arrowLabel %in% c("", NA)) {
      if (verbose) {
         printDebug("arrowAngleScaled (before):", arrowAngleScaled);
      }
      if (flipLabelSrt) {
         arrowAngleScaled <- trimAngle(arrowAngleScaled, minAngle=0, maxAngle=180);
         if (arrowAngleScaled > 90) {
            arrowAngleScaled <- 180 + arrowAngleScaled;
         }
      }
      if (verbose) {
         printDebug("arrowAngleScaled  (after):", arrowAngleScaled);
         printDebug("arrowLabelSrtScaled:", arrowLabelSrtScaled);
      }
      #if (doShadowTextArrows) {
      #   shadowText(x=arrowXlabel, y=arrowYlabel, label=arrowLabel,
      #      cex=arrowLabelCex, adj=0.5, srt=arrowAngleScaled,
      #      col=arrowLabelColor, font=arrowLabelFont);
      #} else {
         text(x=arrowXlabel, y=arrowYlabel, label=arrowLabel,
            cex=arrowLabelCex, adj=0.5, srt=arrowAngleScaled,
            col=arrowLabelColor, font=arrowLabelFont);
      #}
   }
   par(oPar);
   invisible(list(blockArrowSP=blockArrowSP3, col=colorSet, border=borderSet));
}

#' Resize sp Polygon object
#'
#' Resize sp Polygon object
#'
resizePolygon <- function
(obj,
 size=c(1,1),
 offsetXY=c(0,0),
 ...)
{
   ## Purpose is to resize (rescale) a Polygon object from sp package.
   ## Ever so wonderful there isn't an R package to do this function,
   ## like, say, "sp".
   ##
   ## This method was obtained from:
   ## https://r-forge.r-project.org/scm/viewvc.php/branches/S3-classes/Rcaline/R/resize.R?view=markup&root=rcaline&pathrev=58
   if (!suppressPackageStartupMessages(require(sp))) {
      stop("The sp package is required.");
   }
   coords <- coordinates(obj);
   n <- nrow(coords);
   ## Make sure size will match the two dimensions
   #size <- rep(size, length.out=2);

   coords2 <- t(t(coords) * size + offsetXY);

   return(Polygon(coords2));
}

#' Resize sp Polygons object
#'
#' Resize sp Polygons object
#'
resizePolygons <- function
(obj,
 size=c(1,1),
 offsetXY=c(0,0),
 ...)
{
   ## Purpose is to resize (rescale) a Polygons object from sp package, and
   ## the nested internal Polygon objects.
   if (!suppressPackageStartupMessages(require(sp))) {
      stop("The sp package is required.");
   }
   poly <- lapply(obj@Polygons, function(x) {
      resizePolygon(x, size=size, offsetXY=offsetXY);
   });
   return(Polygons(poly, ID=obj@ID));
}

#' Resize sp SpatialPolygons object
#'
#' Resize sp SpatialPolygons object
#'
resizeSpatialPolygons <- function
(obj, size=c(1,1), offsetXY=c(0,0),
 ...)
{
   ## Purpose is to resize (rescale) a SpatialPolygons object, and its
   ## nested internal Polygons of Polygon objects.
   ##
   ## TODO: replace this function with something from rgeos
   if (!suppressPackageStartupMessages(require(sp))) {
      stop("The sp package is required.");
   }
   poly <- lapply(obj@polygons, function(x) {
      resizePolygons(x, size=size, offsetXY=offsetXY);
   });
   return(SpatialPolygons(poly,
      pO=obj@plotOrder,
      proj4string=CRS(proj4string(obj))));
}

#' Trim angle to a fixed range
#'
#' Trim angle to a fixed range
#'
trimAngle <- function
(x,
 minAngle=0,
 maxAngle=NULL,
 unit=c("degrees", "radians"),
 ...)
{
   ## Purpose is to trim angles to be between minAngle and maxAngle
   ## by subtracting the total. For example, to ensure angles
   ## are between 0 and 360 degrees, but any range will suffice.
   unit <- match.arg(unit);
   if (is.null(maxAngle)) {
      if (igrepHas("deg", unit)) {
         maxAngle <- 360;
      } else if (igrepHas("rad", unit)) {
         maxAngle <- pi*2;
      }
   }
   x %% maxAngle;
}
