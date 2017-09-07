volcanoPlot <- function
(x,
 n=NULL,

 ## Data columns
 geneColumn=head(provigrep(geneColGrep, colnames(x)),1),
 lfcCol=head(provigrep(fcColGrep, colnames(x)),1),
 pvCol=head(provigrep(pvColGrep, colnames(x)),1),
 intensityCols=provigrep(intensityColGrep, colnames(x)),

 ## grep expressions
 geneColGrep=c("^gene","symbol"),
 lfcColGrep=c("log.{0,2}fc", "lfc", "log.{0,4}fold", "fold.*ch", "fc", "fold", "log.*ratio", "ratio"),
 pvColGrep=c("^padj|adj.*p","fdr","q.*val","p.*val"),
 intensityColGrep=c("maxgroupmean|maxmean","ave.*expr","^fpkm"),

 ## Stats cutoffs
 pvalueCutoff=0.01,
 fcCutoff=0,
 intensityCutoff=NULL,
 minPvalue=1e-20,
 pvalueFloor=minPvalue,
 maxLog2FC=10,
 xRange=NULL,
 yRange=NULL,
 fcCeiling=NULL,
 fcCutoffInLog=FALSE,

 ## Plot features, title, axis labels, margins
 margins=c(5,6,4,2),
 main="Volcano Plot",
 lineMain=3,
 subMain=NULL,
 font.main=1,
 symmetricAxes=TRUE,
 xlab=paste("Fold change (", fcCol, ")", sep=""),
 ylab=paste("Significance (", pvCol, ")", sep=""),
 nXlabels=13,
 nYlabels=7,
 ablineColor="#000000AA",

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
 usePointColorSet=TRUE,

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
 #smoothColramp=c("white", "lightblue", "lightskyblue3", "royalblue", "darkblue", "orange", "darkorange1", "orangered2"),
 colramp=c("white", "lightblue", "lightskyblue3", "royalblue", "darkblue", "orange", "darkorange1", "orangered2"),
 useRaster=TRUE,

 doBoth=FALSE,
 bothColorOnlyHits=TRUE,
 orientation="up",
 addPlot=FALSE,
 pointSubset=NULL,
 plotOnlySubset=FALSE,
 transFactor=0.24,
 transformation=function(x){x^transFactor}, # use 0.14 for very large datasets
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
   ## pvalueFloor is designed to set extremely low P-values (1e-150) to a minimum
   ## sectionLabelSpacing is the fraction of the y-axis range away from each border
   ##    used to position the highlight hit counts labels (if highlightPoints is not NULL)
   ##
   ## A change was made to use pointColorSet to set the colors, instead of using
   ## hitColorSet, so the colors could independently be manipulated.
   ## To use this scheme, set usePointColorSet=TRUE
   ##
   ## The graph plots fold change on the x-axis, but with log2 scale, the numbers
   ## are technically being reported in normal space. However, to give it a log2 axis label:
   ## xlab=expression(log[2]*"(Fold change)"
   blockArrowCex <- rep(blockArrowCex, length.out=2);
   blockArrowCex <- rep(blockArrowCex, length.out=2);
   if (!usePointColorSet) {
      pointColorSet["up"] <- hitColorSet[2];
      pointColorSet["down"] <- hitColorSet[2];
      pointColorSet["base"] <- hitColorSet[1];
      pointColorSet["base.highlight"] <- hitColorSet[2];
      pointColorSet["up.highlight"] <- hitColorSet[2];
      pointColorSet["down.highlight"] <- hitColorSet[2];
      if (is.null(pointBgColors)) {
         pointAlphas <- (col2alpha(pointColorSet)*2+1)/3;
         if (verbose) {
            printDebug("pointAlphas: ", format(digits=2, pointAlphas), c("orange", "seagreen"));
         }
         pointDFs <- ifelse(col2hcl(pointColorSet)["L",] <= 60, -1.5, 1.5);
         pointBorderSet <- makeColorDarker(pointColorSet, fixAlpha=pointAlphas, darkFactor=pointDFs, sFactor=1.5);
      }
   } else {
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
   }
   ##
   if (doBoth && bothColorOnlyHits) {
      hitColorSet[1] <- "#FFFFFF00";
      pointColorSet["base"] <- hitColorSet[1];
   }

   if (verbose) {
      printDebug("fcCol: ", cPaste(fcCol), c("orange", "lightblue"));
      printDebug("pvCol: ", cPaste(pvCol), c("orange", "lightblue"));
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

   #pvHits <- length(which(x[1:n,pvCol] <= pvalueCutoff &
   #                       (x[1:n,pvCol] <= pvalueCutoff & (x[1:n,fcCol]) > fcCutoff |
   #                        x[1:n,pvCol] <= pvalueCutoff & -(x[1:n,fcCol]) > fcCutoff) ));
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
   pvHitsUpWhich <- which(x[,pvCol] <= pvalueCutoff & (x[,fcCol]) >= fcCutoff & metIntensity);
   pvHitsDownWhich <- which(x[,pvCol] <= pvalueCutoff & -(x[,fcCol]) >= fcCutoff & metIntensity);
   pvHitsBothWhich <- which(x[,pvCol] <= pvalueCutoff & abs(x[,fcCol]) >= fcCutoff & metIntensity);
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

   fcValues <- nameVector(x[,fcCol], rownames(x));

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
      printDebug("   volcanoPlot length(highlightPoints):", formatInt(length(highlightPoints)));
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
         printDebug("   volcanoPlot length(pointSubset):", formatInt(length(pointSubset)));
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
         pvHighlightHitsUp <- length(which(x[1:n,pvCol][pointSubset] <= pvalueCutoff & (x[1:n,fcCol][pointSubset]) >= fcCutoff));
         pvHighlightHitsDown <- length(which(x[1:n,pvCol][pointSubset] <= pvalueCutoff & -(x[1:n,fcCol][pointSubset]) >= fcCutoff));
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

   if (!is.null(fcCeiling)) {
      xValues1[abs(xValues1) > fcCeiling] <- sign(xValues1[abs(xValues1) > fcCeiling]) * fcCeiling;
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
         fcHitsWhich <- which(abs(x[1:n,fcCol]) >= fcCutoff);
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
            printDebug("par(mgp): ", cPaste(par("mgp")), c("orange", "lightblue"));
            printDebug("length(topHistBreaks): ", length(topHistBreaks), c("orange", "lightblue"));
         }
         r1 <- as.integer(length(topHistBreaks)/2);
         topHistCol <- rep(topBarColors, c(r1, r1+1));
         barplot(topHist$counts, axes=FALSE, ylim=c(0, max(topHist$counts)), space=0,
                 col=topHistCol, border=makeColorDarker(topHistCol), horiz=FALSE, las=2, cex.axis=1);
         prettyAt1 <- pretty(c(0, max(topHist$counts)), n=4);
         prettyAt1 <- prettyAt1[prettyAt1 <= max(topHist$counts)];
         axis(2, las=2, cex.axis=1.3, at=prettyAt1, ...);
         title(xlab=paste("Distribution of ", hitType, sep=""), cex.lab=topHistCexSub, xpd=TRUE, outer=FALSE, line=topHistSubLine);
         parList[["postTopHist"]] <- par();
         par("mar"=c(parMar[1], parMar[2], parMar[3]-1.5, parMar[4]));
      }
   }

   ######################
   ## Smooth scatter plot
   if (doSmoothScatter || doBoth) {
      if (orientation %in% c("right")) {
         ## Volcano plot on its side, tiled 90 degrees to the right
         smoothScatterFunc(y=xValues1, x=yValues1, #col=pvColors1, bg=pvBgColors1,
              xaxt="n", yaxt="n", xlab="", ylab="", transformation=transformation, useRaster=useRaster,
              xlim=yRange, ylim=xRange, nbin=nbin, colramp=smoothColramp, add=addPlot, nrpoints=nrpoints, ...);
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
              xlim=xRange, ylim=yRange, nbin=nbin, colramp=smoothColramp, add=addPlot, nrpoints=nrpoints, ...);
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
         pvHitsUp1 <- x[(which(x[1:n,pvCol] <= pvalueCutoff & (x[1:n,fcCol]) >= fcCutoff)),geneColumn];
         multiGenesUp <- tcountPaste(pvHitsUp1[pvHitsUp1 %in% names(tcount(pvHitsUp1)[tcount(pvHitsUp1) > 1])]);
         pvHitsDown1 <- x[(which(x[1:n,pvCol] <= pvalueCutoff & -(x[1:n,fcCol]) >= fcCutoff)),geneColumn];
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
         abline(v=-log10(pvalueCutoff), h=unique(c(-fcCutoff, fcCutoff)), lty="dashed", col=ablineColor);
         if (doBlockArrows) {
            hitCol <- hsv(h=0.06, s=0.75, v=0.9, alpha=1);
            if (doShadowTextArrows) {
               assign("text", get("shadowText"), envir=.GlobalEnv);
            }
            # blockArrowLabelHitColor="#FFFFAA", blockArrowLabelUpColor="#FFFFAA", blockArrowLabelDownColor="#FFFFAA"
            blockArrow(axisPosition="topAxis", xleft=-log10(pvalueCutoff), arrowPosition="right",
                        blockWidthPercent=4*blockArrowCex[2], arrowLengthPercent=4*blockArrowCex[2]/2,
                        arrowLabel=paste(format(big.mark=",", pvHits), " significant ", hitType, sep=""), col=hitCol,
                        labelCex=blockArrowLabelCex, labelFont=blockArrowFont, arrowLabelColor=blockArrowLabelHitColor,
                        arrowLabelBorder=alpha2col(setTextContrastColor(blockArrowLabelHitColor), 0.3));
            blockArrow(axisPosition="rightAxis", ybottom=fcCutoff, arrowPosition="top",
                        blockWidthPercent=4*blockArrowCex[1], arrowLengthPercent=4*blockArrowCex[1],
                        arrowLabel=labelUp, col="#990000FF", labelCex=blockArrowLabelCex,
                        labelFont=blockArrowFont, arrowLabelColor=blockArrowLabelUpColor,
                        arrowLabelBorder=alpha2col(setTextContrastColor(blockArrowLabelUpColor), 0.3));
            blockArrow(axisPosition="rightAxis", ytop=-fcCutoff, arrowPosition="bottom",
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
         abline(h=-log10(pvalueCutoff), v=unique(c(-fcCutoff, fcCutoff)), lty="dashed", col=ablineColor);
         if (doBlockArrows) {
            hitCol <- hsv(h=0.06, s=0.75, v=0.9, alpha=1);
            if (doShadowTextArrows) {
               assign("text", get("shadowText", envir=CBioRUtils), envir=.GlobalEnv);
            }
            if (doTopHist) {
               rightAdj <- 1.2;
            } else {
               rightAdj <- 1;
            }
            blockArrow(axisPosition="rightAxis", ybottom=-log10(pvalueCutoff), arrowPosition="top",
                        blockWidthPercent=4*blockArrowCex[2]*rightAdj, arrowLengthPercent=4*blockArrowCex[2]*rightAdj/2,
                        arrowLabel=paste(format(big.mark=",", pvHits), " significant ", hitType, sep=""), col=blockArrowHitColor,
                        labelCex=blockArrowLabelCex,
                        labelFont=blockArrowFont, arrowLabelColor=blockArrowLabelHitColor,
                        arrowLabelBorder=alpha2col(setTextContrastColor(blockArrowLabelHitColor), 0.3));
            blockArrow(axisPosition="topAxis", xleft=fcCutoff, arrowPosition="right",
                        blockWidthPercent=4*blockArrowCex[1], arrowLengthPercent=4*blockArrowCex[1],
                        arrowLabel=labelUp, col=blockArrowUpColor, labelCex=blockArrowLabelCex,
                        labelFont=blockArrowFont, arrowLabelColor=blockArrowLabelUpColor,
                        arrowLabelBorder=alpha2col(setTextContrastColor(blockArrowLabelUpColor), 0.3));
            blockArrow(axisPosition="topAxis", xright=-fcCutoff, arrowPosition="left",
                        blockWidthPercent=4*blockArrowCex[1], arrowLengthPercent=4*blockArrowCex[1],
                        arrowLabel=labelDown, col=blockArrowDownColor, labelCex=blockArrowLabelCex,
                        labelFont=blockArrowFont, arrowLabelColor=blockArrowLabelDownColor,
                        arrowLabelBorder=alpha2col(setTextContrastColor(blockArrowLabelDownColor), 0.3));
            if (doShadowTextArrows) {
               assign("text", graphics:::text, envir=.GlobalEnv);
            }
         }
      }
      parList[["postBlockArrows"]] <- par();

      ## Display the significance cutoff used
      #subTitle <- paste("Significance cutoff <= ", pvalueCutoff, ", and fold change cutoff > ", round(digits=2, 2^fcCutoff));
      if (!fcCutoffInLog) {
         fcCutoffLabel <- signif(digits=2, 2^fcCutoff);
      } else {
         fcCutoffLabel <- signif(digits=2, fcCutoff);
      }
      if (doStatsSubtitles) {
         if (fcCutoff > 0) {
            #subTitle <- paste(sep="", "Significance cutoff <= ", pvalueCutoff, ", and ",
            #                  ifelse(fcCutoffInLog, "log2 ", ""),
            #                  "fold change cutoff > ", fcCutoffLabel);
            subTitle <- paste(sep="", "Significance cutoff ≤ ", pvalueCutoff, ", and ",
                              ifelse(fcCutoffInLog, "log2 ", ""),
                              "fold change cutoff ≥ ", fcCutoffLabel);
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
         #par(parList[["postTopHist"]]);
         #par("mfg"=c(1,1));
         #par("mar"=parList[["postTopHist"]]$mar);
         #par("plt"=parList[["postTopHist"]]$plt);
         #par("usr"=parList[["postTopHist"]]$usr);
         parMfg(mfg=c(1,1), usePar=parList[["postTopHist"]]);
         ## For some reason, the title gets cropped unless I run it with outer=TRUE first,
         ## thereafter it isn't cropped, although no par() values change at all! Bug.
         title(main="", outer=TRUE, line=0);
         #printDebug(c("par('mar'): ", cPaste(format(digits=2, par('mar')))));
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
      #printDebug(c("pointSubset: ", cPaste(pointSubset)));
      #vp1 <- list(pointSubset=nameVector(pointSubset, names(xValues1[pointSubset])),
      #            pointSubsetX=nameVector(xValues1[pointSubset], names(xValues1[pointSubset])),
      #            pointSubsetY=nameVector(yValues1[pointSubset], names(xValues1[pointSubset])),
      #            pointSubsetCex=nameVector(pointSubsetCex, names(xValues1[pointSubset])),
      #            pointSubsetCol=nameVector(pvColors1[pointSubset], names(xValues1[pointSubset])));
      #vp1DF <- data.frame(lapply(vp1, function(i){i}));
      vp1DF <- data.frame(check.names=FALSE, stringsAsFactors=FALSE,
                          "pointSubset"=pointSubset,
                          "pointSubsetX"=xValues1[pointSubset],
                          "pointSubsetY"=yValues1[pointSubset],
                          "pointSubsetCex"=pointSubsetCex,
                          "pointSubsetCol"=pvColors1[pointSubset]);
      rownames(vp1DF) <- names(pointSubset);
      invisible(list(x=xValues1, y=yValues1, col=pvColors1,
                     pointSubsetDF=vp1DF,
                     parList=parList,
                     labelCoords=labelCoords));
   } else if (length(multiGenesUp) > 0 || length(multiGenesDown) > 0) {
      invisible(list(x=xValues1, y=yValues1, col=pvColors1,
                     multiGenesUp=multiGenesUp, multiGenesDown=multiGenesDown,
                     parList=parList,
                     labelCoords=labelCoords));
   } else {
      invisible(list(x=xValues1, y=yValues1, col=pvColors1,
                     parList=parList,
                     labelCoords=labelCoords));
   }
}
