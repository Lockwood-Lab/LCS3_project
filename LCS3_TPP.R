#### Analysis of Thermal proteome profiling data following Savitsky and Childs et al. (2018)
#### Author: Alvin, version: Aug 13, 2018
#### Uses TPP package version 3.9.3, all potential hits should have their melting curves manually inspected
### Download and install required packages 
source("http://bioconductor.org/biocLite.R")

#biocLite("TPP", dependencies=TRUE)
library(TPP)
citation("TPP")

library(limma)
library(plyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(xlsx)

#Parse raw data off of Proteome Discover in the sample format given
df_dmso <- read.table("isobarquant_dmso.txt", quote=NULL, header=TRUE, fill=TRUE, sep='\t')
df_lcs3 <- read.table("isobarquant_lcs3.txt", quote=NULL, header=TRUE, fill=TRUE, sep='\t')
df_dmso_2 <- read.table("isobarquant_dmso_2.txt", quote=NULL, header=TRUE, fill=TRUE, sep='\t')
df_lcs3_2 <- read.table("isobarquant_lcs3_2.txt", quote=NULL, header=TRUE, fill=TRUE, sep='\t')

#Removing proteins with values of 0 or duplicates
datacleanup <- function(datframe){
  datcolnames <- colnames(datframe)
  datframe <- ddply(datframe,'gene_name',colwise(median,datcolnames[2:length(colnames(datframe))]))
  datframe_noNA <- subset(datframe, datframe$rel_fc_126>0.6) 
  return(datframe_noNA)}

df_dmso_noNA <- datacleanup(df_dmso)
df_lcs3_noNA <- datacleanup(df_lcs3)
df_dmso_noNA_2 <- datacleanup(df_dmso_2)
df_lcs3_noNA_2 <- datacleanup(df_lcs3_2)



dfs <- list(df_dmso_noNA,df_dmso_noNA_2,df_lcs3_noNA,df_lcs3_noNA_2)
names(dfs) <- c("Vehicle_1","Vehicle_2","LCS3_1","LCS3_2")


### Configuration information, formatted according to the Childs et al. (2018) pdf
my_config <- read.xlsx("dup_LCS3_TPP-TR_config.xlsx", sheetName="Sheet 1", check.names=FALSE)
my_config <- my_config[,colSums(is.na(my_config))<nrow(my_config)] #remove all NA cols


resultPath = file.path(getwd(), 'Aug14_New_Graphics')
tpptrDefaultNormReqs()
normreq <- tpptrDefaultNormReqs() 

#Analyze our TPP data and generate QC info, set plotCurves to TRUE to get a list of all the melt curves for all of the proteins ID'd
TRresults <- analyzeTPPTR(configTable = my_config,
                          methods = "meltcurvefit",
                          data = dfs,
                          nCores = 4,
                          resultPath = resultPath,
                          plotCurves = FALSE)

write.table(TRresults,"TRresults.txt"
            ,quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)

####Generating median curves manually; post-normalization
#####
TR_median <- apply(TRresults[,2:41],2, median, na.rm = TRUE)
TR_mean <- apply(TRresults[,2:41],2,mean,na.rm=TRUE)
TR_median <- as.data.frame(TR_median)
TR_mean <- as.data.frame(TR_mean)
temperatures <- rep(c("37","41","44","47","50","53",'56','59','63','67'), times=4)
treatment <- rep(c("DMSO_1","DMSO_2","LCS3_1","LCS3_2"),each=10)
mediann <- cbind(TR_median, TR_mean, temperatures, treatment)
tplot <- qplot(data=mediann, x=temperatures, y=TR_mean, color=treatment, group=treatment)
tplot <- tplot + theme(axis.ticks.x=element_blank(), axis.title.y=element_text(size=18))
tplot
write.table(mediann,"Central_tendencies.txt",quote=FALSE,sep='\t',row.names=TRUE,col.names=NA)

#####

######Using Savitski's filters. Note that the Savitski group optimized these requirements to validate their own targets and is not the most statistically robust method (we prefer NPARC)
#####
tr_targets <- subset(TRresults, fulfills_all_4_requirements)$Protein_ID

#tr_targets <- c("GSTO1","TOE1")

tr_targetdf <- data.frame(tr_targets)
# write.table(tr_targetdf,"tr_targetdf.txt"
#             ,quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)


targets <- grep("",tr_targets,value=TRUE)
cat(targets, sep="\n")

## Step-by-step TPP analysis (and generate curves only for the targets)
trData <- tpptrImport(configTable = my_config, data = dfs)

normResults <- tpptrNormalize(data=trData)

normResults$normData$Vehicle_1@featureData
Biobase::featureNames(normResults$normData$Vehicle_1)

trData_target <- lapply(normResults$normData, function(d)
  d[Biobase::featureNames(d) %in% targets,])

resultPath2 = file.path(getwd(), 'Savitski_Hits') #What are our TPP hits based on the Savitski method

trData_target <- tpptrCurveFit(data=trData_target, resultPath=resultPath, nCores=1)

Biobase::pData(Biobase::featureData(trData_target[["Vehicle_1"]]))[,1:5]


#####Significance Assessment of Melting Pt Shifts

load(file.path(resultPath, "dataObj", "fittedData.RData"), verbose=TRUE)

minR2New <- 0.5 # instead of 0.8
maxPlateauNew <- 0.7 # instead of 0.3
newFilters <- list(minR2 = minR2New,
                   maxPlateau = maxPlateauNew)
TRresultsNew <- tpptrAnalyzeMeltingCurves(data = trDataFitted,
                                          pValFilter = newFilters)
tr_targetsNew <- subset(TRresultsNew, fulfills_all_4_requirements)$Protein_ID
targetsGained <- setdiff(tr_targetsNew,tr_targets)
targetsLost <- setdiff(tr_targets,tr_targetsNew)
cat(targetsGained, sep="\n")
print(targetsLost)

tppExport(tab = TRresultsNew,
          file = file.path(resultPath2, "targets_newFilters.xlsx"))

#### Custom filters 
#####
length(TRresults$Protein_ID) #5593 proteins identified across treatments and replicates to begin with

TRresults$min_slope_LCS3_2_vs_Vehicle_2 #We want the slope of the inflection to be steeper than -0.06
onepass <- filter(TRresults, min_slope_LCS3_1_vs_Vehicle_1< -0.06)
onepass <- filter(onepass, min_slope_LCS3_2_vs_Vehicle_2 < -0.06)
length(onepass$Protein_ID) #2797 passed

### R^2 fit of the normalization curve vs the raw data to be >0.8
twopass <- filter(onepass, passed_filter_LCS3_1_vs_Vehicle_1)
twopass <- filter(twopass, passed_filter_LCS3_2_vs_Vehicle_2)
length(twopass$Protein_ID) #1508  passed

### p-value filters (based on melting point differences)
threepass <- filter(twopass, twopass$pVal_adj_LCS3_1_vs_Vehicle_1<0.05)

length(threepass$Protein_ID) #75 passed

fourpass <- filter(threepass,threepass$pVal_adj_LCS3_2_vs_Vehicle_2<0.05)

length(fourpass$Protein_ID) #11 left in total

cat(fourpass$Protein_ID, sep="\n")


### Make a filtered full TRresults table
#####
fourpassoutput <- data.frame(fourpass$Protein_ID)
colnames(fourpassoutput) <- "Protein_ID"
filt_TRresults <- merge(TRresults,fourpassoutput, by="Protein_ID")

write.table(filt_TRresults,"filt_both_p0.05_TRresults.txt"
            ,quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)

#####


##### NPARC/spline fit, an additional method to identify potential candidates. This method is the only one that DOES NOT require the protein having been melted (all others require this to define a melting point)
#### However, these curves, as all do, must be visually inspected
##### 
resultPath3 = file.path(getwd(), 'Spline_curvesplotted')
TRresultsSp <- analyzeTPPTR(configTable = my_config,
                          methods = "splinefit",
                          data = dfs,
                          resultPath = resultPath3,
                          plotCurves = FALSE)

length(TRresultsSp$Protein_ID)

#p-values based on non-parametric fitting of the DMSO treated and LCS3 treated samples
splinepass <- filter(TRresultsSp, p_adj_NPARC< 0.01)
length(splinepass$Protein_ID)

splinepasscomp <-  splinepass[complete.cases(splinepass),]
length(splinepasscomp$Protein_ID)

splineoutput <- data.frame(splinepasscomp$Protein_ID)
colnames(splineoutput) <- "Protein_ID"
filt_TRresultsSp_comp <- merge(TRresultsSp,splineoutput, by="Protein_ID")

write.table(filt_TRresultsSp_comp,"filt_spline_TRresults_p0.01.txt"
            ,quote=FALSE,sep='\t',col.names=TRUE,row.names=FALSE)

#Get curves for spline fit 
targetspline <- (splinepasscomp$Protein_ID)
targetspline <- grep("",targetspline,value=TRUE) #grep all
cat(targetspline, sep="\n")

trData <- tpptrImport(configTable = my_config, data = dfs)

normResults <- tpptrNormalize(data=trData)

normResults$normData$Vehicle_1@featureData
Biobase::featureNames(normResults$normData$Vehicle_1)

trData_targetspline <- lapply(normResults$normData, function(d)
  d[Biobase::featureNames(d) %in% targetspline,])

resultPathsp2 = file.path(getwd(), 'target_splinep0.01')

trData_targetspline <- tpptrCurveFit(data=trData_targetspline, resultPath=resultPathsp2, nCores=4)

#####

##### Changing the graphics plotted by original TPP package, redefining package Functions
#####
#Redefine the original TPP functions
plotColors <- function(expConditions, comparisonNums){
  ## Create matching color pairs for vehicle and treatment groups.
  groupColors <- FALSE
  
  if (all(!is.na(expConditions))){
    condLevels <- sort(unique(expConditions))
    compLevels  <- sort(unique(comparisonNums))
    condNum <- length(condLevels)
    compNum  <- length(compLevels)
    
    ## Check if several experiments within a comparison got the same condition assigned
    ## (can happen, for example, when comparing two vehicle experiments against each other)
    if (max(table(paste(expConditions, comparisonNums))) == 1){
      ## Check if numbers of conditions and comparisons do not exceed maximum
      if (condNum==2 && compNum<=8  && compNum>0){ # brewer pal can only produce up to 8 color pairs with 2 intensities each
        if(all.equal(condLevels, c("Treatment", "Vehicle")) & all.equal(compLevels, 1:compNum)){
          groupColors <- TRUE
        }
      }
    }
  }
  
  plotClrsLight <-  c('#727272',"000000")                
  plotClrsDark  <-  c('#EF3535',"000000")  
 
  NAsInComparisonNums <- is.na(comparisonNums)
  numNAsInComparisonNums <- sum(NAsInComparisonNums)
  
  if (groupColors==TRUE){
    lightCols = colorRampPalette(plotClrsLight)(length(expConditions) + numNAsInComparisonNums)
    darkCols = colorRampPalette(plotClrsDark)(length(expConditions) + numNAsInComparisonNums)
    
    colorVec <- rep(NA, length(expConditions))
    for (r in compLevels){
      iT <- which(expConditions=="Treatment" & comparisonNums==r)
      iV <- which(expConditions=="Vehicle" & comparisonNums==r)
      colorVec[iT] <- darkCols[r]
      colorVec[iV] <- lightCols[r]
    }
    
    if(numNAsInComparisonNums > 0){
      colorVec[is.na(colorVec)] = darkCols[sum(!is.na(colorVec)) : length(darkCols)]
    }
    
  } else if (groupColors==FALSE){
    colorVec <- colorRampPalette(brewer.pal(n=8, name="Dark2"))(length(expConditions))
  }
  return(colorVec)
}

plotLineTypes <- function (expConditions) 
{
  typeVec <- rep(1, length(expConditions))
  typeVec[expConditions == "Vehicle"] <- 1
  return(typeVec)
}


plotMeltingCurve <- function (modelList, listUpper, listLower, xMat, fcMat, curvePars, 
          protID, filename, plotTheme, expConditions, expComps, addLegend, 
          useCI) 
{
  meltPoint <- NULL
  if (all(is.na(fcMat))) {
    return(NULL)
  }
  else {
    grNames <- names(modelList)
    yMax <- 1.25
    theme_set(plotTheme)
    xLen <- 100
    xMatLarge <- sapply(grNames, function(g) seq(from = min(xMat[g, 
                                                                 ]), to = max(xMat[g, ]), length.out = xLen))
    yPred <- sapply(grNames, function(g) robustNlsPredict(model = modelList[[g]], 
                                                          newdata = list(x = xMatLarge[, g])))
    if (useCI) {
      yUpper <- sapply(grNames, function(g) robustNlsPredict(model = listUpper[[g]], 
                                                             newdata = list(x = xMatLarge[, g])))
      yLower <- sapply(grNames, function(g) robustNlsPredict(model = listLower[[g]], 
                                                             newdata = list(x = xMatLarge[, g])))
    }
    compNumber <- assignCompNumber_to_expName(compDF = expComps, 
                                              expNames = grNames)
    plotCols <- plotColors(expConditions, compNumber)
    plotlTypes <- plotLineTypes(expConditions)
    names(plotCols) <- grNames
    names(plotlTypes) <- grNames
    grOrder <- c()
    if (all(!is.na(compNumber))) {
      for (r in unique(compNumber)) {
        iT <- which(expConditions == "Treatment" & compNumber == 
                      r)
        iV <- which(expConditions == "Vehicle" & compNumber == 
                      r)
        grOrder <- c(grOrder, grNames[iT], grNames[iV])
      }
    }
    if (length(grOrder) == 0) 
      grOrder <- grNames
    plotCols <- plotCols[grOrder]
    plotlTypes <- plotlTypes[grOrder]
    groupCol1 <- factor(rep(grNames, each = xLen), levels = grOrder)
    plotData1 <- data.frame(Group = groupCol1, Temperature = numeric(length(grNames) * 
                                                                       xLen), FoldChange = numeric(length(grNames) * xLen), 
                            CiUp = rep(NA_real_, length(grNames) * xLen), CiLow = rep(NA_real_, 
                                                                                      length(grNames) * xLen), DataType = "Model")
    for (g in grNames) {
      if (useCI) {
        plotData1[plotData1$Group == g, c("Temperature", 
                                          "FoldChange", "CiUp", "CiLow")] = cbind(xMatLarge[, 
                                                                                            g], yPred[, g], yUpper[, g], yLower[, g])
      }
      else {
        plotData1[plotData1$Group == g, c("Temperature", 
                                          "FoldChange")] = cbind(xMatLarge[, g], yPred[, 
                                                                                       g])
      }
    }
    groupCol2 <- factor(rep(grNames, each = ncol(fcMat)), 
                        levels = grOrder)
    plotData2 <- data.frame(Group = groupCol2, Temperature = numeric(length(grNames) * 
                                                                       ncol(fcMat)), FoldChange = numeric(length(grNames) * 
                                                                                                            ncol(fcMat)), DataType = "Measured")
    for (g in grNames) {
      plotData2[plotData2$Group == g, c("Temperature", 
                                        "FoldChange")] = cbind(xMat[g, ], fcMat[g, ])
    }
    xMP <- subset(curvePars, select = meltPoint)
    yMP <- sapply(grNames, function(g) {
      robustNlsPredict(modelList[[g]], newdata = list(x = xMP[g, 
                                                              ]))
    })
    groupCol3 <- factor(names(yMP), levels = grOrder)
    plotData3 <- data.frame(Group = groupCol3, yMP = yMP, 
                            xMP = xMP[, "meltPoint"])
    tableDF <- data.frame(condition = factor(grNames, levels = grOrder), 
                          meltPoint = round(curvePars[grNames, "meltPoint"], 
                                            2), slope = signif(curvePars[grNames, "slope"], 
                                                               2), plateau = round(curvePars[grNames, "plateau"], 
                                                                                   2), R2 = round(curvePars[grNames, "R_sq"], 2))
    p <- ggplot()
    p <- p + scale_color_manual(values = plotCols, breaks=c("Vehicle_1","Vehicle_2","LCS3_1","LCS3_2"), labels=c("Vehicle 1", "Vehicle 2","LCS3 1","LCS3 2"))
    p <- p + scale_linetype_manual(values = plotlTypes, breaks=c("Vehicle_1","Vehicle_2","LCS3_1","LCS3_2"), labels=c("Vehicle 1", "Vehicle 2","LCS3 1","LCS3 2"))
    p <- p + scale_y_continuous(limits = c(0, yMax), breaks=seq(0,2,0.5))
    if (addLegend) {
      p <- p + theme(legend.title = element_blank())
    }
    else {
      p <- p + theme(legend.position = "none")
    }
    p <- p + ggtitle(paste(protID,"\n",sep=""))  #+ theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    p <- p + xlab(paste("\nTemperature (°", "C)", sep = "")) + 
      ylab("Fraction non-denatured\n")
    if (useCI) {
      p <- p + geom_ribbon(data = plotData1, aes_string(x = "Temperature", 
                                                        ymax = "CiUp", ymin = "CiLow", fill = "Group"), 
                           alpha = 0.3) + scale_fill_manual(values = plotCols,breaks=c("Vehicle_1","Vehicle_2","LCS3_1","LCS3_2"), labels=c("Vehicle 1", "Vehicle 2","LCS3 1","LCS3 2"))
    }
    p <- p + geom_line(data = plotData1, size = 1, na.rm = TRUE, 
                       aes_string(x = "Temperature", y = "FoldChange", 
                                  colour = "Group", linetype = "Group"))
    p <- p + geom_point(data = plotData2, na.rm = TRUE, 
                        aes_string(x = "Temperature", y = "FoldChange", 
                                   colour = "Group"))
    p <- p + geom_point(data = plotData3, shape = 4, size = 5, 
                        show.legend = FALSE, na.rm = TRUE, aes_string(x = "xMP", 
                                                                      y = "yMP", colour = "Group"))
    p <- addTableToPlot(plotObj = p, tableDF = tableDF, 
                        meltVar = "condition", clrs = plotCols )
    return(p)
  }
}

#Reload the TPP package to ensure fidelity of the environment

detach("package:TPP", unload=TRUE)
library(TPP)

#Set the environment of the new user defined functions to the functions being replaced in the package

tmpfun <- get("plotMeltingCurve", 
              envir = asNamespace("TPP"))
environment(plotMeltingCurve) <- environment(tmpfun)

#Replace old functions with new functions temporarily 
assignInNamespace(x="plotColors", value="plotColors", ns=asNamespace("TPP"))
assignInNamespace(x="plotLineTypes", value="plotLineTypes", ns=asNamespace("TPP"))
assignInNamespace(x="plotMeltingCurve", value="plotMeltingCurve", ns=asNamespace("TPP"))

#Test curve fit to see if the graphics generated are optimal
trData_target <- tpptrCurveFit(data=dfs, resultPath=resultPath2, nCores=1)
#####  