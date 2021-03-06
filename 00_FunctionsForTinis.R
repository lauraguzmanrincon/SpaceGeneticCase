# ------------------------------ #
# Load basic functions, modified for Tinis
# Copied from /home/laura/Dropbox/Laura/PhD_Year2/04_AboutApproaches/RCode/2018_04/00_Functions.R
# Changes:  uploadPackages install in local library
#           add Tinis directory
# Copy made on: 22.11.2019
# ------------------------------ #

# 0. Load Workspace ----
#setwd("/home/laura/Dropbox/Laura/PhD_Year2/")
namePackages <- c("ggplot2",
                  "data.table")
sapply(namePackages, function(x) ifelse(!x %in% installed.packages(), install.packages(x), 0))
sapply(namePackages, function(x) library(x, character.only = TRUE))

# TINIS add local library to R libraries:
.libPaths(c(.libPaths(), "/home/maths/matsgk/Rlibrary"))

# 1. Save plots in PDF format ----
#' Save figures as pdfs
#' 25.08.2019: added the option of having a list of plots
savePDF <- function(fig, fileName, width = 7, height = 7, ..., figList = NULL){
  # Requires: dirPath
  out <- tryCatch(
    {
      if(is.null(figList)){
        pdf(file = paste(allPlotsPath, "/Plots/", fileName, ".pdf", sep = ""), width = width, height = height,...)
        print(fig)
        dev.off()
      }else{
        pdf(file = paste(allPlotsPath, "/Plots/", fileName, ".pdf", sep = ""), width = width, height = height,...)
        for(i in 1:length(figList)){
          print(figList[[i]])
        }
        dev.off()
      }
      
  }, error = function(e) cat(paste("Error in savePDF: plot not printed.\n", e)))
}

savePNG <- function(fig, fileName, width, height,...){
  # Requires: dirPath
  out <- tryCatch(
    {
      png(file = paste(allPlotsPath, "/Plots/", fileName, ".png", sep = ""), width = width, height = height,...)
      print(fig)
      dev.off()
    }, error = function(e) cat(paste("Error in savePNG: plot not printed.\n", e)))
}

# 2. Default ggplot theme for articles ----
# b. Lauras theme (Capp file ThemeFn.R) ----
theme_laura <- function(size = 11){
  theme_set(theme_grey(base_size = size))
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "gray", fill = NA))
}

theme_laura2 <- function(){
  theme(panel.grid.major = element_line(colour = "gray", size = 0.3)) + #"#9AB5B5"
    theme(panel.grid.minor = element_line(colour = "gray", size = 0.1)) + #"#C9D6D6"
    theme(panel.background = element_blank()) +
    theme(panel.border = element_rect(colour = "gray", fill = NA))
}

# 3. Multiplot ----
# Multiple plot function
# http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# 4. Read Nexus files from PubMLST ----
readNexusPubmlst <- function(pubmlstNexFile){
  # Load data
  X <- scan(file =pubmlstNexFile,
            what = character(), sep = "\n", quiet = TRUE, 
            comment.char = "[", strip.white = TRUE)
  startingLine <- which(sapply(1:length(X), function(x) grep("\\bmatrix\\b", X[x], ignore.case = TRUE)) == 1) + 1
  endingLine <- length(X) - 2
  endingLine - startingLine + 1 == numIsolates
  
  # Create two final lists
  linesList <- lapply(startingLine:endingLine, function(i) strsplit(X[i], "\t"))
  isolatesList <- data.table(idIsolateDist = 1:length(linesList),
                             isolate = as.numeric(lapply(linesList, function(x) strsplit(x[[1]], "\\|")[[1]][1])))
  distancesList <- lapply(1:length(linesList), function(i) cbind(i, 1:i, as.numeric(linesList[[i]][[1]][2:(i + 1)])))
  distancesMatrix <- as.data.table(do.call(rbind, distancesList))
  colnames(distancesMatrix) <- c("from", "to", "distance")
  distancesMatrix[, idDist := 1:nrow(distancesMatrix)]
  #return(list(isolatesList,distancesMatrix[from > to]))
  
  # (opt.) Output with only one list: translate idIsolateDist to isolates
  colnames(isolatesList)[colnames(isolatesList) == "isolate"] <- "isolateFrom"
  setkey(isolatesList, idIsolateDist)
  setkey(distancesMatrix, from)
  temp <- isolatesList[distancesMatrix]
  rm(distancesMatrix)
  colnames(isolatesList)[colnames(isolatesList) == "isolateFrom"] <- "isolateTo"
  setkey(isolatesList, idIsolateDist)
  setkey(temp, to)
  distancesMatrix <- isolatesList[temp]
  return(distancesMatrix[isolateFrom > isolateTo, .(isolateTo, isolateFrom, distance)])
}

# 5. New paste with another default ----
pasteL <- function(...){
  pastedText <- paste(..., sep = "")
  return(pastedText)
}

# 6. Upload packages ----
uploadPackages <- function(namePackages){
  sapply(namePackages, function(x) if(!x %in% installed.packages()){
    cat("Installing ", x, "...\n", sep = "")
    install.packages(x, lib = "/home/maths/matsgk/Rlibrary")
  })
  sapply(namePackages, function(x){
    cat("Loading ", x, "...\n", sep = "")
    library(x, character.only = TRUE)
  })
  # TODO fix error when "trying to use CRAN without setting a mirror" (check contrib.url)
  # Recall: /home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201908/InstallEpiclustRCmd.R
  #install.packages("epiclustR_dir/epiclustR_0.0.0.9000.tar.gz", repos = NULL, type = "source", lib = "/home/maths/matsgk/epiclustR_dir")
}

# 7. Plot matrix ----
uploadPackages("reshape")
plotMatrix <- function(matrixToPlot){
  colnames(matrixToPlot) <- NULL
  dimnames(matrixToPlot) <- NULL
  p <- ggplot(data = melt(matrixToPlot), aes(x = X1, y = X2, fill=value)) +
    geom_tile() +
    theme(panel.background = element_blank()) +
    theme(panel.grid.major = element_blank()) +
    theme(panel.grid.minor = element_blank())
  return(p)
}
