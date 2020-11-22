
# code to be run manually (I just wrote it here because where else?)
# Like in 36c.R!!
if(0){
  # a. Paste files
  # CREATE FOLDER
  # (once) mkdir Rlibrary
  # (once) mkdir FilesFromLocal/38_Files
  # COPY FILES
  nameFiles <- "21112019"
  system(paste("scp", "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38_Files/38_*.R", # files
               #"/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38_Files/38_01_21112019_MCMCInputTinis.RData", # data
               "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201912/38_01_18122019_MCMCInput_GOX.RData", # data (see 38_00.R/0)
               "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_202002/38_01_05032020_MCMCInput_TGOX.RData", # data (see 38_00.R/0)
               "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38_Files/38_14i_ScriptTinis_sbatch.sbatch", # script
               "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38_Files/00_FunctionsForTinis.R", # functions
               "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files"))
  # or individually
  system(paste("scp", "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38_Files/38_*.R", # files
               "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files"))
  system(paste("scp", "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38_Files/38_14i_ScriptTinis_sbatch.sbatch", # script
               "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files"))
  system(paste("scp", "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38_Files/00_FunctionsForTinis.R", # functions
               "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files"))
  system(paste("scp", "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38_Files/38_01_21112019_MCMCInputTinis.RData", # data
               "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files"))
  system(paste("scp", "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38_Files/38_00_MCMC_SpaceGenomeTime.R", # 00 files
               "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files"))
  system(paste("scp", "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_202002/38_01_05032020_MCMCInput_TGOX.RData", # input data TG - interval1
               "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files"))
  system(paste("scp", "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_202002/38_01_29032020_MCMCInput_TGOX.RData", # input data TG - interval4
               "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files"))
  system(paste("scp", "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_202003/38_01_30032020_MCMCInput_TG2OX_I4.RData", # input data TG2 - interval1 / 4
               "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files"))
  
  # Also simulated data
  system(paste("scp", "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38_Files/38_13_21112019_SimulatedData.RData", # data simul (NO)
               "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files")) # TODO borrar from Tinis
  system(paste("scp", "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38_Files/38_13_21112019_SimulatedData.RData", # data simul processed
               "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files"))
  
  # Also hclust information for G block update 26.02.2020
  system(paste("scp", "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38V2_II_18122019_ClustersKInfo_OX.RData",
               "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files"))
  # plus other hclust files
  system(paste("scp", "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_202004/38V2_II_14042020_*",
               "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files"))
  # plus input files for ST
  system(paste("scp", "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_202004/38_01_14042020_FilesForST/*",
               "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_01_14042020_FilesForST"))
  # plus TWSG and ALLTG input files
  system(paste("scp", "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_202004/38_01_*.RData*",
               "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files"))
  
  # RUN
  # sbatch -p cnode FilesFromLocal/38_Files/38_14i_ScriptTinis_sbatch.sbatch
  
  # Bring back results
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_21112019_MCMCOutputTinis.RData",
               "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911/38_Files"))
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_21112019_MCMCOutputTinis.RData",
               "/home/laura/Dropbox/Laura/PhD_Year3/07_MixedModelsP2/RCode_201911"))
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_21112019_SimulOutputTinis*",#_a
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_201911")) # output simulated per parameter!
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_03122019_MCMCOutputTinis.RData",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_201911")) # 5000 iterations
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_0601202001*",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202001")) # 5000 iterations, three kernels, inputdata 18122019
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_0701202001_MCMCOutputTinis_blockBeta.RData",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202001")) # 500 iterations, MAT kernels, Bijk, inputdata 18122019
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_2601202003_Tinis*.RData",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202001")) # 5000 iterations, SE/MAT12/MAT32 kernels, Bijk, inputdata 18122019
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_2602202001_Tinis_SE_GBlock.RData",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202002")) # 5000 iterations, SE/MAT12 kernel, Bk, inputdata 18122019, G block update
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_2602202001_Tinis_MAT12_GBlock_rho.RData",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202002")) # 5000 iterations, MAT12 kernel, Bk, inputdata 18122019, rho adaptative
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_0203202002_Tinis_MAT12_GBlock_rho_two.RData",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202002")) # 5000 iterations, MAT12 kernel, Bk, inputdata 18122019, rho adaptative two updates
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_03032020*.RData",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202002")) # 5000 iterations, MAT12 kernel, Bk, inputdata 18122019, rho adaptative two updates
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_0503202001_TinisTG_MAT12.RData",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202002")) # 5000 iterations, MAT12 kernel, Bk, inputdata 05032020, TG
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_1403202001_TinisSG_MAT12_1000it_cuts.RData",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202003")) # 5000 iterations, MAT12 kernel, Bk, inputdata 14032020, SG - cut analysis
  
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_1703202001_TinisTG_MAT12_500it_Sblock.RData",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202003")) # 500 iterations, MAT12 kernel, Bk, inputdata 05032020, TG with S block updates
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/BORRAR",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202004"))
  # Bring back ST results 15.04.2020
  system(paste("scp -r", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_01_14042020_FilesOutForST",
              "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202004/")) # ST 2000it, 200 burnin, for ROC, input in RCode202004/38_01_14042020_FilesForST
  # Bring back TG/TG2 ALL results 15.04.2020
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_21042020_TinisTG_ALL_interval1_autocorr_1000it_gibbsP.RData",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202004/")) # TG ALL 1000it, 200 burnin
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_15042020_TinisSG_NE_MAT12_5000it.RData",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202004/")) # SG TW 5000it, 200 burnin
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_21042020_TinisTG_ALL_interval1_autocorr_1000it_lessOutbreaks.RData",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202004/"))
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_23042020_TinisTG_ALL_interval1_nocorr_1000it_lessOutbreaks.RData",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202004/"))
  system(paste("scp ", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_00_26042020*ZERO.RData",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202004/")) # new resultsTG/TG2 TEMP
  system(paste("scp -r", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_01_05072020_FilesOutForST",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202006/")) # ST 5000it, 500 burnin, ChST IM, input in RCode202006/38_01_05072020_FilesOutForST
  system(paste("scp -r", "matsgk@tinis.csc.warwick.ac.uk:/home/maths/matsgk/FilesFromLocal/38_Files/38_01_05072020_FilesOutForST/38_00_050720200*",
               "/home/laura/Documents/PhD_Year3_NoDropbox/07_MixedModelsP2/RCode_202006/38_01_05072020_FilesOutForSTIM/")) # For the traces
}

# b. Load Workspace
source("FilesFromLocal/38_Files/00_FunctionsForTinis.R")
namePackages <- c("ggplot2", "data.table", "scales", "reshape", #"ggdendro",
                  "MASS", # for multivariate normal sampling
                  "matrixcalc", # is.positive.definite
                  "expm",  # square of a matrix
                  "spam", # SPArseMatrix
                  "devtools") # many things. Here: install from Github
uploadPackages(namePackages)

# c. Create/load data ----
if(simulatedData == 0){
  #load("FilesFromLocal/38_Files/38_01_21112019_MCMCInputTinis.RData") # CONFIG: inputForData OR notsim, dims23, interval3, beta3, cuts50/300, regionOXMSOA
  if(typeSpatial == "OX" & typeModel == "SG" & numClustSpatial == "MSOA"){
    load("FilesFromLocal/38_Files/38_01_18122019_MCMCInput_GOX.RData") # CONFIG: inputForData OR notsim, dims23, interval1?, beta3, cuts10/50, regionOXMSOA
  }else if(typeSpatial == "OX" & typeModel == "TG" & lengthBlockPeriod == 1){
    load("FilesFromLocal/38_Files/38_01_05032020_MCMCInput_TGOX.RData") # CONFIG: inputForData OR notsim,dims13,interval1,beta3,cuts10/50,regionOXMSOA
  }else if(typeSpatial == "OX" & typeModel == "TG" & lengthBlockPeriod == 4){
    load("FilesFromLocal/38_Files/38_01_29032020_MCMCInput_TGOX.RData") # CONFIG: inputForData OR notsim,dims13,interval4,beta3,cuts10/50,regionOXMSOA
  }else if(typeSpatial == "OX" & typeModel == "TG2" & lengthBlockPeriod == 1){
    load("FilesFromLocal/38_Files/38_01_30032020_MCMCInput_TG2OX.RData") # CONFIG: inputForData OR notsim,dims13,interval1,beta3,cuts10/50,regionOXMSOA
  }else if(typeSpatial == "OX" & typeModel == "TG2" & lengthBlockPeriod == 4){
    load("FilesFromLocal/38_Files/38_01_30032020_MCMCInput_TG2OX_I4.RData") # CONFIG: inputForData OR notsim,dims13,interval4,beta3,cuts10/50,regionOXMSOA
  }else if(typeModel == "ST"){
    load(paste0("FilesFromLocal/38_Files/38_01_14042020_FilesForST/38_01_14042020_MCMCInput_ST",
                typeSpatial,"_", numClustSpatial, "_", lengthBlockPeriod,".RData")) # CONFIG: inputForData OR notsim,dims12,interval[*],beta2,cuts10/50,regionOX[*]
  }else if(typeSpatial == "TW" & typeModel == "SG" & numClustSpatial == "MSOA"){
    load("FilesFromLocal/38_Files/38_01_14042020_MCMCInput_SGTW.RData") # CONFIG: inputForData OR notsim,dims23,interval1,beta3,cuts10/50,regionTWMSOA
  }else if(typeSpatial == "ALL" & typeModel == "TG" & lengthBlockPeriod == 1){
    load("FilesFromLocal/38_Files/38_01_14042020_MCMCInput_TGALL.RData") # CONFIG: inputForData OR notsim,dims13,interval1,beta3,cuts10/50,regionALLMSOA
  }else if(typeSpatial == "ALL" & typeModel == "TG2" & lengthBlockPeriod == 1 & sum(!c(21, 353, 45) %in% stGeneticGroups) == 0){
    load("FilesFromLocal/38_Files/38_01_15042020_MCMCInput_TG2ALL_I1.RData") # CONFIG: inputForData OR notsim,dims13,interval1,beta3,cuts10/50,regionALLMSOA
  }else{
    stop("Input file does not exist")
  }
  if(typeSpatial == "OX"){
    load("FilesFromLocal/38_Files/38V2_II_18122019_ClustersKInfo_OX.RData") # add clustering for G block update! 26.02.2020
    # hClustOut, ddata, readme, exploreCutFn, clusterInfoFn, colsDistanceMatrixGroupsNoDups*
  }else if(typeSpatial == "TW"){
    load("FilesFromLocal/38_Files/38V2_II_14042020_ClustersKInfo_TW.RData")
    # hClustOut, ddata, readme, exploreCutFn, clusterInfoFn, colsDistanceMatrixGroupsNoDups
  }else if(typeSpatial == "ALL"){
    load("FilesFromLocal/38_Files/38V2_II_14042020_ClustersKInfo_ALL.RData")
    # hClustOut, ddata, readme, exploreCutFn, clusterInfoFn, colsDistanceMatrixGroupsNoDups
  }
}else{
  # note simulatedData is stored as 0 in the first file!
  load("FilesFromLocal/38_Files/38_01_21112019_MCMCInputTinis.RData") # CONFIG: notsim, dims23, interval3, beta3, cuts50/300, regionOXMSOA
  load("FilesFromLocal/38_Files/38_13_21112019_SimulatedData.RData") # yy, simulatedParam
  load("FilesFromLocal/38_Files/38V2_II_18122019_ClustersKInfo_OX.RData") # add clustering for G block update! 26.02.2020
  # hClustOut, ddata, readme, exploreCutFn, clusterInfoFn, colsDistanceMatrixGroupsNoDups*
  y <- yy
  simulatedData <- 1
}



