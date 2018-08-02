## chromstaR pipeline in R for Moqtaderi data as before but with reduced window size for increased resolution for tRNAs ##

library(chromstaR)
library(stringr)
library(dplyr)

#read in bam files and create experiment table
files <- list.files("../bam_realloc",pattern="Rpc155.*bam$")
exp = data.frame(file=files, mark="PolIII",condition = "K562", replicate = 1:2, pairedEndReads = FALSE)
#get input (control) file names and assign to samples in exp 
input_files <- list.files("../bam_realloc",pattern=".*Input.*bam$")
exp$controlFiles = paste(input_files[1], input_files[2], sep = "|")

#chromstaR
Chromstar(inputfolder="../bam_realloc", experiment.table=exp, outputfolder="./", mode="separate", numCPU=20, binsize=50)

#Adjust fdr to increase peak calling stringency
model = get(load(file.path("./","combined","combined_mode-separate.RData")))
model2 <- changeFDR(model, fdr = 1e-9)
#Adjust posterior cutoff instead of FDR - this adjusts for each bin and selects bins individually whereas changeFDR selects a whole peak if only 1 bin is significant
model_postCutoff_0.9 = changePostCutoff(model, post.cutoff = 0.9)
model_postCutoff_0.99 = changePostCutoff(model, post.cutoff = 0.99)
#compare number and width of peaks 
length(model$segments); mean(width(model$segments))
length(model2$segments); mean(width(model2$segments))
length(model_postCutoff_0.9$segments); mean(width(model_postCutoff_0.9$segments))
length(model_postCutoff_0.99$segments); mean(width(model_postCutoff_0.99$segments))
#export new peak files
exportPeaks(model2,"./BROWSERFILES/combined_mode-separate")
exportCombinations(model2,"./BROWSERFILES/combined_mode-separate")

exportPeaks(model_postCutoff_0.9, "./BROWSERFILES/combined_mode-separate-PostCutoff0.9")
exportCombinations(model_postCutoff_0.9, "./BROWSERFILES/combined_mode-separate-PostCutoff0.9")
exportPeaks(model_postCutoff_0.99, "./BROWSERFILES/combined_mode-separate-PostCutoff0.99")
exportCombinations(model_postCutoff_0.99, "./BROWSERFILES/combined_mode-separate-PostCutoff0.99")
#export count files
exportCounts(model2,"./BROWSERFILES/combined_mode-separate")
exportCounts(model_postCutoff_0.9, "./BROWSERFILES/combined_mode-separate-PostCutoff0.9")
exportCounts(model_postCutoff_0.99, "./BROWSERFILES/combined_mode-separate-PostCutoff0.99")

# manual chromstaR anaylsis to check replicate correlation
#change file paths back to full names
files <- list.files("../bam",pattern="PolIII.*bam$", full.names=TRUE)
input_files <- list.files("../bam",pattern=".*Input.*bam$", full.names=TRUE)

#bin data
binned.data = list()
for (file in files) {
	binned.data[[basename(file)]] = binReads(file, binsizes=1000, experiment.table=exp)
}

binned.data.input = list()
for (file in input_files) {
	binned.data.input[[basename(file)]] <- binReads(file, binsizes=1000)
}

#univariate peak calling
models <- list()
for (i1 in 1:length(binned.data)) {
	models[[i1]] <- callPeaksUnivariate(binned.data[[i1]], max.time=60, input.data = binned.data.input[[str_split_fixed(exp[match(names(binned.data), str_split_fixed(exp$file,"/",3)[,3]),c("controlFiles")][i1], "/",3)[,3]]], read.cutoff.quantile = 0.9) #apply read cutoff here
}

# multivariate peak calling and correlation check
multi.model =  callPeaksReplicates(models, max.time=60, eps=1, num.threads=20) #this just sits and does nothing!!

multi.model <- callPeaksMultivariate(models, use.states=states,eps = 0.1, num.threads = 20) #this works if read cutoff was applied in univariate peak calling
