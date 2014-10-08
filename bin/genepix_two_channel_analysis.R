#! /usr/bin/Rscript 
library('getopt');
library('marray');
library('limma');
library('convert');
library('nnNorm');
library('genefilter');

Sys.setlocale(locale="C");

# ./genepix_two_channel_analysis.R -i ~/workspace/walid/spruce_bud_set_roadmap/rawFiles -t ~/workspace/walid/spruce_bud_set_roadmap/in_Files/targets.txt -l ~/workspace/walid/spruce_bud_set_roadmap/galFiles/GQSR02-S1-0001-genepix.gal -o ~/workspace/walid/spruce_bud_set_roadmap

#get options, using the spec as defined by the enclosed list.
#we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
	'infile_dir','i', 1, "character",
	'target_infile','t', 1, "character",
	'gal_infile','l', 1, "character",
	'outfile_dir','o', 1, "character",
	'help','h', 0, "logical"
), byrow=TRUE, ncol=4);
opt = getopt(spec);

# if help was asked for print a friendly message
# and exit with a non-zero error code
if(!is.null(opt$help)){
	cat(getopt(spec, usage=TRUE));
	q(status=1);
}

# if no parameter was given for any of these print 
# a friendly message and exit with a non-zero error
# code
if(is.null(opt$infile_dir)){ 
	cat(getopt(spec, usage=FALSE));
	q(status=1);
}
if(is.null(opt$target_infile)){ 
	cat(getopt(spec, usage=FALSE));
	q(status=1);
}
if(is.null(opt$gal_infile)){ 
	cat(getopt(spec, usage=FALSE));
	q(status=1);
}
if(is.null(opt$outfile_dir)){ 
	cat(getopt(spec, usage=FALSE));
	q(status=1); 
}

# Initialize directory and file name variables
infile_dir <- opt$infile_dir;
target_infile <- opt$target_infile;
gal_infile <- opt$gal_infile;
outfile_dir <- opt$outfile_dir;

# generate base output directory if it does not exist.
dir.create(outfile_dir, showWarnings = TRUE);

# Make pre, old, and new normalization output directories
pre_norm_outfile_dir <- paste(outfile_dir, "pre_normalization", sep="/");
dir.create(pre_norm_outfile_dir, showWarnings = TRUE);

old_norm_outfile_dir <- paste(outfile_dir, "old_normalization", sep="/");
dir.create(old_norm_outfile_dir, showWarnings = TRUE);

new_norm_outfile_dir <- paste(outfile_dir, "new_normalization", sep="/");
dir.create(new_norm_outfile_dir, showWarnings = TRUE);

# Process input data files using BioConductor
target_dataframe  <- readTargets(target_infile);
num_target_files <- length(target_dataframe$FileName);


cat(paste("Processing", num_target_files, "target files for QA/QC Microarray Analysis...\n"));
flags <- function(x) as.numeric(x$Flags > -99);
RGList <- read.maimages(target_dataframe$FileName, path=infile_dir, "genepix", wt.fun=flags);
RGList$genes <- readGAL(gal_infile);
RGList$printer <- getLayout(RGList$genes);

# Obtain file names
target_filenames <- read.delim(target_infile)
colnames(RGList) <- target_filenames[1:num_target_files,1]

##
raw_data_table <- as(RGList,"marrayRaw");

# DIAGNOSTIC PLOTS
# 1. boxplot for all slides
png(file=paste(pre_norm_outfile_dir, "boxplot_raw.png", sep="/"));
par(mar=c(7,5,1,1))
maBoxplot(raw_data_table, ylab="maM", main="Boxplots on M-value by array before normalization", xlab="", las=2);
abline(h=0, col="grey", lwd=2);
dev.off();

for (slides in 1:num_target_files){
	if( slides < 10){
		index <- paste("0", slides, sep="")
	}
	else{
		index <- slides
	}
  
  #more boxplot (print-tips)
  png(file = paste(pre_norm_outfile_dir, "/", "boxplots_", index, "_" , colnames(RGList)[slides], ".png", sep=""))
  
  maBoxplot(raw_data_table[,slides], xlab="print-tip group", ylab="M (log2-scale)", main=paste("Boxplots on M-value by print-tip group before normalization for slide ", colnames(RGList)[slides], sep=""));
  abline(h=0, col="grey", lwd=2);
  dev.off();
  
  #2. 2D spatial plots
  png(file = paste(pre_norm_outfile_dir, "/", "spatialplots_Rb_", index, "_", colnames(RGList)[slides], ".png", sep=""));
  image(raw_data_table[,slides], xvar="maRb", bar=TRUE, main=paste("Spatial plots on background intensity of red channel for slide ", colnames(RGList)[slides], sep=""));
  dev.off();
  
  png(file = paste(pre_norm_outfile_dir, "/","spatialplots_Gb_", index, "_", colnames(RGList)[slides], ".png", sep=""));
  image(raw_data_table[,slides], xvar="maGb", bar=TRUE, main=paste("Spatial plots on background intensity of green channel for slide ", colnames(RGList)[slides], sep=""));
  dev.off();
  
  #3. MA plots
  png(file = paste(pre_norm_outfile_dir, "/","MAplots_",index,"_",colnames(RGList)[slides],".png",sep=""));
  limma::plotMA(RGList,
         array = slides,
         xlab = "A:Intensity (log2)",
         ylab = "M:Fold Change (log2)",
         main = paste("MAplot before normalization for slide ", colnames(RGList)[slides], sep=""),
         xlim = c(2,16),
         ylim = c(-6,6)
         );
  abline (h=0, col="red");
  dev.off();
}  

# save(RGList, file = file.path(pre_norm_outfile_dir, "RGList_data.txt"))
#more density plot
RGList.b <- backgroundCorrect(RGList,method="none");
png(file=paste(pre_norm_outfile_dir, "densityplots_raw_Overall.png", sep="/"));
plotDensities(RGList.b);
dev.off();

####################### QUALITY ACCESSMENT ON CONTROLS ############################
# pos.buffer <- grep("Buffer",RGList$gene[,4]);
# pos.empty <- grep("Empty", RGList$gene[,4]);
# pos.gfp <- grep("GFP", RGList$gene[,4]);
# pos.gfp1 <- grep("GFP1", RGList$gene[,4]);

# The array has 12 metarows. Each metarow contains 2916 probes
# Boxplot: variation of the control spots intensities per metarow
# png(file=paste(pre_norm_outfile_dir,"Control_buffer_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
# boxplot(MA$M[pos.buffer,1] ~ ((pos.buffer %/% 2916)),
#         main=paste("Buffer per metarow",sep=""),
#         xlab="Metarow",ylab="Log2(Buffer intensity)",
#         col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE"))
# dev.off()
# 
# png(file=paste(pre_norm_outfile_dir, "Control_empty_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
# boxplot(MA$M[pos.empty,1]~((pos.empty %/% 2916)),
#         main=paste("Empty per metarow",sep=""),
#         xlab="Metarow",ylab="Log2(Empty intensity)",
#         col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE"))
# dev.off()
# 
# png(file=paste(pre_norm_outfile_dir, "Control_GFP_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
# boxplot(MA$M[pos.gfp,1]~((pos.gfp %/% 2916)),
#         main=paste("GFP per metarow",sep=""),
#         xlab="Metarow",ylab="Log2(GFP intensity)",
#         col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE"))
# dev.off()
# 
# png(file=paste(pre_norm_outfile_dir, "Control_GFP1_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
# boxplot(MA$M[pos.gfp1,1]~((pos.gfp1 %/% 2916)),
#         main=paste("GFP1 per metarow",sep=""),
#         xlab="Metarow",ylab="Log2(GFP1 intensity)",
#         col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE"))
# dev.off()


#do some operation based on user input.
# cat(paste(infile_dir, target_infile, gal_infile, outfile_dir, old_norm_outfile_dir, new_norm_outfile_dir, pre_norm_outfile_dir, collapse="\n"));
# cat("\n");

#signal success and exit.
q(status=0);