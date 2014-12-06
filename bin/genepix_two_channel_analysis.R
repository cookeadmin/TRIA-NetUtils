#! /usr/bin/Rscript 
library('getopt');
library('limma');
library('marray');
library('convert');
library('nnNorm');
library('genefilter');

Sys.setlocale(locale="C");

# ./genepix_two_channel_analysis.R -i ~/workspace/walid/S038/INPUT_FILES -t ~/workspace/walid/S038/S038-targets.txt -l ~/workspace/walid/S038/GQSR02-S1-0001-genepix.gal -o ~/workspace/walid/S038

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
dir.create(outfile_dir, showWarnings = FALSE);

# Make before and after normalization output directories
before_normalization_dir <- paste(outfile_dir, "BEFORE_NORMALIZATION_DIR", sep="/");
dir.create(before_normalization_dir, showWarnings = FALSE);

after_normalization_dir <- paste(outfile_dir, "AFTER_NORMALIZATION_DIR", sep="/");
dir.create(after_normalization_dir, showWarnings = FALSE);

# Process input data files using BioConductor
target_dataframe  <- readTargets(target_infile);
num_target_files <- length(target_dataframe$FileName);

cat(paste("Processing", num_target_files, "target files for QA/QC Microarray Analysis...\n"));
RGList = read.maimages(target_dataframe$FileName, path=infile_dir, "genepix", wt.fun=wtflags(0.1));

RGList$genes <- readGAL(gal_infile);
RGList$printer <- getLayout(RGList$genes);

# Obtain file names
target_filenames <- read.delim(target_infile);
colnames(RGList) <- target_filenames[1:num_target_files,1];

# NORMALIZATION METHODS
# Before normalization

MAList <- normalizeWithinArrays(RGList, method="none", bc.method="none");
raw_data_table <- as(RGList, "marrayRaw");

# Apply after normalization method: ( withinarray = print-tip loess; between array = quantile norm)
data_table_after_norm <- normalizeWithinArrays(RGList, method="printtiploess");
data_table_scale_after <- normalizeBetweenArrays(data_table_after_norm, method="quantile");
data_table_after_norm_marray <- as(data_table_scale_after, "marrayNorm");

# Write normalized data to file
after_norm_data <- cbind(data_table_scale_after$genes$ID, data_table_scale_after$M);
colnames(after_norm_data)[1] <- "ID";

write.table(after_norm_data, file=paste(after_normalization_dir, "normalized_data.txt", sep="/"));

# DIAGNOSTIC PLOTS
# 1. boxplot for all slides
png(file=paste(before_normalization_dir, "boxplot_raw.png", sep="/"), width=800, height=600, res=NA, units="px");
par(mar=c(10,5,1,1));
maBoxplot(raw_data_table, ylab="maM", main="Boxplots on M-value by array before normalization", xlab="", las=2);
abline(h=0, col="grey", lwd=2);
dev.off();

png(file = paste(after_normalization_dir, "boxplot_within_slide.png", sep="/"), width=800, height=600, res=NA, units="px");
par(mar=c(10,5,1,1));
maBoxplot(as(data_table_after_norm,"marrayNorm"), ylab="maM", main="Boxplots on M-value by array after within-array normalization", xlab="", las=2);
abline(h=0, col="grey", lwd=2);
dev.off();

png(file = paste(after_normalization_dir, "boxplot_between_slides.png", sep="/"), width=800, height=600, res=NA, units="px");
par(mar=c(10,5,1,1));
maBoxplot(data_table_after_norm_marray, ylab="maM", main="Boxplots on M-value by array after between-array normalization", xlab="", las=2);
abline(h=0, col="grey", lwd=2);
dev.off();

for (slides in 1:num_target_files){
	if( slides < 10){
		index <- paste("0", slides, sep="");
	}
	else{
		index <- slides;
	}
  
	#more boxplot (print-tips)
	png(file = paste(before_normalization_dir, "/", paste("boxplots", index, colnames(RGList)[slides], sep="_"), ".png", sep=""), width=800, height=600, res=NA, units="px");
	maBoxplot(raw_data_table[,slides], xlab="print-tip group", ylab="M (log2-scale)", main=paste("Boxplots on M-value by print-tip group before normalization for slide", colnames(RGList)[slides], sep=" "));
	abline(h=0, col="grey", lwd=2);
	dev.off();
	
	png(file = paste(after_normalization_dir,"/", paste("boxplots", index, colnames(RGList)[slides], sep="_"), ".png", sep=""), width=800, height=600, res=NA, units="px");
	maBoxplot(data_table_after_norm_marray[,slides], xlab="print-tip group", ylab="M (log2-scale)", main= paste("Boxplots on M-value by print-tip group after normalization for slide", colnames(RGList)[slides], sep=" "));
	abline(h=0, col="grey", lwd=2);
	dev.off();

	
	#2. 2D spatial plots
	png(file = paste(before_normalization_dir, "/", paste("spatialplots_Rb", index, colnames(RGList)[slides], sep="_"), ".png", sep=""), width=800, height=600, res=NA, units="px");
	image(raw_data_table[,slides], xvar="maRb", bar=TRUE, main=paste("Spatial plots on background intensity of red channel for slide", colnames(RGList)[slides], sep=" "));
	dev.off();

	png(file = paste(before_normalization_dir, "/", paste("spatialplots_Gb", index, colnames(RGList)[slides], sep="_"), ".png", sep=""), width=800, height=600, res=NA, units="px");
	image(raw_data_table[,slides], xvar="maGb", bar=TRUE, main=paste("Spatial plots on background intensity of green channel for slide", colnames(RGList)[slides], sep=" "));
	dev.off();

	#3. MA plots
	png(file = paste(before_normalization_dir, "/", paste("MAplots", index, colnames(RGList)[slides], sep="_"), ".png", sep=""), width=800, height=600, res=NA, units="px");
	limma::plotMA(RGList,
		array = slides,
		xlab = "A:Intensity (log2)",
		ylab = "M:Fold Change (log2)",
		main = paste("MAplot before normalization for slide", colnames(RGList)[slides], sep=" "),
		xlim = c(2,16),
		ylim = c(-6,6)
	);
	abline (h=0, col="red");
	dev.off();
	
	png(file = paste(after_normalization_dir, "/", paste("MAplots", index, colnames(RGList)[slides], sep="_"), ".png", sep=""), width=800, height=600, res=NA, units="px");
	limma::plotMA(data_table_scale_after,
		array = slides,
		xlab = "A:Intensity (log2)",
		ylab = "M:Fold Change (log2)",
		main = paste("MAplot after normalization for slide", colnames(RGList)[slides], sep=" "),
		xlim = c(2,16),
		ylim = c(-6,6)
	);
	abline (h=0, col="red");
	dev.off();

	png(file=paste(before_normalization_dir, "/", paste("densityplots_raw", index, colnames(RGList)[slides], sep="_"), ".png", sep=""), width=800, height=600, res=NA, units="px");
	npoint <- 512;
	colorPalate <- c("red","green");

	RGList.subtractBC <- backgroundCorrect(RGList[,slides], method="subtract");

	RIntensities <- log2(RGList.subtractBC$R);
	GIntensities <- log2(RGList.subtractBC$G);

	RIntensities <- na.omit(RIntensities);
	GIntensities <- na.omit(GIntensities);

	densitiesRX <- densitiesRY <- matrix(0, npoint, ncol(RIntensities));
	for (column in 1:ncol(RIntensities)){
		densitiesR <- density(RIntensities[,column], n=npoint);
		densitiesRX[,column] <- densitiesR$x;
		densitiesRY[,column] <- densitiesR$y;
	}

	densitiesGX <- densitiesGY <- matrix(0, npoint, ncol(GIntensities));
	for (column in 1:ncol(GIntensities)){
		densitiesG <- density(GIntensities[,column], n=npoint);
		densitiesGX[,column] <- densitiesG$x;
		densitiesGY[,column] <- densitiesG$y;
	}

	plot.new();
	plot.window(xlim=range(cbind(densitiesRX, densitiesGX)), ylim=range(cbind(densitiesRY, densitiesGY)), xaxs ="i");
	box();
	axis(1);
	axis(2);

	title(main=paste("Density Plot before within-array normalization for slide", colnames(RGList)[slides], sep=" "), xlab="Intensity", ylab="Density");
	for (column in 1:ncol(RIntensities)){
		lines(densitiesRX[,column], densitiesRY[,column], lty="solid", lwd=2, col="red");
	}

	for (column in 1:ncol(GIntensities)){
		lines(densitiesGX[,column], densitiesGY[,column], lty="solid", lwd=2, col="green");
	}
	dev.off();

	png(file=paste(after_normalization_dir, "/", paste("densityplots_within_array", index, colnames(RGList)[slides], sep="_"), ".png", sep=""), width=800, height=600, res=NA, units="px");
	npoint <- 512;
	colorPalate <- c("red","green");
	RIntensities <- data_table_after_norm[,slides]$A+data_table_after_norm[,slides]$M/2;
	GIntensities <- data_table_after_norm[,slides]$A-data_table_after_norm[,slides]$M/2;

	RIntensities <- na.omit(RIntensities);
	GIntensities <- na.omit(GIntensities);

	densitiesRX <- densitiesRY <- matrix(0, npoint, ncol(RIntensities));
	for (column in 1:ncol(RIntensities)){
		densitiesR <- density(RIntensities[,column], n=npoint);
		densitiesRX[,column] <- densitiesR$x;
		densitiesRY[,column] <- densitiesR$y;
	}

	densitiesGX <- densitiesGY <- matrix(0, npoint, ncol(GIntensities));
	for (column in 1:ncol(GIntensities)){
		densitiesG <- density(GIntensities[,column], n=npoint);
		densitiesGX[,column] <- densitiesG$x;
		densitiesGY[,column] <- densitiesG$y;
	}

	plot.new();
	plot.window(xlim=range(cbind(densitiesRX, densitiesGX)), ylim=range(cbind(densitiesRY, densitiesGY)), xaxs ="i");
	box();
	axis(1);
	axis(2);
	title(main=paste("Density Plot after within-array normalization for slide", colnames(RGList)[slides], sep=" "), xlab="Intensity", ylab="Density");

	for (column in 1:ncol(RIntensities)){
		lines(densitiesRX[,column], densitiesRY[,column], lty="solid", lwd=2, col="red");
	}

	for (column in 1:ncol(GIntensities)){
		lines(densitiesGX[,column], densitiesGY[,column], lty="solid", lwd=2, col="green");
	}
	dev.off();
}  

#more density plot
png(file=paste(before_normalization_dir, "densityplots_raw_overall.png", sep="/"), width=800, height=600, res=NA, units="px");
npoint <- 512;
colorPalate <- c("red","green");

RGList.subtractBC <- backgroundCorrect(RGList, method="subtract");

RIntensities <- log2(RGList.subtractBC$R);
GIntensities <- log2(RGList.subtractBC$G);

RIntensities <- na.omit(RIntensities);
GIntensities <- na.omit(GIntensities);

densitiesRX <- densitiesRY <- matrix(0, npoint, ncol(RIntensities));
for (column in 1:ncol(RIntensities)){
	densitiesR <- density(RIntensities[,column], n=npoint);
	densitiesRX[,column] <- densitiesR$x;
	densitiesRY[,column] <- densitiesR$y;
}

densitiesGX <- densitiesGY <- matrix(0, npoint, ncol(GIntensities));
for (column in 1:ncol(GIntensities)){
	densitiesG <- density(GIntensities[,column], n=npoint);
	densitiesGX[,column] <- densitiesG$x;
	densitiesGY[,column] <- densitiesG$y;
}

plot.new();
plot.window(xlim=range(cbind(densitiesRX, densitiesGX)), ylim=range(cbind(densitiesRY, densitiesGY)), xaxs ="i");
box();
axis(1);
axis(2);

title(main="Density Plot before within-array normalization", xlab="Intensity", ylab="Density");
for (column in 1:ncol(RIntensities)){
	lines(densitiesRX[,column], densitiesRY[,column], lty="solid", lwd=2, col="red");
}

for (column in 1:ncol(GIntensities)){
	lines(densitiesGX[,column], densitiesGY[,column], lty="solid", lwd=2, col="green");
}
dev.off();

png(file=paste(after_normalization_dir, "densityplots_within_array_overall.png", sep="/"), width=800, height=600, res=NA, units="px");
npoint <- 512;
colorPalate <- c("red","green");
RIntensities <- data_table_after_norm$A+data_table_after_norm$M/2;
GIntensities <- data_table_after_norm$A-data_table_after_norm$M/2;

RIntensities <- na.omit(RIntensities);
GIntensities <- na.omit(GIntensities);

densitiesRX <- densitiesRY <- matrix(0, npoint, ncol(RIntensities));
for (column in 1:ncol(RIntensities)){
	densitiesR <- density(RIntensities[,column], n=npoint);
	densitiesRX[,column] <- densitiesR$x;
	densitiesRY[,column] <- densitiesR$y;
}

densitiesGX <- densitiesGY <- matrix(0, npoint, ncol(GIntensities));
for (column in 1:ncol(GIntensities)){
	densitiesG <- density(GIntensities[,column], n=npoint);
	densitiesGX[,column] <- densitiesG$x;
	densitiesGY[,column] <- densitiesG$y;
}

plot.new();
plot.window(xlim=range(cbind(densitiesRX, densitiesGX)), ylim=range(cbind(densitiesRY, densitiesGY)), xaxs ="i");
box();
axis(1);
axis(2);
title(main="Density Plot after within-array normalization", xlab="Intensity", ylab="Density");

for (column in 1:ncol(RIntensities)){
	lines(densitiesRX[,column], densitiesRY[,column], lty="solid", lwd=2, col="red");
}

for (column in 1:ncol(GIntensities)){
	lines(densitiesGX[,column], densitiesGY[,column], lty="solid", lwd=2, col="green");
}
dev.off();

# QUALITY ACCESSMENT ON CONTROLS
pos.buffer <- grep("Buffer",RGList$gene[,4]);
pos.empty <- grep("Empty", RGList$gene[,4]);
pos.gfp <- grep("GFP", RGList$gene[,4]);
pos.gfp1 <- grep("GFP1", RGList$gene[,4]);

# The array has 12 metarows. Each metarow contains 2916 probes
# Boxplot: variation of the control spots intensities per metarow
png(file=paste(before_normalization_dir, "control_buffer_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(MAList$M[pos.buffer,1] ~ ((pos.buffer %/% 2916)),
	main="Buffer per metarow",
	xlab="Metarow", ylab="Log2(Buffer intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(after_normalization_dir, "control_buffer_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(data_table_scale_after$M[pos.buffer,1] ~ ((pos.buffer %/% 2916)),
	main="Buffer per metarow",
	xlab="Metarow", ylab="Log2(Buffer intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(before_normalization_dir, "control_empty_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(MAList$M[pos.empty,1] ~ ((pos.empty %/% 2916)),
	main="Empty per metarow",
	xlab="Metarow",ylab="Log_2 (Empty intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(after_normalization_dir, "control_empty_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(data_table_scale_after$M[pos.empty,1] ~ ((pos.empty %/% 2916)),
	main="Empty per metarow",
	xlab="Metarow",ylab="Log_2 (Empty intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(before_normalization_dir, "control_GFP_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(MAList$M[pos.gfp,1] ~ ((pos.gfp %/% 2916)),
	main="GFP per metarow",
	xlab="Metarow",ylab="Log_2 (GFP intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(after_normalization_dir, "control_GFP_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(data_table_scale_after$M[pos.gfp,1] ~ ((pos.gfp %/% 2916)),
	main="GFP per metarow",
	xlab="Metarow",ylab="Log_2 (GFP intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();


png(file=paste(before_normalization_dir, "control_GFP1_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(MAList$M[pos.gfp1,1] ~ ((pos.gfp1 %/% 2916)),
	main="GFP1 per metarow",
	xlab="Metarow",ylab="Log_2 (GFP1 intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(after_normalization_dir, "control_GFP1_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(data_table_scale_after$M[pos.gfp1,1] ~ ((pos.gfp1 %/% 2916)),
	main="GFP1 per metarow",
	xlab="Metarow",ylab="Log_2 (GFP1 intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();


#signal success and exit.


q(status=0);