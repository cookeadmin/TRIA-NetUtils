#! /usr/bin/Rscript 
library('getopt');
library('marray');
library('limma');
library('convert');
library('nnNorm');
library('genefilter');
options(warn=1);

Sys.setlocale(locale="C");

# ./genepix_two_channel_analysis.R -i ~/workspace/walid/32kspruce_earlylatestage_kevin/rawFiles -t ~/workspace/walid/32kspruce_earlylatestage_kevin/in_Files/targets.txt -l ~/workspace/walid/32kspruce_earlylatestage_kevin/galFiles/GQSR02-S1-0001-genepix.gal -o ~/workspace/walid/32kspruce_earlylatestage_kevin
# ./genepix_two_channel_analysis.R -i ~/workspace/walid/extreme_spruce_trees/INPUT_FILES -t ~/workspace/walid/extreme_spruce_trees/targets_Feb-19-2014.txt -l ~/workspace/walid/extreme_spruce_trees/GQSR02-S1-0001-genepix.gal -o ~/workspace/walid/extreme_spruce_trees/OUTPUT_FILES

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

# Make pre, old, and new normalization output directories
pre_normalized_dir <- paste(outfile_dir, "PRE_NORMALIZATION_DIR", sep="/");
dir.create(pre_normalized_dir, showWarnings = FALSE);

old_normalized_dir <- paste(outfile_dir, "OLD_NORMALIZATION_DIR", sep="/");
dir.create(old_normalized_dir, showWarnings = FALSE);

new_normalized_dir <- paste(outfile_dir, "NEW_NORMALIZATION_DIR", sep="/");
dir.create(new_normalized_dir, showWarnings = FALSE);

# Process input data files using BioConductor
target_dataframe  <- readTargets(target_infile);
num_target_files <- length(target_dataframe$FileName);

cat(paste("Processing", num_target_files, "target files for QA/QC Microarray Analysis...\n"));
# RGList <- read.maimages(target_dataframe$FileName, path=infile_dir, "genepix", wt.fun=wtflags(0.1));
# RG = read.maimages("slide 149.gpr", source="genepix", columns=list(R="F633 Mean", G="F543 Mean", Rb="B633 Median", Gb="B543 Median"))
#  <- read.maimages(files,columns=list(R="F635 Mean",G="F532 Mean",Rb="B635 Median",Gb="B532 Median"),annotation=c("Block","Row","Column","ID","Name");
# RGList = read.maimages(target_dataframe$FileName, path=infile_dir, "genepix", columns=list(R="F635 Mean",G="F532 Mean",Rb="B635 Median",Gb="B532 Median"), annotation=c("Block","Row","Column","ID","Name"), sep="\t", wt.fun=wtflags(0.1));
RGList = read.maimages(target_dataframe$FileName, path=infile_dir, "genepix", columns=list(R="F635 Mean",G="F532 Mean",Rb="B635 Median",Gb="B532 Median"), annotation=c("Block","Row","Column","ID","Name"), sep="\t", wt.fun=wtflags(0));

# RGList <- read.maimages(target_dataframe$FileName, path=infile_dir, "genepix", wt.fun=wtflags(0.1));
# R="F650 Mean"
# G="F550 Mean"
# Rb="B650 Mean"
# Gb="B550 Mean"
RGList$genes <- readGAL(galfile=gal_infile, header=TRUE, skip=53, sep="\t", quote="\"");
RGList$printer <- getLayout(RGList$genes);
# show(RGList$genes);
# Obtain file names
target_filenames <- read.delim(target_infile);
colnames(RGList) <- target_filenames[1:num_target_files,1];

# show(RGList);



# Get RGList Red and Green forefround and background intensities.
write.table(RGList$R,
	file = paste(outfile_dir,"RGList_R.txt", sep="/"),
	sep = "\t"
);

write.table(RGList$G,
	file = paste(outfile_dir,"RGList_G.txt", sep="/"),
	sep = "\t"
);

write.table(RGList$Rb,
	file = paste(outfile_dir,"RGList_Rb.txt", sep="/"),
	sep = "\t"
);


write.table(RGList$Gb,
	file = paste(outfile_dir,"RGList_Gb.txt", sep="/"),
	sep = "\t"
);

write.table(RGList$genes,
	file = paste(outfile_dir,"RGList_genes.txt", sep="/"),
	sep = "\t"
);

write.table(RGList$printer$status,
	file = paste(outfile_dir,"RGList_printer.txt", sep="/"),
	sep = "\t"
);

write.table(RGList$targets,
	file = paste(outfile_dir,"RGList_targets.txt", sep="/"),
	sep = "\t"
);

write.table(RGList$weights,
	file = paste(outfile_dir,"RGList_weights.txt", sep="/"),
	sep = "\t"
);

# NORMALIZATION METHODS
# Before normalization
# Correct for infinite values.
# RGList.b$R[is.infinite(RGList.b$R)] <- 0.1;
# RGList.b$R[RGList.b$R <= 0] <- 0.1;
# RGList.b$G[is.infinite(RGList.b$G)] <- 0.1;
# RGList.b$G[RGList.b$G <= 0] <- 0.1;
# RGList.b$Rb[is.infinite(RGList.b$Rb)] <- 0.1;
# RGList.b$Rb[RGList.b$Rb <= 0] <- 0.1;
# RGList.b$Rg[is.infinite(RGList.b$Gb)] <- 0.1;
# RGList.b$Rg[RGList.b$Rg <= 0] <- 0.1;

RGList.b <- backgroundCorrect(RGList, method="normexp");
# MAList <- normalizeWithinArrays(RGList, method="none", bc.method="none");
MAList <- normalizeWithinArrays(RGList.b, method="none", bc.method="none");
raw_data_table <- as(RGList, "marrayRaw");

# Apply old normalization method: ( withinarray = print-tip loess; between array = quantile norm)
data_table_old_norm <- normalizeWithinArrays(RGList.b);
data_table_scale_old <- normalizeBetweenArrays(data_table_old_norm, method="quantile");
data_table_old_norm_marray <- as(data_table_scale_old, "marrayNorm");

# Apply new normalization method (withinarray = 2D spatial and print-tip loess; between array = quantile norm)
data_table_scale_new <- normalizeBetweenArrays(RGList, method="quantile");
data_table_new_marray <- as(data_table_scale_new, "marrayNorm");

# Write normalized data to file
old_norm_data <- cbind(data_table_scale_old$genes$ID, data_table_scale_old$M);
new_norm_data <- cbind(data_table_scale_new$genes$ID, data_table_scale_new$M);
colnames(old_norm_data)[1] <- "ID";
colnames(new_norm_data)[1] <- "ID";

write.table(old_norm_data, file=paste(old_normalized_dir, "normalized_data.txt", sep="/"));
write.table(new_norm_data, file=paste(new_normalized_dir, "normalized_data.txt", sep="/"));

# DIAGNOSTIC PLOTS
# 1. boxplot for all slides
png(file=paste(pre_normalized_dir, "boxplot_raw.png", sep="/"));
par(mar=c(7,5,1,1));
maBoxplot(raw_data_table, ylab="maM", main="Boxplots on M-value by array before normalization", xlab="", las=2);
abline(h=0, col="grey", lwd=2);
dev.off();

png(file = paste(old_normalized_dir, "boxplot_within_slide.png", sep="/"));
par(mar=c(7,5,1,1));
maBoxplot(as(data_table_old_norm,"marrayNorm"), ylab="maM", main="Boxplots on M-value by array after within-array normalization", xlab="", las=2);
abline(h=0, col="grey", lwd=2);
dev.off();

png(file = paste(old_normalized_dir, "boxplot_between_slides.png", sep="/"));
par(mar=c(7,5,1,1));
maBoxplot(data_table_old_norm_marray, ylab="maM", main="Boxplots on M-value by array after between-array normalization", xlab="", las=2);
abline(h=0, col="grey", lwd=2);
dev.off();

png(file = paste(new_normalized_dir, "boxplot_between_slides.png", sep="/"));
par(mar=c(7,5,1,1));
maBoxplot(data_table_new_marray, ylab="maM", main="Boxplots on M-value by array after between-arrays normalization", xlab="", las=2);
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
	png(file = paste(pre_normalized_dir, "/", paste("boxplots", index, colnames(RGList)[slides], sep="_"), ".png", sep=""))
	maBoxplot(raw_data_table[,slides], xlab="print-tip group", ylab="M (log2-scale)", main=paste("Boxplots on M-value by print-tip group before normalization for slide", colnames(RGList)[slides], sep=" "));
	abline(h=0, col="grey", lwd=2);
	dev.off();
	
	png(file = paste(old_normalized_dir,"/", paste("boxplots", index, colnames(RGList)[slides], sep="_"), ".png", sep=""));
	maBoxplot(data_table_old_norm_marray[,slides], xlab="print-tip group", ylab="M (log2-scale)", main= paste("Boxplots on M-value by print-tip group after normalization for slide", colnames(RGList)[slides], sep=" "));
	abline(h=0, col="grey", lwd=2);
	dev.off();

	png(file = paste(new_normalized_dir, "/", paste("boxplots", index, colnames(RGList)[slides], sep="_"), ".png", sep=""));
	maBoxplot(data_table_new_marray[,slides], xlab="print-tip group", ylab="M (log2-scale)", main= paste("Boxplots on M-value by print-tip group after normalization for slide", colnames(RGList)[slides], sep=" "));
	abline(h=0, col="grey", lwd=2);
	dev.off();
	
	#2. 2D spatial plots
	png(file = paste(pre_normalized_dir, "/", paste("spatialplots_Rb", index, colnames(RGList)[slides], sep="_"), ".png", sep=""));
	image(raw_data_table[,slides], xvar="maRb", bar=TRUE, main=paste("Spatial plots on background intensity of red channel for slide", colnames(RGList)[slides], sep=" "));
	dev.off();

	png(file = paste(pre_normalized_dir, "/", paste("spatialplots_Gb", index, colnames(RGList)[slides], sep="_"), ".png", sep=""));
	image(raw_data_table[,slides], xvar="maGb", bar=TRUE, main=paste("Spatial plots on background intensity of green channel for slide", colnames(RGList)[slides], sep=" "));
	dev.off();

	#3. MA plots
	png(file = paste(pre_normalized_dir, "/", paste("MAplots", index, colnames(RGList)[slides], sep="_"), ".png", sep=""));
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
	
	png(file = paste(old_normalized_dir, "/", paste("MAplots", index, colnames(RGList)[slides], sep="_"), ".png", sep=""));
	limma::plotMA(data_table_scale_old,
		array = slides,
		xlab = "A:Intensity (log2)",
		ylab = "M:Fold Change (log2)",
		main = paste("MAplot after normalization for slide", colnames(RGList)[slides], sep=" "),
		xlim = c(2,16),
		ylim = c(-6,6)
	);
	abline (h=0, col="red");
	dev.off();


	png(file = paste(new_normalized_dir, "/", paste("MAplots", index, colnames(RGList)[slides], sep="_"), ".png", sep=""));
	limma::plotMA(data_table_scale_new,
		array = slides,
		xlab = "A:Intensity (log2)",
		ylab = "M:Fold Change (log2)",
		main = paste("MAplot after normalization for slide", colnames(RGList)[slides], sep=" "),
		xlim = c(2,16),
		ylim = c(-6,6)
	);
	abline (h=0, col="red");
	dev.off();    
}  

#more density plot
png(file=paste(pre_normalized_dir, "densityplots_raw_overall.png", sep="/"));
plotDensities(RGList.b, log=TRUE, arrays=c(1:6), groups=c(rep(1,6), rep(2,6)), col=c("red","green"));
dev.off();

summary(data_table_old_norm);
png(file=paste(old_normalized_dir, "densityplots_within_array_overall.png", sep="/"));
plotDensities(data_table_old_norm, log=TRUE, arrays=c(1:6), groups=c(rep(1,6), rep(2,6)), col=c("red","green"));
dev.off();

png(file=paste(old_normalized_dir, "densityplots_between_array_overall.png", sep="/"));
plotDensities(data_table_scale_old, log=TRUE, arrays=c(1:6), groups=c(rep(1,6), rep(2,6)), col=c("red","green"));
dev.off();

png(file=paste(new_normalized_dir, "densityplots_between_array_overall.png", sep="/"));
plotDensities(data_table_scale_new, log=TRUE, arrays=c(1:6), groups=c(rep(1,6), rep(2,6)), col=c("red","green"));
dev.off();

# QUALITY ACCESSMENT ON CONTROLS
pos.buffer <- grep("Buffer",RGList$gene[,4]);
pos.empty <- grep("Empty", RGList$gene[,4]);
pos.gfp <- grep("GFP", RGList$gene[,4]);
pos.gfp1 <- grep("GFP1", RGList$gene[,4]);

# The array has 12 metarows. Each metarow contains 2916 probes
# Boxplot: variation of the control spots intensities per metarow
png(file=paste(pre_normalized_dir, "control_buffer_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(MAList$M[pos.buffer,1] ~ ((pos.buffer %/% 2916)),
	main="Buffer per metarow",
	xlab="Metarow", ylab="Log2(Buffer intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(old_normalized_dir, "control_buffer_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(data_table_scale_old$M[pos.buffer,1] ~ ((pos.buffer %/% 2916)),
	main="Buffer per metarow",
	xlab="Metarow", ylab="Log2(Buffer intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(new_normalized_dir, "control_buffer_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(data_table_scale_new$M[pos.buffer,1] ~ ((pos.buffer %/% 2916)),
	main="Buffer per metarow",
	xlab="Metarow",ylab="Log2(Buffer intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(pre_normalized_dir, "control_empty_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(MAList$M[pos.empty,1] ~ ((pos.empty %/% 2916)),
	main="Empty per metarow",
	xlab="Metarow",ylab="Log_2 (Empty intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(old_normalized_dir, "control_empty_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(data_table_scale_old$M[pos.empty,1] ~ ((pos.empty %/% 2916)),
	main="Empty per metarow",
	xlab="Metarow",ylab="Log_2 (Empty intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(new_normalized_dir, "control_empty_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(data_table_scale_new$M[pos.empty,1] ~ ((pos.empty %/% 2916)),
	main="Empty per metarow",
	xlab="Metarow",ylab="Log_2 (Empty intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(pre_normalized_dir, "control_GFP_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(MAList$M[pos.gfp,1] ~ ((pos.gfp %/% 2916)),
	main="GFP per metarow",
	xlab="Metarow",ylab="Log_2 (GFP intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(old_normalized_dir, "control_GFP_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(data_table_scale_old$M[pos.gfp,1] ~ ((pos.gfp %/% 2916)),
	main="GFP per metarow",
	xlab="Metarow",ylab="Log_2 (GFP intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(new_normalized_dir, "control_GFP_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(data_table_scale_new$M[pos.gfp,1] ~ ((pos.gfp %/% 2916)),
	main="GFP per metarow",
	xlab="Metarow",ylab="Log_2 (GFP intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(pre_normalized_dir, "control_GFP1_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(MAList$M[pos.gfp1,1] ~ ((pos.gfp1 %/% 2916)),
	main="GFP1 per metarow",
	xlab="Metarow",ylab="Log_2 (GFP1 intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(old_normalized_dir, "control_GFP1_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(data_table_scale_old$M[pos.gfp1,1] ~ ((pos.gfp1 %/% 2916)),
	main="GFP1 per metarow",
	xlab="Metarow",ylab="Log_2 (GFP1 intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

png(file=paste(new_normalized_dir, "control_GFP1_by_metarow.png", sep="/"), width=800, height=600, res=NA, units="px");
boxplot(data_table_scale_new$M[pos.gfp1,1] ~ ((pos.gfp1 %/% 2916)),
	main="GFP1 per metarow",
	xlab="Metarow",ylab="Log_2 (GFP1 intensity)",
	col=c("RED","GREEN","BLUE","CYAN","YELLOW","PURPLE")
);
dev.off();

# ANNOTATION 
# Build annotation object
# annotated_data <- read.csv(paste(infile_dir, "annotated_spruce3.tab", sep="/"),
# 	header = TRUE,
# 	sep="\t",
# 	fill = TRUE
# );
# 
# annotated_data <- unique(annotated_data);
# rownames(annotated_data) <- annotated_data$Probe;
# 
# # STATISTICAL ANALYSIS FOR SPRUCES
# pos.spruce <- grep("Spruce_",RGList$gene[,4]);
# myDesign <- cbind(Time1 = c(1,1,1,1,1,1,0,0,0,0,0), Time2 = c(0,0,0,0,0,0,1,1,1,1,1));
# print(myDesign);
# 
# # Non-specific filter 1:: remove all the control spots
# non_specific_filter_old <- data_table_scale_old[pos.spruce,];
# non_specific_filter_new <- data_table_scale_new[pos.spruce,];
# 
# # Fit into a linear model
# fit_old <- lmFit(non_specific_filter_old, design=myDesign, ndups=1,weights=non_specific_filter_old$weights);
# fit_new <- lmFit(non_specific_filter_new, design=myDesign, ndups=1,weights=non_specific_filter_old$weights);
# 
# # Non-specific filter 2:: remove those with low intensity value (ie Amean < 7)
# non_specific_filter_old_2 <- fit_old[fit_old$Amean >= 7,];
# non_specific_filter_new_2 <- fit_new[fit_new$Amean >= 7,];
# 
# # Apply empirical bayes method as moderated t-statistic
# ebayes_old <- eBayes(non_specific_filter_old_2, proportion=0.01);
# ebayes_new <- eBayes(non_specific_filter_new_2, proportion=0.01);
# 
# # Without filter 2
# ebayes_old_full <- eBayes(fit_old, proportion = 0.01);
# ebayes_new_full <- eBayes(fit_new, proportion = 0.01);
# 
# 
# # ## fdr adjust for multiple comparisons
# n_comparison <- ncol(myDesign);
# 
# for (i in 1:n_comparison){
# 	comparison <- colnames(myDesign)[i];
# 	top_non_specific_filter_old <- topTable(ebayes_old,   			
# 		coef = i,
# 		number = nrow(ebayes_old),
# 		adjust.method="BH",
# 		sort.by = "P"
# 	);
# 	top_all_old <- topTable(ebayes_old_full,
# 		coef = i,
# 		number = nrow(ebayes_old_full),
# 		adjust.method="BH",
# 		sort.by = "P"
# 	);
# 
# 	up_001_old <- nrow(top_non_specific_filter_old[top_non_specific_filter_old$adj.P.Val < 0.01 & top_non_specific_filter_old$logFC > 0,]);
# 	up_005_old <- nrow(top_non_specific_filter_old[top_non_specific_filter_old$adj.P.Val < 0.05 & top_non_specific_filter_old$logFC > 0,]);
# 	down_001_old <- nrow(top_non_specific_filter_old[top_non_specific_filter_old$adj.P.Val < 0.01 & top_non_specific_filter_old$logFC < 0,]);
# 	down_005_old <- nrow(top_non_specific_filter_old[top_non_specific_filter_old$adj.P.Val < 0.05 & top_non_specific_filter_old$logFC < 0,]);
# 	exclude_gene_names_old <- setdiff(top_all_old$ID, top_non_specific_filter_old$ID);
# 
# 	filter_genelist_old <- subset(top_non_specific_filter_old, select=c(ID, logFC, AveExpr, t, P.Value, adj.P.Val, B));
# 	exclude_genelist_old <- subset(top_all_old[top_all_old$ID %in% exclude_gene_names_old,], select=c(ID, logFC, AveExpr, t, P.Value, adj.P.Val, B));
# 	exclude_genelist_old$t <- "NA";
# 	exclude_genelist_old$P.Value <- "NA";
# 	exclude_genelist_old$adj.P.Val <- "NA";
# 	exclude_genelist_old$B <- "NA";
# 
# 	all_genelist_old <- rbind(filter_genelist_old, exclude_genelist_old);
# 	rownames(all_genelist_old) <- all_genelist_old$ID;
# 	all_genelist_old <- merge(all_genelist_old, annotated_data, by="row.names", all.x=TRUE, sort = FALSE);
# 
# 	top_non_specific_filter_new <- topTable(ebayes_new,     		
# 		coef = i,
# 		number = nrow(ebayes_new),
# 		adjust.method="BH",
# 		sort.by = "P"
# 	);
# 	top_all_new <- topTable(ebayes_new_full,
# 		coef = i,
# 		number = nrow(ebayes_new_full),
# 		adjust.method="BH",
# 		sort.by = "P"
# 	);
# 
# 	up_001_new <- nrow(top_non_specific_filter_new[top_non_specific_filter_new$adj.P.Val < 0.01 & top_non_specific_filter_new$logFC > 0,]);
# 	up_005_new <- nrow(top_non_specific_filter_new[top_non_specific_filter_new$adj.P.Val < 0.05 & top_non_specific_filter_new$logFC > 0,]);
# 	down_001_new <- nrow(top_non_specific_filter_new[top_non_specific_filter_new$adj.P.Val < 0.01 & top_non_specific_filter_new$logFC < 0,]);
# 	down_005_new <- nrow(top_non_specific_filter_new[top_non_specific_filter_new$adj.P.Val < 0.05 & top_non_specific_filter_new$logFC < 0,]);
# 	exclude_gene_names_new <- setdiff(top_all_new$ID, top_non_specific_filter_new$ID);
# 
# 	filter_genelist_new <- subset(top_non_specific_filter_new, select=c(ID, logFC, AveExpr, t, P.Value, adj.P.Val, B));
# 	exclude_genelist_new <- subset(top_all_new[top_all_new$ID %in% exclude_gene_names_new,], select=c(ID, logFC, AveExpr, t, P.Value, adj.P.Val, B));
# 	exclude_genelist_new$t <- "NA";
# 	exclude_genelist_new$P.Value <- "NA";
# 	exclude_genelist_new$adj.P.Val <- "NA";
# 	exclude_genelist_new$B <- "NA";
# 
# 	all_genelist_new <- rbind(filter_genelist_new, exclude_genelist_new);
# 	rownames(all_genelist_new) <- all_genelist_new$ID;
# 	all_genelist_new <- merge(all_genelist_new, annotated_data, by="row.names", all.x=TRUE, sort = FALSE);
# 
# 
# 	## significant genes count
# 	count <- rbind(c(up_001_old, down_001_old, up_005_old, down_005_old), c(up_001_new, down_001_new, up_005_new, down_005_new));
# 	colnames(count) <- c("0.01 UP","0.01 DOWN","0.05 UP", "0.05 DOWN");
# 	rownames(count) <- c("Old Normalization", "New Normalization");
# # 	count
# 
# 	## output result
# 	write.table(all_genelist_old,
# 		file = paste(old_outPath,"/",comparison,"_genelist.xls", sep=""),
# 		sep = "\t"
# 	);
# 
# 	write.table(all_genelist_new,
# 		file = paste(new_outPath,"/",comparison,"_genelist.xls", sep=""),
# 		sep = "\t"
# 	);
# 
# 	write.table(count,
# 		file = paste(pre_outPath,"/",comparison,"_genecount.txt",sep=""),
# 		sep="\t",
# 		row.names = FALSE,
# 		col.names = FALSE
# 	);
# 
# 
# 	######## VENN DIAGRAM TO COMPARE TWO METHODS ##########################
# 	universe <- annotated_data$Probe;
# 	old_005 <- filter_genelist_old[filter_genelist_old$adj.P.Val < 0.05, 1];
# 	new_005 <- filter_genelist_new[filter_genelist_new$adj.P.Val < 0.05, 1];
# 
# 	old_001 <- filter_genelist_old[filter_genelist_old$adj.P.Val < 0.01, 1];
# 	new_001 <- filter_genelist_new[filter_genelist_new$adj.P.Val < 0.01, 1];
# 
# 	# the count matrix
# 	count_005 <- matrix(0, nrow= length(universe), ncol=2);
# 	colnames(count_005) <- c("OLD","NEW");
# 
# 	for (i in 1:length(universe)){
# 		count_005[i,1] <- universe[i] %in% old_005;
# 		count_005[i,2] <- universe[i] %in% new_005;
# 	}
# 	
# 	venn_005 <- vennCounts(count_005)
# 	png( file = paste(pre_outPath,"/",comparison,"_venn_diagram_005.png", sep=""), height=600, width=800);
# 	vennDiagram(venn_005,
# 		cex = 1, lwd = 2,
# 		counts.col = "red", circle.col = c("Blue", "Orange")
# 	);
# 	dev.off();
# 
# 	count_001 <- matrix(0, nrow= length(universe), ncol=2);
# 	colnames(count_001) <- c("OLD","NEW");
# 
# 	for (i in 1:length(universe)){
# 		count_001[i,1] <- universe[i] %in% old_001;
# 		count_001[i,2] <- universe[i] %in% new_001;
# 	}
# 	
# 	venn_001 <- vennCounts(count_001);
# 	png( file = paste(pre_outPath,"/",comparison,"_venn_diagram_001.png", sep=""), height=600, width=800);
# 	vennDiagram(venn_001,
# 		cex = 1, lwd = 2,
# 		counts.col = "red", circle.col = c("Blue", "Orange")
# 	);
# 	dev.off();
# }

#signal success and exit.
q(status=0);