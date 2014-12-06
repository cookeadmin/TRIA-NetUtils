#! /usr/bin/Rscript
library("ggplot2");

infile = "/home/cookeadmin/workspace/all.fasta.gbs_adaptor_blastn.tsv";
data_frame = read.delim(infile, header = TRUE, sep = "\t");
# data_frame;

min_align_length = min(data_frame$align_length);
# min_align_length;
max_align_length = max(data_frame$align_length);
# max_align_length;

# break_points = seq(min(data_frame$align_length), max(data_frame$align_length), by=((max(data_frame$align_length) - min(data_frame$align_length))/(length(data_frame$align_length) - 1)));
# break_points;
png(file="/home/cookeadmin/workspace/test_input_r.png");
hist(data_frame$align_length, col="red", right=FALSE, main="Histogram of GBS Common Adaptor Length", xlab="GBS Common Adaptor Length (bp)", ylab="Frequency (# Observed)");
# histogram_plot <- ggplot(data_frame, aes(x=align_length));
# histogram_plot + geom_histogram(binwidth=1, colour="black", fill="red") + ggtitle("Histogram of GBS Common Adaptor Length") + xlab("GBS Common Adaptor Length (bp)") + ylab("Frequency (# Observed)");

dev.off();
