#!/bin/sh

STAR=/usr/local/bin/staticSTAR

REFERENCE_GENOME_INPUT_DIR=$1

echo "Using reference genome input directory: $REFERENCE_GENOME_INPUT_DIR\n"

NUM_FILES=$2

FILENAME=$3

echo "Using reference genome filename: $FILENAME\n"

CDNA_FASTQ_FILE=$4

echo "Using cDNA fastq input file: $CDNA_FASTQ_FILE\n"

NUM_CPU_THREADS=$5
echo "Using $NUM_CPU_THREADS CPU cores\n"

STAR_OUTPUT_DIR=$6

echo "Using STAR output directory: $STAR_OUTPUT_DIR\n"
mkdir -p $STAR_OUTPUT_DIR

for i in $(seq 1 $NUM_FILES)
do

REFERENCE_FASTA=$REFERENCE_GENOME_INPUT_DIR/$FILENAME-scaffolds.part$i.fa

echo "Executing $STAR using $REFERENCE_FASTA\n";


# 1) STAR uses genome index files that must be saved in unique directories.

genomeBaseDir=$STAR_OUTPUT_DIR/$FILENAME-part$i-dir

echo "Creating the $genomeBaseDir directory.\n"
mkdir -p $genomeBaseDir

genomeDir1=$genomeBaseDir/$FILENAME-part$i-refgen-pass1-dir

echo "Creating the $genomeDir1 directory.\n"
mkdir -p $genomeDir1

echo $STAR --runMode genomeGenerate --genomeDir $genomeDir1 --genomeFastaFiles $REFERENCE_FASTA  --runThreadN $NUM_CPU_THREADS


# 2) Alignment jobs were executed as follows:

runDir1=$genomeBaseDir/$FILENAME-part$i-alignment-pass1-dir

echo "Creating the $runDir1 directory."
mkdir -p $runDir1

echo "Changing directory to the $runDir1 directory.\n"
cd $runDir1

echo $STAR --genomeDir $genomeDir1 --readFilesIn $CDNA_FASTQ_FILE --outSAMmode Full --outSAMattributes All --outFileNamePrefix $FILENAME-part$i-alignment-pass1 --runThreadN $NUM_CPU_THREADS

# 3) For the 2-pass STAR, a new index is then created using splice junction information contained in the file SJ.out.tab from the first pass:

genomeDir2=$genomeBaseDir/$FILENAME-part$i-refgen-pass2-dir

echo "Creating the $genomeDir2 directory.\n"
mkdir -p $genomeDir2

echo $STAR --runMode genomeGenerate --genomeDir $genomeDir2 --genomeFastaFiles $REFERENCE_FASTA --sjdbFileChrStartEnd $runDir1/$FILENAME-part$i-alignment-pass2SJ.out.tab --sjdbOverhang 75 --runThreadN $NUM_CPU_THREADS


# 4) The resulting index is then used to produce the final alignments as follows:

runDir2=/$genomeBaseDir/$FILENAME-part$i-alignment-pass2-dir

echo "Creating the $runDir2 directory."
mkdir -p $runDir2

echo "Changing directory to the $runDir2 directory.\n"
cd $runDir2

echo $STAR --genomeDir $genomeDir2 --readFilesIn $CDNA_FASTQ_FILE --outSAMmode Full --outSAMattributes All --outFileNamePrefix $FILENAME-part$i-alignment-pass2 --runThreadN $NUM_CPU_THREADS

done
