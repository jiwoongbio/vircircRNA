# Author: Jiwoong Kim (jiwoongbio@gmail.com)
#!/bin/bash

directory=`dirname $0`
fastaFile=$1
gffFile=$2
outputPrefix=$3
fastqFiles=${@:4}

if [ -z "$fastqFiles" ]; then
	echo 'Usage: ./vircircRNA.sh <chromosome.fasta> <gene.gff> <output.prefix> <input.1.fastq> [input.2.fastq]' 1>&2
	echo 1>&2
	exit 1
fi

perl $directory/vircircRNA_chromosome.pl $fastaFile > $outputPrefix.vircircRNA_chromosome.fasta

bwa index $outputPrefix.vircircRNA_chromosome.fasta
if [ -x "$(command -v samtools)" ]; then
	bwa mem -T 19 -Y $outputPrefix.vircircRNA_chromosome.fasta $fastqFiles | samtools view -S -h -F 4 - > $outputPrefix.vircircRNA_chromosome.sam
else
	bwa mem -T 19 -Y $outputPrefix.vircircRNA_chromosome.fasta $fastqFiles > $outputPrefix.vircircRNA_chromosome.sam
fi

perl $directory/vircircRNA_junction.pl -g $gffFile -A $outputPrefix.vircircRNA_junction.alignment.html $outputPrefix.vircircRNA_chromosome.sam $outputPrefix.vircircRNA_chromosome.fasta > $outputPrefix.vircircRNA_junction.txt

if perl -MGD -MGD::Text::Align -e '' 2> /dev/null; then
	perl $directory/vircircRNA_diagram.pl -g $gffFile $outputPrefix.vircircRNA_junction.txt $outputPrefix.vircircRNA_chromosome.fasta > $outputPrefix.vircircRNA_diagram.png
fi
