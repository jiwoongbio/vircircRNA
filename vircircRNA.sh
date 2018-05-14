# Author: Jiwoong Kim (jiwoongbio@gmail.com)
#!/bin/bash

directory=`dirname $0`
referenceFastaFile=$1
outputPrefix=$2
threads=$3
fastqFiles=${@:4}

if [ -z "$fastqFiles" ]; then
	echo 'Usage: ./vircircRNA.sh <reference.fasta> <output.prefix> <threads> <input.1.fastq> [input.2.fastq]' 1>&2
	echo 1>&2
	exit 1
fi
perl -MBio::DB::Fasta -e '' || exit 1

perl $directory/vircircRNA_genome.pl $referenceFastaFile > $outputPrefix.vircircRNA_genome.fasta

bwa index $outputPrefix.vircircRNA_genome.fasta
bwa mem -t $threads -T 19 -Y $outputPrefix.vircircRNA_genome.fasta $fastqFiles | awk -F'\t' '($o !~ /^@/ && and($2, 4) == 0)' > $outputPrefix.vircircRNA_genome.sam

perl $directory/vircircRNA_junction.pl -c '*' -A $outputPrefix.vircircRNA_junction.alignment.html $outputPrefix.vircircRNA_genome.sam $outputPrefix.vircircRNA_genome.fasta > $outputPrefix.vircircRNA_junction.txt
