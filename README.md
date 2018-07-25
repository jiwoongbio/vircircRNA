# vircircRNA

De novo circular RNA detection from circular viral genome


## Requirements

1. Perl - https://www.perl.org
2. BWA - http://bio-bwa.sourceforge.net
3. Common linux commands: bash, awk, ...


## Install

If you already have Git (https://git-scm.com) installed, you can get the latest development version using Git.
```
git clone https://github.com/jiwoongbio/vircircRNA.git
```

## Usages

```
vircircRNA.sh <reference.fasta> <output.prefix> <threads> <input.1.fastq> [input.2.fastq]
```


## Usages (perl)

1. Generate concatenated sequence of circular chromosome
   ```
   perl vircircRNA_genome.pl -c <circular_chromosome> <reference.fasta> > <reference.concatenated_sequence_of_circular_chromosome.fasta>
   ```

2. Mapping reads to concatenated sequence of circular chromosome
   ```
   bwa index <reference.concatenated_sequence_of_circular_chromosome.fasta>
   bwa mem -T 19 -Y <reference.concatenated_sequence_of_circular_chromosome.fasta> <input.1.fastq> [input.2.fastq] > <mapped_read.sam>
   ```

3. Back-splice junction identification
   ```
   perl vircircRNA_junction.pl -c <circular_chromosome> -A <junction.alignment.html> <mapped_read.sam> <reference.concatenated_sequence_of_circular_chromosome.fasta> > <junction.txt>
   ```
