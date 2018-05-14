# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

my @circularChromosomeList = ();
GetOptions(
	'h' => \(my $help = ''),
	'c=s' => \@circularChromosomeList,
	'l=i' => \(my $sequenceLineLength = 70),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl vircircRNA_genome.pl [options] reference.fasta circular.chromocome [...]

Options: -h       display this help message
         -c STR   circular chromosome
         -l INT   sequence line length [$sequenceLineLength]

EOF
}

my ($referenceFastaFile) = @ARGV;

my @chromosomeList = ();
my %chromosomeSequenceHash = ();
{
	my $chromosome = '';
	open(my $reader, $referenceFastaFile);
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(\S*)/) {
			push(@chromosomeList, $chromosome = $1);
		} else {
			$chromosomeSequenceHash{$chromosome} .= $line;
		}
	}
	close($reader);
}

foreach my $chromosome (@circularChromosomeList ? @circularChromosomeList : @chromosomeList) {
	$chromosomeSequenceHash{$chromosome} = $chromosomeSequenceHash{$chromosome} x 2;
}

foreach my $chromosome (@chromosomeList) {
	my $sequence = $chromosomeSequenceHash{$chromosome};
	print ">$chromosome\n";
	for(my $index = 0; $index < length($sequence); $index += $sequenceLineLength) {
		print substr($sequence, $index, $sequenceLineLength), "\n";
	}
}
