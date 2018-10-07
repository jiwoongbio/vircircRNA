# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use Getopt::Long qw(:config no_ignore_case);

my $sequenceLineLength = 70;
GetOptions(
	'h' => \(my $help = ''),
	'l' => \(my $isLinearChromosome = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl vircircRNA_chromosome.pl [options] chromosome.fasta [...]

Options: -h       display this help message
         -l       is linear chromosome

EOF
}

my ($fastaFile) = @ARGV;

my @chromosomeList = ();
my %chromosomeSequenceHash = ();
{
	my $chromosome = '';
	open(my $reader, $fastaFile);
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

foreach my $chromosome (@chromosomeList) {
	my $sequence = $chromosomeSequenceHash{$chromosome};
	$sequence = $sequence x 2 if($isLinearChromosome eq '');
	print ">$chromosome\n";
	for(my $index = 0; $index < length($sequence); $index += $sequenceLineLength) {
		print substr($sequence, $index, $sequenceLineLength), "\n";
	}
}
