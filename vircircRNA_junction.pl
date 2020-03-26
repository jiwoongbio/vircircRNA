# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use List::Util qw(min max sum any all);
use Scalar::Util qw(looks_like_number);
use Getopt::Long qw(:config no_ignore_case);

my @blockRegionList = ();
GetOptions(
	'h' => \(my $help = ''),
	'l' => \(my $isLinearChromosome = ''),
	'p' => \(my $printAll = ''),
	'q=i' => \(my $minimumMappingQuality = 0),
	'm=i' => \(my $minimumAlignmentLength = 16),
	's=s' => \(my $stranded = ''),
	'b=s' => \@blockRegionList,
	'g=s' => \(my $gffFile = ''),
	'f=s' => \(my $feature = 'gene'),
	'a=s' => \(my $attribute = 'gene'),
	'A=s' => \(my $alignmentFile = ''),
	'L=i' => \(my $alignmentLength = 100),
	'F=i' => \(my $fontSize = 10),
	'C=s' => \(my $colors = '#880000,#000088'),
	'R=s' => \(my $readNameFile = ''),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl vircircRNA_junction.pl [options] softclipping.sam chromosome.fasta > junction.txt

Options: -h       display this help message
         -l       is linear chromosome
         -p       print both forward-splice and back-splice junctions
         -q INT   minimum mapping quality [$minimumMappingQuality]
         -m INT   minimum alignment length [$minimumAlignmentLength]
         -s STR   stranded, "f" or "r"
         -b STR   block region
         -g STR   GTF file
         -f STR   GTF feature [$feature]
         -a STR   GTF attribute [$attribute]
         -A FILE  alignment HTML file
         -L INT   alignment length [$alignmentLength]
         -F INT   alignment font size [$fontSize]
         -C STR   alignment colors [$colors]
         -R FILE  read name file

EOF
}
@blockRegionList = map {($_ =~ /^(.+):([0-9]+)-([0-9]+)$/) ? [$1, $2, $3] : ()} @blockRegionList;
my ($color1, $color2) = split(/,/, $colors);

my ($samFile, $fastaFile) = @ARGV;

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
my %chromosomeLengthHash = map {$_ => length($chromosomeSequenceHash{$_})} @chromosomeList;

if($isLinearChromosome eq '') {
	$chromosomeLengthHash{$_} = $chromosomeLengthHash{$_} / 2 foreach(@chromosomeList);
	push(@blockRegionList, map {[$_->[0], $_->[1] + $chromosomeLengthHash{$_->[0]}, $_->[2] + $chromosomeLengthHash{$_->[0]}]} @blockRegionList);
}

my %junctionReadCountHash = ();
my %junctionAlignmentListHash = ();
my %junctionReadNameListHash = ();
open(my $reader, $samFile);
my ($readName, @tokenHashList) = ('');
while(my $line = <$reader>) {
	chomp($line);
	next if($line =~ /^@/);
	my %tokenHash = ();
	(@tokenHash{'qname', 'flag', 'rname', 'pos', 'mapq', 'cigar', 'rnext', 'pnext', 'tlen', 'seq', 'qual'}, my @tagTypeValueList) = split(/\t/, $line);
	$tokenHash{"$_->[0]:$_->[1]"} = $_->[2] foreach(map {[split(/:/, $_, 3)]} @tagTypeValueList);
	next if($tokenHash{'flag'} & 4);
	if($tokenHash{'qname'} ne $readName) {
		addJunctionList() if(scalar(@tokenHashList) > 0);
		($readName, @tokenHashList) = ($tokenHash{'qname'});
	}
	push(@tokenHashList, \%tokenHash);
}
addJunctionList() if(scalar(@tokenHashList) > 0);
close($reader);

my @junctionCountList = ();
my %chromosomePositionCountHash = ();
foreach(sort {compare($a, $b)} map {[split(/\t/, $_)]} keys %junctionReadCountHash) {
	my ($chromosome, $position1, $position2, $strand) = @$_;
	my $junction = join("\t", $chromosome, $position1, $position2, $strand);
	my $count = sum(values %{$junctionReadCountHash{$junction}});
	if($position2 <= $position1 || $printAll) {
		if(max(map {length($_->[0])} @{$junctionAlignmentListHash{$junction}}) >= $minimumAlignmentLength && max(map {length($_->[1])} @{$junctionAlignmentListHash{$junction}}) >= $minimumAlignmentLength) {
			push(@junctionCountList, [$junction, $chromosome, $position1, $position2, $strand, $count]);
		}
	}
	$chromosomePositionCountHash{$chromosome}->{"+$position1"} += $count;
	$chromosomePositionCountHash{$chromosome}->{"-$position2"} += $count;
}
if($gffFile ne '') {
	my %junctionGeneListHash = ();
	open(my $reader, $gffFile);
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^#/);
		next if($line eq '');
		my %tokenHash = ();
		@tokenHash{'chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'} = split(/\t/, $line);
		my @attributeList = split(/; */, $tokenHash{'attribute'});
		my %attributeHash = ();
		if(all {/^(\S+)=(.*)$/} @attributeList) {
			foreach(@attributeList) {
				/^(\S+)=(.*)$/;
				$attributeHash{$1} = $2;
			}
		} elsif(all {/^(\S+) "(.*)"$/} @attributeList) {
			foreach(@attributeList) {
				/^(\S+) "(.*)"$/;
				$attributeHash{$1} = $2;
			}
		}
		if($tokenHash{'feature'} eq $feature) {
			my @startEndList = ([@tokenHash{'start', 'end'}]);
			push(@startEndList, [map {$_ + $chromosomeLengthHash{$tokenHash{'chromosome'}}} @tokenHash{'start', 'end'}]) if($isLinearChromosome eq '');
			foreach(@junctionCountList) {
				my ($junction, $chromosome, $position1, $position2, $strand, $count) = @$_;
				my ($start, $end) = sort {$a <=> $b} ($position1, $position2);
				if($tokenHash{'chromosome'} eq $chromosome && $tokenHash{'strand'} eq $strand && any {$_->[0] <= $end && $start <= $_->[1]} @startEndList) {
					if(defined(my $gene = $attributeHash{$attribute})) {
						push(@{$junctionGeneListHash{$junction}}, $gene);
					}
				}
			}
		}
	}
	close($reader);

	@junctionCountList = grep {defined($_->[-1])} map {[@$_, $junctionGeneListHash{$_->[0]}]} @junctionCountList;
	foreach(@junctionCountList) {
		my ($junction, $chromosome, $position1, $position2, $strand, $count, $geneList) = @$_;
		my $ratio = ($count + $count) / ($chromosomePositionCountHash{$chromosome}->{"+$position1"} + $chromosomePositionCountHash{$chromosome}->{"-$position2"});
		my %geneHash = ();
		my $gene = join(',', map {defined($geneHash{$_}) ? () : ($geneHash{$_} = $_)} @$geneList);
		print join("\t", $chromosome, $position1, $position2, $strand, $count, $ratio, $gene), "\n";
	}
} else {
	foreach(@junctionCountList) {
		my ($junction, $chromosome, $position1, $position2, $strand, $count) = @$_;
		my $ratio = ($count + $count) / ($chromosomePositionCountHash{$chromosome}->{"+$position1"} + $chromosomePositionCountHash{$chromosome}->{"-$position2"});
		print join("\t", $chromosome, $position1, $position2, $strand, $count, $ratio), "\n";
	}
}
if($alignmentFile ne '' && @junctionCountList) {
	open(my $writer, "> $alignmentFile");
	print $writer <<EOF;
<!DOCTYPE html>
<html>
<head>
<style>
* {
	font-size: ${fontSize}px;
	line-height: ${fontSize}px;
}
</style>
</head>
<body>
EOF
	foreach(@junctionCountList) {
		my ($junction, $chromosome, $position1, $position2, $strand, $count) = @$_;
		if(defined(my $alignmentList = $junctionAlignmentListHash{$junction})) {
			my @alignmentList = sort {length($b->[0]) <=> length($a->[0]) || length($a->[1]) <=> length($b->[1]) || "$a->[0]$a->[1]" cmp "$b->[0]$b->[1]"} @$alignmentList;
			my $header1 = (' ' x ($alignmentLength - length("$chromosome:$position1"))) . "$chromosome:$position1";
			my $header2 = $position2 . (' ' x ($alignmentLength - length($position2)));
			my $sequence1 = substr($chromosomeSequenceHash{$chromosome}, $position1 - $alignmentLength, $alignmentLength);
			my $sequence2 = substr($chromosomeSequenceHash{$chromosome}, $position2 - 1, $alignmentLength);
			foreach(@alignmentList) {
				my ($alignment1, $alignment2) = @$_;
				$alignment1 = substr($alignment1, -$alignmentLength);
				$alignment1 = (' ' x ($alignmentLength - length($alignment1))) . $alignment1;
				$alignment2 = substr($alignment2, 0, $alignmentLength);
				$alignment2 = $alignment2 . (' ' x ($alignmentLength - length($alignment2)));
				@$_ = ($alignment1, $alignment2);
			}
			print $writer "<p><pre>\n";
			foreach(["<b>$header1</b>", "<b>$header2</b>"], ["<b>$sequence1</b>", "<b>$sequence2</b>"], @alignmentList) {
				print $writer join(' ', "<span style=\"color: $color1;\">$_->[0]</span>", "<span style=\"color: $color2;\">$_->[1]</span>"), "\n";
			}
			print $writer "</pre></p>\n";
		}
	}
	print $writer <<EOF;
</body>
</html>
EOF
	close($writer);
}
if($readNameFile ne '' && @junctionCountList) {
	open(my $writer, "> $readNameFile");
	foreach(@junctionCountList) {
		my ($junction, $chromosome, $position1, $position2, $strand, $count) = @$_;
		foreach(@{$junctionReadNameListHash{$junction}}) {
			my ($readName, $number, $readStrand) = @$_;
			print $writer join("\t", $chromosome, $position1, $position2, $strand, $readName, $number, $readStrand), "\n";
		}
	}
	close($writer);
}

sub addJunctionList {
	my %numberTokenHashListHash = ();
	foreach(@tokenHashList) {
		my %tokenHash = %$_;
		my $number = ($tokenHash{'flag'} & 192) / 64;
		push(@{$numberTokenHashListHash{$number}}, \%tokenHash);
	}
	foreach my $number (sort {$a <=> $b} keys %numberTokenHashListHash) {
		my @tokenHashList = @{$numberTokenHashListHash{$number}};
		my %sequenceHash = ();
		foreach(@tokenHashList) {
			my %tokenHash = %$_;
			my $sequence = $tokenHash{'seq'};
			$sequence = getReverseComplementarySequence($sequence) if($tokenHash{'flag'} & 16);
			$sequenceHash{$sequence} = 1;
		}
		next if(scalar(my ($sequence) = keys %sequenceHash) > 1);
		my %indexBreakpointListHash = ();
		foreach(@tokenHashList) {
			my %tokenHash = %$_;
			{
				my ($chromosome, $position, $cigar, $isReverse) = (@tokenHash{'rname', 'pos', 'cigar'}, ($tokenHash{'flag'} & 16));
				my @positionList = getPositionList($position, $cigar);
				my @baseList = split(//, $isReverse ? getReverseComplementarySequence($sequence) : $sequence);
				my %positionBaseHash = map {$_->[0] => $_->[1]} grep {$_->[0] ne ''} map {[$positionList[$_], $baseList[$_]]} (0 .. $#positionList);
				@positionList = reverse(@positionList) if($isReverse);
				if($cigar =~ /^([0-9]+)S/) {
					my $index = $isReverse ? $1 - length($tokenHash{'seq'}) : $1;
					push(@{$indexBreakpointListHash{$index}}, [$chromosome, '-', \@positionList, \%positionBaseHash]);
				}
				if($cigar =~ /([0-9]+)S$/) {
					my $index = $isReverse ? $1 : $1 - length($tokenHash{'seq'});
					push(@{$indexBreakpointListHash{$index}}, [$chromosome, '+', \@positionList, \%positionBaseHash]);
				}
			}
			while(defined($tokenHash{'XA:Z'}) && $tokenHash{'XA:Z'} =~ /([^,;]+),([-+][0-9]+),([^,;]+),([0-9]+);/g) {
				my ($chromosome, $position, $cigar, $isReverse) = ($1, abs($2), $3, $2 < 0);
				my @positionList = getPositionList($position, $cigar);
				my @baseList = split(//, $isReverse ? getReverseComplementarySequence($sequence) : $sequence);
				my %positionBaseHash = map {$_->[0] => $_->[1]} grep {$_->[0] ne ''} map {[$positionList[$_], $baseList[$_]]} (0 .. $#positionList);
				@positionList = reverse(@positionList) if($isReverse);
				if($cigar =~ /^([0-9]+)S/) {
					my $index = $isReverse ? $1 - length($tokenHash{'seq'}) : $1;
					push(@{$indexBreakpointListHash{$index}}, [$chromosome, '-', \@positionList, \%positionBaseHash]);
				}
				if($cigar =~ /([0-9]+)S$/) {
					my $index = $isReverse ? $1 : $1 - length($tokenHash{'seq'});
					push(@{$indexBreakpointListHash{$index}}, [$chromosome, '+', \@positionList, \%positionBaseHash]);
				}
			}
		}
		foreach my $index1 (sort {$a <=> $b} grep {$_ > 0} keys %indexBreakpointListHash) {
			foreach my $index2 (sort {$a <=> $b} grep {$_ < 0} keys %indexBreakpointListHash) {
				if(my @indexList = $index1 .. abs($index2)) {
					foreach(getJunctionList(@indexBreakpointListHash{$index1, $index2}, @indexList)) {
						my ($chromosome, $position1, $position2, $strand, $readStrand, $alignment1, $alignment2) = @$_;
						my $junction = join("\t", $chromosome, $position1, $position2, $strand);
						my $readStrandNumber = "$readStrand$number";
						if($stranded eq '') {
							$junctionReadCountHash{$junction}->{$readStrandNumber} += 1;
							if($position2 <= $position1 || $printAll) {
								push(@{$junctionAlignmentListHash{$junction}}, [$alignment1, $alignment2]);
								push(@{$junctionReadNameListHash{$junction}}, [$readName, $number, $readStrand]) if($readNameFile ne '');
							}
						} elsif($stranded eq 'f' && $strand eq '+' && any {$_ eq $readStrandNumber} ('+0', '+1', '-2')) {
							$junctionReadCountHash{$junction}->{$readStrandNumber} += 1;
							if($position2 <= $position1 || $printAll) {
								push(@{$junctionAlignmentListHash{$junction}}, [$alignment1, $alignment2]);
								push(@{$junctionReadNameListHash{$junction}}, [$readName, $number, $readStrand]) if($readNameFile ne '');
							}
						} elsif($stranded eq 'f' && $strand eq '-' && any {$_ eq $readStrandNumber} ('-0', '-1', '+2')) {
							$junctionReadCountHash{$junction}->{$readStrandNumber} += 1;
							if($position2 <= $position1 || $printAll) {
								push(@{$junctionAlignmentListHash{$junction}}, [$alignment1, $alignment2]);
								push(@{$junctionReadNameListHash{$junction}}, [$readName, $number, $readStrand]) if($readNameFile ne '');
							}
						} elsif($stranded eq 'r' && $strand eq '+' && any {$_ eq $readStrandNumber} ('-0', '-1', '+2')) {
							$junctionReadCountHash{$junction}->{$readStrandNumber} += 1;
							if($position2 <= $position1 || $printAll) {
								push(@{$junctionAlignmentListHash{$junction}}, [$alignment1, $alignment2]);
								push(@{$junctionReadNameListHash{$junction}}, [$readName, $number, $readStrand]) if($readNameFile ne '');
							}
						} elsif($stranded eq 'r' && $strand eq '-' && any {$_ eq $readStrandNumber} ('+0', '+1', '-2')) {
							$junctionReadCountHash{$junction}->{$readStrandNumber} += 1;
							if($position2 <= $position1 || $printAll) {
								push(@{$junctionAlignmentListHash{$junction}}, [$alignment1, $alignment2]);
								push(@{$junctionReadNameListHash{$junction}}, [$readName, $number, $readStrand]) if($readNameFile ne '');
							}
						}
					}
				}
			}
		}
	}
}

sub getJunctionList {
	my ($breakpointList1, $breakpointList2, @indexList) = @_;
	my %differenceJunctionListHash = ();
	foreach my $index (@indexList) {
		foreach my $breakpoint1 (@$breakpointList1) {
			my ($chromosome1, $strand1, $positionList1, $positionBaseHash1) = @$breakpoint1;
			my $position1 = $positionList1->[$index];
			next if($position1 eq '');
			foreach my $breakpoint2 (@$breakpointList2) {
				my ($chromosome2, $strand2, $positionList2, $positionBaseHash2) = @$breakpoint2;
				my $position2 = $positionList2->[$index - 1];
				next if($position2 eq '');
				if($chromosome1 eq $chromosome2) {
					next if(all {$chromosomeLengthHash{$chromosome1} < min(keys %$_)} $positionBaseHash1, $positionBaseHash2);
					next if(any {$_->[0] eq $chromosome1 && (($position1 <= $_->[2] && $_->[1] <= $position2) || ($position2 <= $_->[2] && $_->[1] <= $position1))} @blockRegionList);
					my $site1 = '';
					$site1 = uc(substr($chromosomeSequenceHash{$chromosome1}, $position1, 2)) if($strand1 eq '+');
					$site1 = uc(substr($chromosomeSequenceHash{$chromosome1}, $position1 - 3, 2)) if($strand1 eq '-');
					my $site2 = '';
					$site2 = uc(substr($chromosomeSequenceHash{$chromosome2}, $position2, 2)) if($strand2 eq '+');
					$site2 = uc(substr($chromosomeSequenceHash{$chromosome2}, $position2 - 3, 2)) if($strand2 eq '-');
					my $alignment1 = getAlignment($strand1, $position1, %$positionBaseHash1);
					my $alignment2 = getAlignment($strand2, $position2, %$positionBaseHash2);
					if($strand1 eq '+' && $strand2 eq '-' && ("$site1$site2" eq 'GTAG')) {
						my $difference = abs($position2 - $position1);
						push(@{$differenceJunctionListHash{$difference}}, [$chromosome1, $position1, $position2, '+', '-', $alignment1, $alignment2]);
					}
					if($strand2 eq '+' && $strand1 eq '-' && ("$site2$site1" eq 'GTAG')) {
						my $difference = abs($position1 - $position2);
						push(@{$differenceJunctionListHash{$difference}}, [$chromosome1, $position2, $position1, '+', '+', $alignment2, $alignment1]);
					}
					if($strand1 eq '+' && $strand2 eq '-' && ("$site1$site2" eq 'CTAC')) {
						my $difference = abs($position2 - $position1);
						push(@{$differenceJunctionListHash{$difference}}, [$chromosome1, $position1, $position2, '-', '-', $alignment1, $alignment2]);
					}
					if($strand2 eq '+' && $strand1 eq '-' && ("$site2$site1" eq 'CTAC')) {
						my $difference = abs($position1 - $position2);
						push(@{$differenceJunctionListHash{$difference}}, [$chromosome1, $position2, $position1, '-', '+', $alignment2, $alignment1]);
					}
				}
			}
		}
	}
	if(my @differenceList = keys %differenceJunctionListHash) {
		my $difference = min(@differenceList);
		return @{$differenceJunctionListHash{$difference}};
	} else {
		return ();
	}
}

sub getAlignment {
	my ($strand, $position, %positionBaseHash) = @_;
	return join('', map {defined($_) ? $_ : '-'} @positionBaseHash{min(keys %positionBaseHash) .. $position}) if($strand eq '+');
	return join('', map {defined($_) ? $_ : '-'} @positionBaseHash{$position .. max(keys %positionBaseHash)}) if($strand eq '-');
	return '';
}

sub getPositionList {
	my ($position, $cigar) = @_;
	my @positionList = ();
	my $index = 0;
	while($cigar =~ s/^([0-9]+)([MIDNSHP=X])//) {
		my ($length, $operation) = ($1, $2);
		if($operation eq 'M') {
			@positionList[$index .. $index + $length - 1] = $position .. $position + $length - 1;
			$index += $length;
			$position += $length;
		} elsif($operation eq 'I') {
			@positionList[$index .. $index + $length - 1] = ('') x $length;
			$index += $length;
		} elsif($operation eq 'D') {
			$position += $length;
		} elsif($operation eq 'N') {
			$position += $length;
		} elsif($operation eq 'S') {
			@positionList[$index .. $index + $length - 1] = ('') x $length;
			$index += $length;
		}
	}
	return @positionList;
}

sub getReverseComplementarySequence {
	my ($sequence) = @_;
	($sequence = reverse($sequence)) =~ tr/ACGT/TGCA/;
	return $sequence;
}

sub compare {
	my ($a, $b) = @_;
	my @a = @$a;
	my @b = @$b;
	if(scalar(@a) > 0 && scalar(@b) > 0) {
		$a = shift @a;
		$b = shift @b;
		if(looks_like_number($a) && looks_like_number($b)) {
			return $a <=> $b || compare(\@a, \@b);
		} else {
			return $a cmp $b || compare(\@a, \@b);
		}
	} elsif(scalar(@a) > 0) {
		return 1;
	} elsif(scalar(@b) > 0) {
		return -1;
	} else {
		return 0;
	}
}
