# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use GD;
use GD::Text::Align;
use List::Util qw(all);
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
	'l' => \(my $isLinearChromosome = ''),
	'r=s' => \(my $targetRegion = ''),
	't=s' => \(my $ttfFile = ''),
	'g=s' => \(my $gffFile = ''),
	'f=s' => \(my $feature = 'gene'),
	'a=s' => \(my $attribute = 'gene'),
	'X=f' => \(my $multiplex = 0.1),
	'Y=f' => \(my $multipley = 50),
	'T=i' => \(my $boxExtraThick = 0),
	'F=i' => \(my $boxFontSize = 10),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl vircircRNA_diagram.pl [options] junction.txt chromosome.fasta > diagram.png

Options: -h       display this help message
         -l       is linear chromosome
         -r STR   target region
         -t STR   TTF file
         -g STR   GTF file
         -f STR   GTF feature [$feature]
         -a STR   GTF attribute [$attribute]
         -X FLOAT multiple of x [$multiplex]
         -Y FLOAT multiple of y [$multipley]
         -T INT   box extra thick [$boxExtraThick]
         -F INT   box font size [$boxFontSize]

EOF
}

my ($junctionFile, $fastaFile) = @ARGV;

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

my $targetChromosome = $chromosomeList[0];
my ($targetStart, $targetEnd) = (1, $chromosomeLengthHash{$targetChromosome});
($targetChromosome, $targetStart, $targetEnd) = ($1, $2, $3) if($targetRegion =~ /^(.*):([0-9]+)-([0-9]+)$/);
if($isLinearChromosome eq '') {
	$chromosomeLengthHash{$_} = $chromosomeLengthHash{$_} / 2 foreach(@chromosomeList);
}
my $chromosomeLength = $chromosomeLengthHash{$targetChromosome};

my @geneStartEndListList = ();
if($gffFile) {
	my @geneStartEndList = ();
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
		if($tokenHash{'chromosome'} eq $targetChromosome && $tokenHash{'feature'} eq $feature) {
			if(defined(my $gene = $attributeHash{$attribute})) {
				{
					my ($start, $end) = @tokenHash{'start', 'end'};
					push(@geneStartEndList, [$gene, $start, $end]);
				}
				if($isLinearChromosome eq '') {
					my ($start, $end) = map {$_ + $chromosomeLength} @tokenHash{'start', 'end'};
					push(@geneStartEndList, [$gene, $start, $end]);
				}
			}
		}
	}
	close($reader);
	@geneStartEndList = sort {$a->[1] <=> $b->[1] || $a->[2] <=> $b->[2]} @geneStartEndList;

	foreach(@geneStartEndList) {
		my ($gene, $start, $end) = @$_;
		my $index = 0;
		for(; $index < scalar(@geneStartEndListList); $index += 1) {
			last if($geneStartEndListList[$index]->[-1]->[2] < $start);
		}
		push(@{$geneStartEndListList[$index]}, [$gene, $start, $end]);
	}
}

my $image = new GD::Image(getPixel($targetEnd) + 1 + $boxExtraThick, $multipley * (2 + scalar(@geneStartEndListList)) + 1 + $boxExtraThick);
my $white = getColor($image, '#FFFFFF');
my $black = getColor($image, '#000000');

my $y = $multipley * 2;
if(@geneStartEndListList) {
	foreach(@geneStartEndListList) {
		my ($y1, $y2) = ($y, $y + $multipley);
		foreach(@$_) {
			my ($gene, $start, $end) = @$_;
			my ($x1, $x2) = map {getPixel($_)} ($start, $end);
			boxText($image, $x1, $y1, $x2, $y2, $black, $gene);
		}
		$y = $y + $multipley;
	}
}

{
	open(my $reader, $junctionFile);
	while(my $line = <$reader>) {
		chomp($line);
		my ($chromosome, $position1, $position2, $strand, $thick, $color) = split(/\t/, $line);
		$thick = log($thick) / log(2);
		$color = getColor($image, $color);
		if($chromosome eq $targetChromosome) {
			{
				my ($x1, $x2) = map {getPixel($_)} sort {$a <=> $b} ($position1, $position2);
				my ($y1, $y2) = ($multipley, $multipley * 2);
				if($position1 < $position2) {
					angleThick($image, $x1, $y1, $x2, $y2, $thick, $color);
				} else {
					arcThick($image, $x1, $y1, $x2, $y2, $thick, $color);
				}
			}
			if($isLinearChromosome eq '' && $position1 <= $chromosomeLength && $position2 <= $chromosomeLength) {
				my ($x1, $x2) = map {getPixel($_)} map {$_ + $chromosomeLength} sort {$a <=> $b} ($position1, $position2);
				my ($y1, $y2) = ($multipley, $multipley * 2);
				if($position1 < $position2) {
					angleThick($image, $x1, $y1, $x2, $y2, $thick, $color);
				} else {
					arcThick($image, $x1, $y1, $x2, $y2, $thick, $color);
				}
			}
			if($isLinearChromosome eq '' && $position1 > $chromosomeLength && $position2 > $chromosomeLength) {
				my ($x1, $x2) = map {getPixel($_)} map {$_ - $chromosomeLength} sort {$a <=> $b} ($position1, $position2);
				my ($y1, $y2) = ($multipley, $multipley * 2);
				if($position1 < $position2) {
					angleThick($image, $x1, $y1, $x2, $y2, $thick, $color);
				} else {
					arcThick($image, $x1, $y1, $x2, $y2, $thick, $color);
				}
			}
		}
	}
	close($reader);
}

binmode STDOUT;
print $image->png(0);

sub angleThick {
	my ($image, $x1, $y1, $x2, $y2, $thick, $color) = @_;
	my ($cx, $cy) = (int(($x1 + $x2) / 2), $y2);
	my $width = ($cx - $x1) * 2 + 1;
	foreach(0 .. $thick) {
		$image->line($x1, $cy, $cx, $y1 - $_, $color);
		$image->line($x2, $cy, $cx, $y1 - $_, $color);
	}
	foreach(0 .. $thick) {
		$image->line($x1, $cy, $cx, $y1 + $_, $color);
		$image->line($x2, $cy, $cx, $y1 + $_, $color);
	}
}

sub arcThick {
	my ($image, $x1, $y1, $x2, $y2, $thick, $color) = @_;
	my ($cx, $cy) = (int(($x1 + $x2) / 2), $y2);
	foreach(0 .. $thick) {
		my $width = ($cx - $x1) * 2 + 1;
		my $height = ($y2 - ($y1 - $_)) * 2 + 1;
		$image->arc($cx, $cy, $width, $height, 180, 360, $color);
	}
	foreach(0 .. $thick) {
		my $width = ($cx - $x1) * 2 + 1;
		my $height = ($y2 - ($y1 + $_)) * 2 + 1;
		$image->arc($cx, $cy, $width, $height, 180, 360, $color);
	}
}

sub boxText {
	my ($image, $x1, $y1, $x2, $y2, $color, $text) = @_;
	foreach(0 .. $boxExtraThick) {
		$image->rectangle($x1 - $_, $y1 - $_, $x2 + $_, $y2 + $_, $color);
		$image->rectangle($x1 + $_, $y1 + $_, $x2 - $_, $y2 - $_, $color);
	}

	my $align = GD::Text::Align->new($image, valign => 'center', halign => 'center', color => $color);
	$align->set_font($ttfFile, $boxFontSize);
	$align->set_text($text);
	$align->draw(int(($x1 + $x2) / 2), int(($y1 + $y2) / 2), 0);
}

sub getColor {
	my ($image, $color) = @_;
	if($color =~ /^#([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})([0-9A-Fa-f]{2})$/) {
		return $image->colorAllocate(map {hex} ($1, $2, $3));
	} else {
		my $red = 255 - int((1 - $color) * 256);
		$red = 0 if($red < 0);
		$red = 255 if($red > 255);
		return $image->colorAllocate($red, 0, 0);
	}
}

sub getPixel {
	return sprintf('%.0f', ($_[0] - $targetStart) * $multiplex + $boxExtraThick);
}
