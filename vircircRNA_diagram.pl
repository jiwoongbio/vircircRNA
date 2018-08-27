# Author: Jiwoong Kim (jiwoongbio@gmail.com)
use strict;
use warnings;
local $SIG{__WARN__} = sub { die $_[0] };

use GD;
use GD::Text::Align;
use Getopt::Long qw(:config no_ignore_case);

GetOptions(
	'h' => \(my $help = ''),
	'c' => \(my $isCircularChromosome = ''),
	't=s' => \(my $ttfFile = ''),
	'g=s' => \(my $gtfFile = ''),
	'f=s' => \(my $feature = 'exon'),
	'a=s' => \(my $attribute = 'gene_id'),
	'X=f' => \(my $multiplex = 0.1),
	'Y=f' => \(my $multipley = 50),
	'L=f' => \(my $multipleLineThick = 0.01),
	'T=i' => \(my $boxExtraThick = 0),
	'F=i' => \(my $boxFontSize = 10),
);
if($help || scalar(@ARGV) == 0) {
	die <<EOF;

Usage:   perl vircircRNA_diagram.pl [options] junction.txt reference.fasta chromosome f|forward|r|reverse > diagram.png

Options: -h       display this help message
         -c       is circular chromosome
         -t STR   TTF file
         -g STR   GTF file
         -f STR   GTF feature [$feature]
         -a STR   GTF attribute [$attribute]
         -X FLOAT multiple of x [$multiplex]
         -Y FLOAT multiple of y [$multipley]
         -L FLOAT multiple of line thick [$multipleLineThick]
         -T INT   box extra thick [$boxExtraThick]
         -F INT   box font size [$boxFontSize]

EOF
}
my ($junctionFile, $referenceFastaFile, $targetChromosome, $targetStrand) = @ARGV;
$targetStrand = '+' if($targetStrand eq 'f' || $targetStrand eq 'forward');
$targetStrand = '-' if($targetStrand eq 'r' || $targetStrand eq 'reverse');

my $chromosomeLength = 0;
{
	my $chromosome = '';
	open(my $reader, $referenceFastaFile);
	while(my $line = <$reader>) {
		chomp($line);
		if($line =~ /^>(\S*)/) {
			$chromosome = $1;
		} else {
			$chromosomeLength += length($line) if($chromosome eq $targetChromosome);
		}
	}
	close($reader);
}
$chromosomeLength = $chromosomeLength / 2 if($isCircularChromosome);

my @geneStartEndListList = ();
if($gtfFile) {
	my @geneStartEndList = ();
	open(my $reader, $gtfFile);
	while(my $line = <$reader>) {
		chomp($line);
		next if($line =~ /^#/);
		my %tokenHash = ();
		@tokenHash{'chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'} = split(/\t/, $line);
		my %attributeHash = ();
		$attributeHash{$1} = $2 while($tokenHash{'attribute'} =~ m/([^"; ]+) +"([^"]+)";/g);
		if($tokenHash{'chromosome'} eq $targetChromosome && $tokenHash{'strand'} eq $targetStrand && $tokenHash{'feature'} eq $feature) {
			if(defined(my $gene = $attributeHash{$attribute})) {
				{
					my ($start, $end) = @tokenHash{'start', 'end'};
					push(@geneStartEndList, [$gene, $start, $end]);
				}
				if($isCircularChromosome) {
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

my $image = new GD::Image(getPixel($isCircularChromosome ? $chromosomeLength * 2 : $chromosomeLength) + 1 + $boxExtraThick * 2, $multipley * (2 + scalar(@geneStartEndListList)) + 1 + $boxExtraThick);
my $white = $image->colorAllocate(255, 255, 255);
my $black = $image->colorAllocate(0, 0, 0);
my $red = $image->colorAllocate(255, 0, 0);

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
		my ($chromosome, $position1, $position2, $strand, $count, $ratio) = split(/\t/, $line);
		my $thick = $count * $multipleLineThick;
		if($chromosome eq $targetChromosome && $strand eq $targetStrand) {
			{
				my ($x1, $x2) = map {getPixel($_)} sort {$a <=> $b} ($position1, $position2);
				my ($y1, $y2) = ($multipley, $multipley * 2);
				if($position1 < $position2) {
					angleThick($image, $x1, $y1, $x2, $y2, $black, $thick);
				} else {
					arcThick($image, $x1, $y1, $x2, $y2, $red, $thick);
				}
			}
			if($isCircularChromosome && $position1 <= $chromosomeLength && $position2 <= $chromosomeLength) {
				my ($x1, $x2) = map {getPixel($_)} map {$_ + $chromosomeLength} sort {$a <=> $b} ($position1, $position2);
				my ($y1, $y2) = ($multipley, $multipley * 2);
				if($position1 < $position2) {
					angleThick($image, $x1, $y1, $x2, $y2, $black, $thick);
				} else {
					arcThick($image, $x1, $y1, $x2, $y2, $red, $thick);
				}
			}
			if($isCircularChromosome && $position1 > $chromosomeLength && $position2 > $chromosomeLength) {
				my ($x1, $x2) = map {getPixel($_)} map {$_ - $chromosomeLength} sort {$a <=> $b} ($position1, $position2);
				my ($y1, $y2) = ($multipley, $multipley * 2);
				if($position1 < $position2) {
					angleThick($image, $x1, $y1, $x2, $y2, $black, $thick);
				} else {
					arcThick($image, $x1, $y1, $x2, $y2, $red, $thick);
				}
			}
		}
	}
	close($reader);
}

binmode STDOUT;
print $image->png(0);

sub angleThick {
	my ($image, $x1, $y1, $x2, $y2, $color, $thick) = @_;
	my ($cx, $cy) = (int(($x1 + $x2) / 2), $y2);
	my $width = ($cx - $x1) * 2 + 1;
	foreach(map {$y1 - $_} 0 .. $thick) {
		$image->line($x1, $cy, $cx, $_, $color);
		$image->line($x2, $cy, $cx, $_, $color);
	}
	foreach(map {$y1 + $_} 0 .. $thick) {
		$image->line($x1, $cy, $cx, $_, $color);
		$image->line($x2, $cy, $cx, $_, $color);
	}
}

sub arcThick {
	my ($image, $x1, $y1, $x2, $y2, $color, $thick) = @_;
	my ($cx, $cy) = (int(($x1 + $x2) / 2), $y2);
	my $width = ($cx - $x1) * 2 + 1;
	foreach(map {$y1 - $_} 0 .. $thick) {
		my $height = ($y2 - $_) * 2 + 1;
		$image->arc($cx, $cy, $width, $height, 180, 360, $color);
	}
	foreach(map {$y1 + $_} 0 .. $thick) {
		my $height = ($y2 - $_) * 2 + 1;
		$image->arc($cx, $cy, $width, $height, 180, 360, $color);
	}
}

sub boxText {
	my ($image, $x1, $y1, $x2, $y2, $color, $text) = @_;
	foreach(0 .. $boxExtraThick) {
		$image->rectangle($x1 - $_, $y1 - $_, $x2 + $_, $y2 + $_, $black);
		$image->rectangle($x1 + $_, $y1 + $_, $x2 - $_, $y2 - $_, $black);
	}

	my $align = GD::Text::Align->new($image, valign => 'center', halign => 'center', color => $color);
	$align->set_font($ttfFile, $boxFontSize);
	$align->set_text($text);
	$align->draw(int(($x1 + $x2) / 2), int(($y1 + $y2) / 2), 0);
}

sub getPixel {
	return sprintf('%.0f', ($_[0] - 1) * $multiplex + $boxExtraThick);
}
