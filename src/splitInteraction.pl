use warnings;
use strict;
use File::Basename;

my $pdbfileIn=$ARGV[0];
my $dirOut=$ARGV[1];

my $basename=basename($pdbfileIn);
my $filename;
if($basename=~/(.+)\.pdb/){
	$filename=$1;
}

my @outNames = ("A", "B");
my $outCounter = 0;
my $flag=0;

open(INFILE, $pdbfileIn);
my @inlines = <INFILE>;
close(INFILE);

open(OUTFILE,">".$dirOut."/".$filename."_".$outNames[$outCounter].".pdb");
foreach my $inline (@inlines){
	if($flag==1){
		close(OUTFILE);
		open(OUTFILE,">".$dirOut."/".$filename."_".$outNames[$outCounter].".pdb");
		$flag=0;
	}
	
	print OUTFILE $inline;
	
	if($inline=~/^TER/){
		$flag++;
		$outCounter++;
		if($outCounter>1){
			last;
		}
	}
}
close(OUTFILE);