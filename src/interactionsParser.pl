use warnings;
use strict;

my $interactionsFile = $ARGV[0];
my $int3dDir=$ARGV[1];

print STDERR "Caching set of protein/interactions to read...\n";

open(INTERACTIONS, $interactionsFile);
my @intlines = <INTERACTIONS>;
close(INTERACTIONS);

my %hcol;
chomp($intlines[0]);
my @headerFields = split(/\t/, $intlines[0]);
for(my $j=0;$j<scalar(@headerFields);$j++){
	$hcol{$headerFields[$j]}=$j;
}

my %chainDict;
#Pre-reading all interactions and proteins
my @allinterFiles;
my @allprotFiles;
my %visitedInterFile;
my %visitedProtFile;
for(my $i=1;$i<scalar(@intlines);$i++){
	my $line = $intlines[$i];
	chomp($line);
	my @fields=split("\t", $line);
	
	my $intFile = $int3dDir."/interactions/rsas/".$fields[$hcol{"FILENAME"}];
	if($intFile=~/^(.+)\.pdb$/){
		$intFile=$1.".rsa";
	}
	my $p1File;
	my $p2File;
	my $fakeFile = $int3dDir."/proteins/rsas/".$fields[$hcol{"FILENAME"}];
	if($fakeFile=~/^(.+)\.pdb$/){
		$p1File=$1."_A.rsa";
		$p2File=$1."_B.rsa";
	}
  ${${$chainDict{$fields[$hcol{"FILENAME"}]}}{$fields[$hcol{"PROT1"}]}}{"A"}=$fields[$hcol{"CHAIN1"}];        
  ${${$chainDict{$fields[$hcol{"FILENAME"}]}}{$fields[$hcol{"PROT2"}]}}{"B"}=$fields[$hcol{"CHAIN2"}];
	if(-e $intFile){
		#files that exist and are not empty
		if(-f $intFile && -s $intFile){	
			if(! defined($visitedInterFile{$intFile})){
				push(@allinterFiles, $intFile);
				$visitedInterFile{$intFile}=1;
			}
		}
	}
	if(-e $p1File){
		#files that exist and are not empty
		if(-f $p1File && -s $p1File){
			if(! defined($visitedProtFile{$p1File})){
				push(@allprotFiles, $p1File);
				$visitedProtFile{$p1File}=1;
			}
		}
	}
	if(-e $p2File){
		#files that exist and are not empty
		if(-f $p2File && -s $p2File){
			if(! defined($visitedProtFile{$p2File})){
				push(@allprotFiles, $p2File);
				$visitedProtFile{$p2File}=1;
			}
		}
	}	
}

####################

print STDERR "Reading Protein files...\n";
my %proteins = %{batchReadProteins(\@allprotFiles)};

####################

print STDERR "Reading Interactions files...\n";
my %interactions = %{batchReadInteractions(\@allinterFiles)};

####################

print STDERR "Printing...\n";

my @header = ();
push(@header,"PROTEIN");
push(@header, "PARTNER");
push(@header, "AA");
push(@header, "POS");
push(@header, "ACCDIFFPERCENT");
push(@header, "STRUCT");
push(@header, "PDB");
push(@header, "CHAIN");
push(@header, "POS_PDB");
push(@header, "CHAIN_PARTNER");
push(@header, "INTERFACE");
print join("\t",@header)."\n";

for(my $i=1;$i<scalar(@intlines);$i++){
	my $line = $intlines[$i];
	chomp($line);
	
	my @fields=split("\t", $line);	
	#Interaction file
	my $intFile = $int3dDir."/interactions/rsas/".$fields[$hcol{"FILENAME"}];
	if($intFile=~/^(.+)\.pdb$/){
		$intFile=$1.".rsa";
	}
	my $p1File;
	my $p2File;
	my $fakeFile = $int3dDir."/proteins/rsas/".$fields[$hcol{"FILENAME"}];
	if($fakeFile=~/^(.+)\.pdb$/){
		$p1File=$1."_A.rsa";
		$p2File=$1."_B.rsa";
	}
	##PRINTING###############
	printAll($intFile, $p1File, $p2File, \%proteins, \%interactions, $fields[$hcol{"PROT1"}], $fields[$hcol{"PROT2"}], \%hcol, \@fields);
	printAll($intFile, $p2File, $p1File, \%proteins, \%interactions, $fields[$hcol{"PROT2"}], $fields[$hcol{"PROT1"}], \%hcol, \@fields);		
}



####################

#FUNCTIONS

sub readRSA{
	my $proteinFile=$_[0];
	
	open(INFILE, $proteinFile);
	my @inlines = <INFILE>;
	close(INFILE);
	
	my %interaction;
	foreach my $line (@inlines){
		chomp($line);
		if($line=~/^RES/){
			my $chain = substr($line, 8, 1);
			my $position = substr($line, 9, 4);
			if($position=~/^(\s+)?(.+)(\s+)?$/){
				$position=$2;
			}
			my $residue = substr($line, 4, 3);
			my $relativeAcc = substr($line, 22, 6);
			if($relativeAcc=~/^(\s+)?(.+)(\s+)?$/){
				$relativeAcc=$2;
			}
			push(@{${$interaction{$chain}}{"positions"}}, $position);
			${${$interaction{$chain}}{"aminoacids"}}{$position}=translateAA($residue);
			${${$interaction{$chain}}{"relAcc"}}{$position}=$relativeAcc;	
		}
	}
	return(\%interaction);
}

sub getInterface{
	my %protein=%{$_[0]};
	my %complex=%{$_[1]};
	
	my @proteinkeys = keys(%protein);
	if(scalar(@proteinkeys)>1){
		die "Multichain protein";
	}
	my %different;
	my $chain = $proteinkeys[0];
	foreach my $pos (@{${$protein{$chain}}{"positions"}}){
		my $diffAcc = ${${$protein{$chain}}{"relAcc"}}{$pos} - ${${$complex{$chain}}{"relAcc"}}{$pos}."\t";
		if($diffAcc != 0){
			$different{$pos}=$diffAcc;
		}		
	}
	return(\%different);
}

sub getSequence{
	my %protein = %{$_[0]};
	
	my @proteinkeys = keys(%protein);
	if(scalar(@proteinkeys)>1){
		die "Multichain protein";
	}
	my $chain = $proteinkeys[0];
	my @sequence;
	my @allsequences;
	my @allstarts;
	my $previous;
	my @residues;
	my $start ='';
	foreach my $pos (@{${$protein{$chain}}{"positions"}}){
		if(!defined($previous)){
			$previous=$pos;
		}else{
			if($pos != ($previous +1)){
				my $sequence=join("", @residues);
				push(@allsequences, $sequence);
				push(@allstarts, $start);
				
				$start='';
				@residues=();
			}
			$previous=$pos;
		}
		if($start eq ''){
			$start=$pos;
		}
		push(@residues, ${${$protein{$chain}}{"aminoacids"}}{$pos});
	}
	my $sequence=join("", @residues);
	push(@allsequences, $sequence);
	push(@allstarts, $start);
	
	my %seqInfo = ("sequences"=>\@allsequences,"starts"=>\@allstarts);
	
	return(\%seqInfo);

}

#it returns 0 when all the fragments are not found or a given fragment matches in more than one position. 
sub getPDBEnspTranslator{
	my %sequenceObject=%{$_[0]};
	my $ensSequence=$_[1];
	
	my %translator;
	
	my $match=0;
	for(my $j=0;$j<scalar(@{$sequenceObject{"sequences"}});$j++){
		my $fragment = ${$sequenceObject{"sequences"}}[$j];
		my $start = ${$sequenceObject{"starts"}}[$j];
		if($ensSequence=~/$fragment/){
			#if we find more than one match of the fragment we ignore the protein
			my @fragmentMatches =($ensSequence =~ /$fragment/g);
			if(scalar(@fragmentMatches)>1){
				return(0);
			}
			$match++;
			my @residues = split('', $fragment);
			for (my $i=0;$i<scalar(@residues);$i++){
				# print (($i+${$sequenceObject{"starts"}}[$j])."\t".($-[0]+$i+1)."\t".$residues[$i]."\n");
				$translator{$i+$start}=$-[0]+$i+1;
			}
			# print $fragment."\t".$start."\n";
		}
	}
	if(scalar(@{$sequenceObject{"sequences"}} == $match)){
		return(\%translator);
	}else{
		return(0);
	}
}

# It gets the percentage of protein that is covered by the pdb/rsa file
sub getCoverage{
	my %sequenceObject=%{$_[0]};
	my $ensemblSequence=$_[1];
	
	my @ensresidues = split("", $ensemblSequence);
	my $total = scalar(@ensresidues);
	
	my $totalInModel=0;
	foreach my $fragment (@{$sequenceObject{"sequences"}}){
		my @fragresidues = split("", $fragment);
		$totalInModel=$totalInModel+scalar(@fragresidues);
	}
	return($totalInModel/$total*100);
}

# It reads an rsa protein file and saves its important content
sub readProtein{
	my $proteinFile=$_[0];
	
	open(INFILE, $proteinFile);
	my @inlines = <INFILE>;
	close(INFILE);
	
	my %protein;
	foreach my $line (@inlines){
		chomp($line);
		if($line=~/^RES/){
			my @fields=split(/\s+/, $line);
			push(@{$protein{"positions"}}, $fields[3]);
			push(@{$protein{"aminoacids"}}, $fields[1]);
			push(@{$protein{"relAcc"}}, $fields[5]);			
		}
	}
	return(\%protein);
}

# Translate the three letter aminoacids intro one letter
sub translateAA{
	my $threeLetters = $_[0];
	my %translator = ('VAL'=>'V', 'ILE'=>'I', 'LEU'=>'L', 'GLU'=>'E', 'GLN'=>'Q',
	'ASP'=>'D', 'ASN'=>'N', 'HIS'=>'H', 'TRP'=>'W', 'PHE'=>'F', 'TYR'=>'Y',
	'ARG'=>'R', 'LYS'=>'K', 'SER'=>'S', 'THR'=>'T', 'MET'=>'M', 'ALA'=>'A',
	'GLY'=>'G', 'PRO'=>'P', 'CYS'=>'C', 'ASX'=>'B', 'GLX'=>'Z', 'XLE'=>'J',
	'SEC'=>'U', 'UNK'=>'X');
	my $oneLetter=$translator{$threeLetters};
	if(!defined($oneLetter)){
		die "Unknown aminoacid ".$threeLetters;
	}
	return($oneLetter);
}

# It reads all the Interaction RSA files and save the content on a hash in memory
sub batchReadInteractions{
	my @allIntFiles = @{$_[0]};
	
	my %interactions;
	foreach my $intFile (@allIntFiles){
		$interactions{$intFile} = readRSA($intFile);
	}
	return(\%interactions);
}

# It reads all the protein RSA files and save the content on a hash in memory
sub batchReadProteins{
	my @allProtFiles = @{$_[0]};
	my %proteins;
	foreach my $protFile (@allProtFiles){
		$proteins{$protFile} = readRSA($protFile);
	}
	return(\%proteins);
}

#Prints all the output
sub printAll{
	my $intFile=$_[0];
	my $p1File=$_[1];
	my $p2File=$_[2];
	my %proteins=%{$_[3]};
	my %interactions=%{$_[4]};
	my $p1=$_[5];
	my $p2=$_[6];
	my %hcol=%{$_[7]};
	my @fields=@{$_[8]};
		
	if(-e $intFile){
		#files that exist and are not empty
		if(-f $intFile && -s $intFile){
			my %complex = %{$interactions{$intFile}};
			my %p1protein = %{$proteins{$p1File}};
			my %p2protein = %{$proteins{$p2File}};  
			#Protein interface on protein based on the given interface
			my %interface = %{getInterface(\%p1protein,\%complex)};
	
			#Chain
			my @chains1 = keys(%p1protein);
			my $chain1 = $chains1[0];
			my @chains2 = keys(%p2protein);
			my $chain2 = $chains2[0];
      
			# print STDERR "seq\n";

			#For all the elements in the interface
			foreach my $intRes (keys %interface){
				if(defined($interface{$intRes})){
					my @toprint;
					push(@toprint, $p1);
					push(@toprint, $p2);
					push(@toprint, ${${$p1protein{$chain1}}{"aminoacids"}}{$intRes});
					push(@toprint, $intRes);
					push(@toprint, sprintf("%.2f", $interface{$intRes}));
					push(@toprint, $fields[$hcol{"TYPE"}]);
					push(@toprint, $fields[$hcol{"PDB_ID"}]);
					push(@toprint, ${${$chainDict{$fields[$hcol{"FILENAME"}]}}{$p1}}{$chain1});
					push(@toprint, $intRes);
					push(@toprint, ${${$chainDict{$fields[$hcol{"FILENAME"}]}}{$p2}}{$chain2});
					push(@toprint, $fields[$hcol{"FILENAME"}]);
		
					#Printing
					print join("\t", @toprint)."\n";
				}
			}
		}	
	}
}
