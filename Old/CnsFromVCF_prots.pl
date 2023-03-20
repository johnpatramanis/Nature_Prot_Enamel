######################################
# Last update Oct 11th 2018
# perl CnsFromVCF_prots.pl -i Aymara_22.vcf -s Aymara -chr 22
#
# IMPORTANT NOTES:
# creates a fasta file from a VCF format: VCFv4.2
# does not take into account INDELS
# keep sites with up to 3 alternative alleles 
######################################

# chr17   48245650        .       G       .       112.23  .       AN=282;DP=3650;set=ReferenceInAll       GT:DP:MQ:MQ0:Q  0/0:29:60.00:0:112.23

use Getopt::Long;
my %opts=();
GetOptions (\%opts,'i=s', 'chr=s', 's=s');

%hetmap=();
$hetmap{"AC"}="M";
$hetmap{"AG"}="R";
$hetmap{"AT"}="W";
$hetmap{"CA"}="M";
$hetmap{"CG"}="S";
$hetmap{"CT"}="Y";
$hetmap{"GA"}="R";
$hetmap{"GC"}="S";
$hetmap{"GT"}="K";
$hetmap{"TA"}="W";
$hetmap{"TC"}="Y";
$hetmap{"TG"}="K";


$chr=$opts{chr};
$GREPCMD="grep -m 1 ID=".$chr.", ".$opts{i};

$contigline=`$GREPCMD`;
if($contigline=~/length=(\d+)/){
	$contiglength=$1;
}
#print $contiglength;

@CNS=();
@RAS=();
for($i=0; $i<$contiglength; $i++){
	$CNS[$i]="N";
	$RAS[$i]="N";
}

open(IN, $opts{i});
while(<IN>){
	if($_ =~ /#/){
		next;
	}
	$line=$_;
	chomp($line);
	@l=split(/\t/, $line);
	$gt=$l[9];
	@gtline=split(/:/, $gt);
	$gt=$gtline[0];
	$k=$l[1];
	$ref=$l[3];
	$alt=$l[4];

	# if genotype missing
	if($gt eq "./."){
		next;
	}

	# this part is to skip insertions and deletions:
	$lenR=length $ref;
	$lenA=length $alt;
	#print "ref $lenR $ref\t";
	#print "alt $lenA $alt";

	##count the number of alternative alleles
	$nalleles=0;
	if($lenR==1 && $lenA==1){
		#print "len1 $lenR \n";
		$nalleles=1;
	}elsif($alt=~/,/ && $lenR==1){
		#print "len2 $lenR \n";
		@altline=split(/,/, $alt);
		foreach $site (@altline){
			$nalleles++;
			$lenSite=length $site;
			if($lenSite>1){
				$nalleles+=5; #add 5 so that the number of alt alleles is >3 and skipped
			}
		}
	}

	$r=rand();

	if($nalleles==1){
		#print "passed 1\n";
		# homozygous reference (0/0)
		if($gt eq "0/0"){
			$CNS[$k-1]=$ref;
			$RAS[$k-1]=$ref;
			#print "$ref \n"; 
		# homozygous alternative (1/1)
		}elsif($gt eq "1/1"){
			$CNS[$k-1]=$alt;
			$RAS[$k-1]=$alt;
			#print "$alt \n"; 
		# hererozygous (0/1)
		}elsif($gt eq "0/1"){
			$CNS[$k-1]=$hetmap{$ref.$alt};
			if($r>.5){
				$RAS[$k-1]=$ref;
			}else{
				$RAS[$k-1]=$alt;			
			}
			#print "$hetmap{$ref.$alt} \n"; 
		# Unknown genotype, to check if there is something wrong
		}else{
			#$CNS[$k-1]="-";
			#print "$line"; 
			#print $gt; 
			#print " genotype unknown 1\n";
		}

	}elsif($nalleles>0 && $nalleles<=3){
		#print "passed 2\n";
		# homozygous reference (0/0)
		if($gt eq "0/0"){
			$CNS[$k-1]=$ref;
			$RAS[$k-1]=$ref;
			#print "$ref \n"; 
		# homozygous alternative 1 (1/1)
		}elsif($gt eq "1/1"){
			$CNS[$k-1]=$altline[0];
			$RAS[$k-1]=$altline[0];
			#print "$altline[0] \n"; 
		# homozygous alternative 2 (2/2)
		}elsif($gt eq "2/2"){
			$CNS[$k-1]=$altline[1];
			$RAS[$k-1]=$altline[1];
			#print "$altline[1] \n"; 
		# homozygous alternative 3 (3/3)
		}elsif($gt eq "3/3"){
			$CNS[$k-1]=$altline[2];
			$RAS[$k-1]=$altline[2];
			#print "$altline[2] \n"; 
		# heterozygous (0/1)
		}elsif($gt eq "0/1"){
			$CNS[$k-1]=$hetmap{$ref.$altline[0]};
			if($r>.5){
				$RAS[$k-1]=$ref;
			}else{
				$RAS[$k-1]=$altline[0];			
			}
			#print "$hetmap{$ref.$altline[0]} \n"; 
		# heterozygous with alternative 2 (0/2)
		}elsif($gt eq "0/2"){
			$CNS[$k-1]=$hetmap{$ref.$altline[1]};
			if($r>.5){
				$RAS[$k-1]=$ref;
			}else{
				$RAS[$k-1]=$altline[1];			
			}
			#print "$hetmap{$ref.$altline[1]} \n"
		# heterozygous with alternative 3 (0/3)
		}elsif($gt eq "0/3"){
			$CNS[$k-1]=$hetmap{$ref.$altline[2]};
			if($r>.5){
				$RAS[$k-1]=$ref;
			}else{
				$RAS[$k-1]=$altline[2];			
			}
			#print "$hetmap{$ref.$altline[2]} \n";
		# heterozygous with alternatives 1 and 2 ( (1/2)
		}elsif($gt eq "1/2"){
			$CNS[$k-1]=$hetmap{$altline[0].$altline[1]};
			if($r>.5){
				$RAS[$k-1]=$altline[0];
			}else{
				$RAS[$k-1]=$altline[1];			
			}
			#print "$hetmap{$altline[0].$altline[1]} \n";
		# heterozygous with alternatives 1 and 3 ( (1/3)
		}elsif($gt eq "1/3"){
			$CNS[$k-1]=$hetmap{$altline[0].$altline[2]};
			if($r>.5){
				$RAS[$k-1]=$altline[0];
			}else{
				$RAS[$k-1]=$altline[2];			
			}
			#print "$hetmap{$altline[0].$altline[2]} \n";
		# heterozygous with alternatives 2 and 3 ( (2/3)
		}elsif($gt eq "2/3"){
			$CNS[$k-1]=$hetmap{$altline[1].$altline[2]};
			if($r>.5){
				$RAS[$k-1]=$altline[1];
			}else{
				$RAS[$k-1]=$altline[2];			
			}
			#print "$hetmap{$altline[1].$altline[2]} \n";
		# Unknown genotype, to check if there is something wrong
		}else{
			#$CNS[$k-1]="-";
			#print "$line"; 
			#print $gt; 
			#print " genotype unknown 2\n";
		}
	}else{
		#print "dont passed\n";	
	}
}

close(IN);

open(HET, ">".$opts{s}."_".$chr."_"."HET".".fa");
open(RAND, ">".$opts{s}."_".$chr."_"."RAND".".fa");

print HET ">".$chr."\n";
print RAND ">".$chr."\n";
$cont=0;
for($i=0; $i<$contiglength; $i++){
	print HET $CNS[$i];
	print RAND $RAS[$i];
	$cont++;
	if($cont==60){
		print HET "\n";
		print RAND "\n";
		$cont=0;
	}
}
print HET "\n";
print RAND "\n";

close(HET);
close(RAND);
