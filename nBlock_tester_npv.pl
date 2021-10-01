#!/usr/sbin/perl
####################################################################################################################################
# Script by Aimer G. Diaz (2018)
# Assited by Clara I Bermudez and Steve Hoffmann 
# 
#
# Run as 
# perl blockbuster_Gaussian_test.pl  12d4_S7.mapped.ncRNA.bed 12d4_S7.mapped.ncRNA.bam
# This is the main code non-paralel version of NBlock_tester, try it with your own non-redundant or tag reduced bam files and an 
# annotation file in bed format, to run in parallel the perl package  Parallel::ForkManager is neded, in that case delete commnetsi
#
# The main purpose is discriminate ncRNA and intergenic source expression from small RNA patterns derived from the same locus
####################################################################################################################################
use strict;
use warnings;
use List::Util qw(sum);
use List::Util qw( reduce );
use Getopt::Long;
#use Data::Dumper;
#use List::MoreUtils qw(uniq);
#use Parallel::ForkManager;

#########################################################################################################################
##---------------------------------------Data Parallelism----------------------------------------------------------------
##------------------------------To read and save  a heavy file into a  single hash using several cpus per fork  ---------
#########################################################################################################################
#uncomment to parallelize 
#my $forks = shift @ARGV;
#my $forks = $ARGV[0]; 
#print "$forks and $ARGV[0]\n"; 
##---------------------------------------To read and save  a heavy file into a  single hash using several cpus  ----------------------------------
#by division of read lines 
#my $line_count = `wc -l  $ARGV[0] | awk '{print \$1}'`;
#my $parts = ( int($line_count/$forks) + 1) ;  

#my %result; 
#my %Bed; 
#my @sub_hashes;
#my @q; 
#my $pm = Parallel::ForkManager->new($forks);
#$pm->run_on_finish( sub {
#my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
#print "** $ident finish with PID $pid and exit code: $exit_code\n";
#print "** $ident started, pid: $pid\n";
#print "** $ident finish with PID $pid and exit code: $exit_code\n";
# my $q =  $data_structure_reference->{input} ;
#  $Bed{$q} = $data_structure_reference->{output};
# %Bed = ( %$data_structure_reference , %Bed) ;
# push @sub_hashes, $data_structure_reference;
#}
#);

#my $inf_limit = 0; 
#inf=0; for f  in 1681 3200 ; do awk -v lim="$f" -v inf="$inf" '{if (NR > inf && NR < lim) print $0}' LP115-1.R1.ncRNA.bed | wc  ; inf=$f; done 
#########################################################################################################################
##---------------------------------------Data Feed----------------------------------------------------------------------
##------------------------------Open the coordinates of intergenic or ncRNA expressions---------------------------------
#########################################################################################################################

my $name = $ARGV[0];
 $name =~ s/.*\/(.*)\.bed/$1/;
my $annotfile  = $ARGV[1];; 

$annotfile  =~ s/.*\/(.*).bam/$1/;
print "\nWorking on ".$annotfile." bam file with non-redundant reads (or read tags) mapped against human genome version hg19. Reads alignment tool used segemehl.\nThe coordinates to check block patterns of expression with small RNA like features is in ".$name." bed file\n\n";
#my $denv_lib =~  s/.*\/(DENV.*)\/.*/$1/; 
#print $denv_lib."\n";  
my %Bed; 
#for my $i (1 .. $forks) {
#my $sup_limit =  $i * $parts;
#push @coordinates, $inf_limit."_".$sup_limit;
#$inf_limit = $sup_limit + 1 ;
#$sup_limit =  $parts / $i;
#} 
my %new_start ;
my $wd = "./SmallRNA_patterns_unnanotated/";
open my $OUT, '>>', $wd.$name.".sm-blocks.bed";

#Parallel:
#for my $i (1 .. $forks) {
#my $pid = $pm->start($i) and next Parallel;
#print "loop $i\n";
#my @coor =  split /_/,  $coordinates[$i -1];
#my $inf_limit = $coor[0]  ;
#my $sup_limit = $coor[1]; 
#print $inf_limit." vs  $sup_limit\n";
#open(BED, "awk -v lim_inf=$inf_limit -v lim_sup=$sup_limit '{if (NR > lim_inf && NR < lim_sup ) print }' $ARGV[0] | ");
open(BED, "$ARGV[0]"); 
#my out = (); 
	while (<BED>){
	my $bed = $_;
	chomp $bed;
	my @split = split /\t/, $bed;
	my $chr = $split[0];
	my $start = $split[1];
        my $end = $split[2];
        my $feature = $split[3];
        my $S_count = $split[4];
        my $sense_o = $split[5];
 
	my $sense_filter ; 
	if ($sense_o eq "+") { 
	$sense_filter = "-F 20";
	} elsif ($sense_o eq "-") { 
	$sense_filter = "-f 16";
	}

my (%mean_ordered, %single_line);       
#########################################################################################################################
##-------------------------------------- opening of BAM file following region of BED open hash---------------------------
##------------------------------Use Samttols to recover, from non redundant tag count to  total read count and store ----
#########################################################################################################################
 
	open(BAM, "samtools view $sense_filter $ARGV[1] $chr:$start-$end |  "); 
	while (<BAM>){
        my $bam = $_;
        chomp $bam;
        my @sam = split /\t/, $bam;  ## splitting SAM line into array
        my $bam_end = $sam[5];
        $bam_end =~ s/[^0-9]/ /g ;
        my @numbers = split / /, $bam_end;
        my $sum = sum(@numbers );
	$bam_end = $sam[3] + $sum;
        my @real_count = split /_/, $sam[0];
        my $mu_i = ($bam_end+$sam[3])/2;
	my $read_name =  join('_', @real_count[0 .. ($#real_count-1)]);
	my $read = "$sam[2]\t$sam[3]\t$bam_end\t$read_name\t$S_count\t$real_count[$#real_count]\n"; 
	push @{ $mean_ordered{$mu_i} } , $read;
                }  #close while BAM
         
#warn Dumper \%mean_ordered;
foreach my $keys ( keys %mean_ordered) {
	foreach my $bed_line ( @{ $mean_ordered{$keys} }  ) {
	$single_line{$keys."\t".$bed_line} = "" ;
							    } # close foreach bed lines ordered by means 

					   } # close foreach bed 				
next if((scalar keys %single_line <= 1)); #discard ncRNA read groups with only one mean  
#warn Dumper \%single_line; 
#########################################################################################################################
##---------------------------------------Non redundant mapping re-count -------------------------------------------------
##------------------------------Take each ncRNA region and discriminate between source and smallRNA derivend pattern-----
#########################################################################################################################
my ($mu, @block, $dev, @dev, $range_rigth,  @control, @total_mu, @block_mean, @mu, %block_starts, %block_ends, $rstart, $rend );
my $count = 0; 
my $total_count = 0 ;
my $tag_count = 1;
my $last_line = (scalar keys %single_line)   ;
my $current_line=0 ;  
foreach my $keys (sort keys  %single_line ) {
	$current_line++; 
	my @correction = split /\t/, $keys;
	$mu = $correction[0];
	push @mu, $mu; 
	$rstart = $correction[2];
	$rend = $correction[3];
	$count = $correction[6];
	push  @total_mu, ($rstart, $rend);
	push @block_mean, ($rstart, $rend);
	$dev =  stdev(\@total_mu);
	push @dev, $dev;
	$range_rigth = $mu[0] + ($dev[0]/2); #range simmetry calculation 
	push @control, $range_rigth;
#	 print "$rstart - $rend"."\t".$mu."\t".$dev."\t".($dev/4)."\t".($rend-$rstart)."\t".$range_rigth."\n"; #uncomment to see the first read 
	
		if ($#control > 0) { # Nex if the read is the first
			if ( $mu < $range_rigth && $current_line <= $last_line ) { #check if the new mean is lower than the previous defined right limit 
			push @{ $block_starts{$rstart} } , $count;  #count how many times is the start present in the entire tags set 
			push @{ $block_ends{$rend} } , $count; 
			$total_count = $total_count +  $count;  
			my $real_mean = mean(\@block_mean); 
			push @mu , $real_mean; 
			$dev =  stdev(\@block_mean);
			push @dev, $dev;
			$tag_count ++  ; 
	#uncomment next lines to see the iterative process  
#	print "IN $rstart:$rend mean is inside the range $mu <$range_rigth. The new sd is $dev ($dev[0]) and ".($mu[2] +  ($dev[1]/3))." count is  $total_count\n";
#			if ($current_line == $last_line ) { 
#			print "The last mean is the cluster was included \n "; 
#			}  
			@total_mu =  @total_mu[ 2 .. 3 ]; #erase the initial start and end coordenate 
		        @control = $control[1]; #erase the initial rigth limit
			@dev = $dev[2];	
			@mu = $mu[2];
			}  else  { #the mean tested has a higher value than the rigth limit of the previous distribution 
		   	my $real_mean = sprintf("%.0f",$mu[0]);
#			push 
				if ($total_count >= 50) { #Filter distributions with lower counts than 50
				my (%final_start, %final_end) ;       
				for my $i (sort keys  %block_starts) {
				$final_start{$i} = sum(@{ $block_starts{$i} }); #sum the frequency of a given start                                                    
				} 
				my $final_start = reduce { $final_start{$a} > $final_start{$b} ? $a : $b } keys %final_start; #extract the start coordinate with the highest value 
			for my $i (sort keys  %block_ends) {
				$final_end{$i} = sum(@{ $block_ends{$i}  });	
	                	                          }
				my $final_final = reduce { $final_end{$a} > $final_end{$b} ? $a : $b } keys %final_end;
 
				my $raw_block_mean = $final_start+abs($final_start-$final_final)/2;
				@block_mean = sort  @block_mean;
				my $cute_sd= sprintf("%.2f",$dev[0]);
				if  ( ($final_final-$final_start) <= 16  ) {
 print $OUT "$correction[1]\t".int($real_mean-15)."\t".int($real_mean+15)."\t$feature\t$total_count\t$sense_o\tsl=$start:$end;dl=$block_mean[0]:$block_mean[$#block_mean];sc=$final_start{$final_start};ec=$final_end{$final_final};soc=$S_count;tc=$tag_count;rm=$real_mean;bm=$raw_block_mean;sd=$cute_sd\t".(int($real_mean+15)-int($real_mean-15))."\n";
				}  else {
print $OUT "$correction[1]\t$final_start\t$final_final\t$feature\t$total_count\t$sense_o\tsl=$start:$end;dl=$block_mean[0]:$block_mean[$#block_mean];sc=$final_start{$final_start};ec=$final_end{$final_final};soc=$S_count;tc=$tag_count;rm=$real_mean;bm=$raw_block_mean;sd=$cute_sd\t".($final_final-$final_start)."\n";
				        } 
							} #close if heigth filter 
			$total_count = 0; #re-start the counter
			undef(@block_mean);  #re-start the means  array
			undef(%block_starts);  #re-start the start and end array 
			undef(%block_ends);
#			undef(@range_couples); 
			$tag_count = 1;
		#print "-----------OUT $rstart-$rend  $mu > $range_rigth (  $dev vs $dev[0] ). New test with  (".($mu +  ($dev[1]/3)).")\n";#." < ".($dev[$#dev]/2)."\n";
		@total_mu =  @total_mu[ 2 .. 3 ]; #keep only the last added start and end 
		@control = $control[1];
		@dev = $dev[1];
		@mu = $mu[1];
#
			} # close if check  $mu < $range_couples[0] 	
#
		} #close if for test if is the first occurrence 
} #close foreach explotarion of the bed file 
} #CLOSE INITIAL BED
# 
##------------------------ Subroutines --------------
##
sub mean {
    my ($arrayref) = @_;
    my $result;
    foreach (@$arrayref) { $result += $_ }
    return $result / @$arrayref;
}

sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &mean($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

