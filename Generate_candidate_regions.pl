=cut 
	function:
	1. split CpG island
	2. The region out of CpG island containing dense CpGs  
	Input: 1. ref; 2. CpG island inf; 3. paraments distance k
	Ouput: CpG cluster
=cut

use strict;

my $CpG_island_inf=shift;
my $hg_38=shift;
my $k_flanking=20;
my $out=$k_flanking."_Final_CpG_cluster.txt";
open(OUT,">$out")or die $!;

my %CpG_Cluster=();
open(IN,$CpG_island_inf)or die $!;
while(my $line=<IN>){
	chomp $line;
	my @temp=split /\t/,$line;
	my $sub=int($temp[6]/10);
	my $left=$temp[6]%10;
	if($left>6){
		$sub++;
	}
	my $region=int(($temp[3]-$temp[2]+1)/$sub);
	for(my $i=0;$i<$sub-1;$i++){
		$CpG_Cluster{$temp[1]}{($temp[2]+$i*$region)}=$temp[2]+($i+1)*$region-1;
		#print $temp[1]."\t".($temp[2]+$i*$region)."\t".($temp[2]+($i+1)*$region-1)."\n";
	}
	$CpG_Cluster{$temp[1]}{($temp[2]+($sub-1)*$region)}=$temp[3];
}
close(IN);

my $chr_CF_count=();
my @region_start=();
my $chr="";
my $pos=0;
my $last_base="";
my $last_CF=-1;
my $CF_id=0;

my %new_cluster=();
my %cpg_island_cluster_flag=();
my $new_id=0;
my $former_id=0;
my $chr_end=0;
open(IN,$hg_38)or die $!;
while(my $line=<IN>){
	chomp $line;
	if($line=~/^>(.+)$/){
		print $line;
		if($chr ne ""){
			print $chr."\t".$chr_end."\n";
			&print_CF($chr_end); 
		}
		$chr=$1;
		print $chr."\n";
		$pos=0;
		$chr_CF_count=0;		
		$last_base="";
		@region_start=sort{$a<=>$b}keys(%{$CpG_Cluster{$chr}});
		print "CpG cluster_num:\t".@region_start."\n";
		%new_cluster=();
		%cpg_island_cluster_flag=();
		$CF_id=0;
		$new_id=0;
		$former_id=0;
		$chr_end=0;

	}else{
		
		my @temp=split //,$line;
		my $i=0;
		my $two_base="";
		my $cpg_pos=0;
		$chr_end+=@temp;
		if($last_base ne ""){
			$two_base=$last_base.$temp[$i];
			if($two_base=~/CG/i){
			#	$CF_pos{$pos+$i}=1;
			#	if($N>0){print $two_base."\t".($pos+$i)."\t".($pos+$i-$last_CF)."\n";$N--;}
				$cpg_pos=$pos-1;
				$chr_CF_count++;
				
				while($CF_id<@region_start && $cpg_pos>$CpG_Cluster{$chr}{$region_start[$CF_id]}){
					$CF_id++;
				}
				if($cpg_pos>=$region_start[$CF_id] && $cpg_pos<=$CpG_Cluster{$chr}{$region_start[$CF_id]}){
					
					$last_CF=-1;
				
					#print ($pos+$i-1)."\t".$region_start[$CF_id]."\t".$CpG_Cluster{$chr}{$region_start[$CF_id]}."\n"; 
				}else{
					if($last_CF!=-1){
						if($cpg_pos-$last_CF<=$k_flanking){
							push @{$new_cluster{$new_id}},$cpg_pos;
						}else{
							push @{$new_cluster{++$new_id}},$cpg_pos;
						}
					}else{
						push @{$new_cluster{++$new_id}},$cpg_pos;
					}
					$last_CF=$cpg_pos;
				}
				
				$i++;
			}
			
		}
		for($i;$i<@temp-1;$i++){
			$two_base=$temp[$i].$temp[$i+1];
			if($two_base=~/CG/i){
			#	$CF_pos{$pos+$i}=1;
			#	if($N--){print $two_base."\t".($pos+$i)."\t".($pos+$i-$last_CF)."\n";}
				$chr_CF_count++;
				$cpg_pos=$pos+$i;
				while($CF_id<@region_start && $cpg_pos>$CpG_Cluster{$chr}{$region_start[$CF_id]}){
					$CF_id++;
				}
				if($cpg_pos>=$region_start[$CF_id] && $cpg_pos<=$CpG_Cluster{$chr}{$region_start[$CF_id]}){
					
					$last_CF=-1;
					
					#print $cpg_pos."\n"; 
				}else{
					if($last_CF!=-1){
						if($cpg_pos-$last_CF<=$k_flanking){
							push @{$new_cluster{$new_id}},$cpg_pos;
						}else{
							push @{$new_cluster{++$new_id}},$cpg_pos;
						}
					}else{
						
						push @{$new_cluster{++$new_id}},$cpg_pos;
					}
					$last_CF=$cpg_pos;
				}
			}
		}
		$last_base=$temp[-1];
		$pos=$pos+@temp;
		
	}
	
}
if($chr ne ""){
	print $chr."\t".$chr_end."\n";
	&print_CF($chr_end); 
}

sub print_CF{
	my $end=shift;
	my @cluster=sort{$a<=>$b}keys(%new_cluster);
	my $gene_body_region_id=0;
	my $glob_id=0;
	my $CF_inCluster=0;
	my $chr_CF_CpG_count=0;
	my $start=0;
	my $last=0;
	for my $id(@cluster){
		my @size=sort{$a<=>$b}@{$new_cluster{$id}};
		if(@size>6){
			$CpG_Cluster{$chr}{$size[0]}=$size[-1]+1;
		}
	}
	
	my @cluster=sort{$a<=>$b}keys(%{$CpG_Cluster{$chr}});
	for my $start(@cluster){
		print OUT (++$glob_id)."\t".$chr."\t".$start."\t".$CpG_Cluster{$chr}{$start}."\t".($CpG_Cluster{$chr}{$start}-$start)."\n";
	}
}
close(IN);
close(OUT);
close(GR);

