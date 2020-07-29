=cut header
	function: 1. Extract the methylation state of CpGs in reads 
			2. Reads falling into clusters
			3. Contruct haplotype for each clusters
			4. Label Clusters and Reads
			
	INPUT: sorted sam file, CpG_cluster_file
	OUTPUT：
			chromosome	marker_start	marker_end	Hap_start	Hap_end	Haplotype1	Haplotype 2
			chr1	23425332	23425432	23425332	23425385 1111000	000000
			
	Test: chr22_sorted_sam_file,chr22_cluster_file(waston/crick?)
	Process: 1. store the start and end coordinates of each cluster
			 2. classify each read to each cluster (may be a read fall into two cluster)
			 3. ouput the mehtylation states of reads fall in the clusters
	#First_version: 2019/1/30
	#Modified_version:
=cut

use strict;
use warnings;
use Getopt::Long;
use FileHandle;
my $help;
my $sam_file;
my $CF_file;

GetOptions(
	'help|man|h' => \$help,
	'sam:s' => \$sam_file,
	'cls:s' => \$CF_file,
);###Get input paramters

##### end options #####

#### input arguments ####
if($help){
	print_helpfile();
	die "\n";
}
if(!$sam_file){
	die "You need to specify the sam file.\n";
}
if(!$CF_file){
	die "You need to specify the cytosine cluster file.\n";
}

$CF_file=~/^(.{6})/;
my $sufix=$1;

#variables for cluster and reads
my %CF_Cluster=();
my %chr_hash=();
my %CG_Num=();
my %CG_Id=();
my @CFstart=();
my $next_i=0;
my %cluster_read=();
my $chr="";

#variables for constructing haplotypes
my %ploy_site=();
my %group=();
my %haplotype=();
my %heter=();
my $heter_depth=2;
my $homo_depth=4;

my $homo_mode_perc=0.7;
my $sec_to_max_perc=0.4;
my $out="Hap_".$sufix."_".$sam_file;
open(OUT,">>$out")or die $!;
#The sam file is sorted by single reads or paired reads.
open(SM,$sam_file)or die $!;
while(my $line=<SM>){
	if($line=~/^@/){
		next;
	}
	my @map_inf=split /\s+/,$line;
	if($map_inf[2] ne "*" && $map_inf[4]>=10){
		if(!exists($chr_hash{$map_inf[2]})){
			if(@CFstart>0 && exists($cluster_read{$CFstart[-1]}) && $cluster_read{$CFstart[-1]} ne ""){
				&haplotype($cluster_read{$CFstart[-1]},$CFstart[$CFstart[-1]]);
				$cluster_read{$CFstart[-1]}="";
			}
			
			$chr=$map_inf[2];
			%CF_Cluster=();
			%cluster_read=();
			%CG_Id=();
			%CG_Num=();
			load_cluster($map_inf[2]);
			@CFstart=sort{$a<=>$b}keys(%CF_Cluster);
			$chr_hash{$map_inf[2]}=1;
			$next_i=0;
			#print $map_inf[2]."\t".$CFstart[0]."\t".$CFstart[-1]."\n";
			
		}
		$map_inf[13]=~/XM:Z:(.*)$/;
		my $meth_str=$1;
		$map_inf[14]=~/XR:Z:(.*)$/;
		my $XR= $1;
		$map_inf[15]=~/XG:Z:(.*)$/;
		my $XG = $1;
		my $start=$map_inf[3];
		my $cgir=$map_inf[5];
		#print "cgir_methstr:".$cgir."\t".$meth_str."\n";
		my ($end,$rew_meth_str)=&read_mapping_length($cgir,$meth_str,$start,$XG);
		
		$end=$end+$start;
		#print "end_newmethstr:".$end."\t".$rew_meth_str."\n";
		my @pos=split /\t/,$rew_meth_str;
		if(@CFstart<1){
		#	print "No.... $map_inf[2]\n";
			next;
		}
		if($start>$CF_Cluster{$CFstart[-1]} or $end < $CFstart[0]){
			
		}else{
			for (my $cluster_i=$next_i;$cluster_i<@CFstart;$cluster_i++){
			
				if($start>$CF_Cluster{$CFstart[$cluster_i]}){
					if(exists($cluster_read{$cluster_i}) && $cluster_read{$cluster_i} ne ""){
					#	print $cluster_read{$cluster_i}."\n";
						&haplotype($cluster_read{$cluster_i},$CFstart[$cluster_i]);
						$cluster_read{$cluster_i}="";
					}
					$next_i=$cluster_i+1;
				}elsif($end < $CF_Cluster{$CFstart[$cluster_i]}){
					last;
				}else{
					#print "cluster_start_end:\t".$CFstart[$cluster_i]."\t".$CF_Cluster{$CFstart[$cluster_i]}."\n";
			
					
					my $nline="";
					
					for my $ps(@pos){
						
						my @spos=split /\//,$ps;
						#print $ps.":\t".$spos[0]."\t".$spos[1]."\n";
						if($spos[0]>=$CFstart[$cluster_i] && $spos[0]<=$CF_Cluster{$CFstart[$cluster_i]}){
							$nline.="\t".$ps
						}
						
					}
					if($nline ne ""){
						$nline=$map_inf[0].$nline."\n";
						if(!exists($cluster_read{$cluster_i})){
							$cluster_read{$cluster_i}=$nline;					
						}else{
							$cluster_read{$cluster_i}.=$nline;
						}
					}
					
				}
			}
		}
	}
	

}
if(exists($cluster_read{$CFstart[-1]}) && $cluster_read{$CFstart[-1]} ne ""){
	&haplotype($cluster_read{$CFstart[-1]},$CFstart[$CFstart[-1]]);
	$cluster_read{$CFstart[-1]}="";
}
close(SM);

sub haplotype{
	my ($cluster,$start)=@_;
	my $hap=&construct_haplotype($cluster);
	if($hap ne ""){
		print OUT $chr."\t".$CG_Id{$start}."\t".$CG_Num{$start}."\t".$hap;
	}
}

sub construct_haplotype{
	my $strand_read=shift;
	my @reads=split /\n/,$strand_read;
	
	%ploy_site=();
	%group=();
	%haplotype=();
	%heter=();
	for my $line(@reads){
		my @temp=split /\t/,$line;
		for(my $i=1;$i<@temp;$i++){
			my @temp1=split /\//,$temp[$i];
			if($temp1[0] ne ""){
				
				if(exists($heter{$temp1[0]}{$temp1[1]})){
					$heter{$temp1[0]}{$temp1[1]}++;
				}else{
					$heter{$temp1[0]}{$temp1[1]}=1;
				}

				if(exists($ploy_site{$temp1[0]})){
					$ploy_site{$temp1[0]}++;
				}else{
					$ploy_site{$temp1[0]}=1;
				}
				$group{$temp[0]}{$temp1[0]}=$temp1[1];
				#one bias: when the paired reads are overlaped and have conflicts on CpG calling. 
				#this version is not suitable for cfDNA
			}
		}
	}
	my $hap_str=&process_fragment();
	return($hap_str);
}

sub process_fragment{
	my @key_ploy_site=sort{$a<=>$b}keys(%ploy_site);	
	my @valid_site=();
	my %homo=();
	my %heter_sites=();
#	print "All sites ".@key_ploy_site.":\t @key_ploy_site\n";
	for my $site(@key_ploy_site){
		if($site ne ""){
			my @site_type=sort{${$heter{$site}}{$b}<=>${$heter{$site}}{$a}}keys(%{$heter{$site}});
			if(@site_type>1){  #control error, however, if  $ploy_site{$site}==1 it has been in homo
				if($heter{$site}{$site_type[0]}>=$heter_depth && $heter{$site}{$site_type[1]}>=$heter_depth){
					push @valid_site,$site;
				}elsif($heter{$site}{$site_type[0]}>=$homo_depth){
						$homo{$site}=$chr."\t".$site."\t".$site_type[0]."/".($heter{$site}{$site_type[0]}/2)."\t".$site_type[0]."/".($heter{$site}{$site_type[0]}/2)."\n";
				}else{
#					print "Discard 1:$site\t".($homo_mode_perc*($heter{$site}{$site_type[0]}+$heter{$site}{$site_type[1]}))."$site_type[0]\t$heter{$site}{$site_type[0]}\t$site_type[1]\t$heter{$site}{$site_type[1]}\n";
				}
			}elsif(@site_type==1 && $ploy_site{$site}>=$homo_depth){
				$homo{$site}=$chr."\t".$site."\t".$site_type[0]."/".($ploy_site{$site}/2)."\t".$site_type[0]."/".($ploy_site{$site}/2)."\n";
			}else{
#				print "Discard 2:$site\t$site_type[0]\t$heter{$site}{$site_type[0]}\t$site_type[1]\t$heter{$site}{$site_type[1]}\n";
			}
		 }else{
#			print "Discard 3:$site\n";
		 }
	}
	
	
	if(@valid_site==1){
		my @site_type=sort{${$heter{$valid_site[0]}}{$b}<=>${$heter{$valid_site[0]}}{$a}}keys(%{$heter{$valid_site[0]}});
		if($heter{$valid_site[0]}{$site_type[1]}>=$sec_to_max_perc*$heter{$valid_site[0]}{$site_type[0]}){
			$haplotype{$valid_site[0]}=$chr."\t".$valid_site[0]."\t".$site_type[0]."/".$heter{$valid_site[0]}{$site_type[0]}."\t".$site_type[1]."/".$heter{$valid_site[0]}{$site_type[1]}."\n";
					
		}else{
			$haplotype{$valid_site[0]}=$chr."\t".$valid_site[0]."\t".$site_type[0]."/".($heter{$valid_site[0]}{$site_type[0]}/2)."\t".$site_type[0]."/".($heter{$valid_site[0]}{$site_type[0]}/2)."\n";
				
		}
=cut		
		if($heter{$valid_site[0]}{$site_type[0]}>$heter{$valid_site[0]}{$site_type[1]}){
			$haplotype{$valid_site[0]}=$chr."\t".$valid_site[0]."\t".$site_type[0]."/".$heter{$valid_site[0]}{$site_type[0]}."\t".$site_type[1]."/".$heter{$valid_site[0]}{$site_type[1]}."\n";
					
		}elsif($neighboring_pattern{$max_p}>0.7*$mode_depth){
			$haplotype{$valid_site[0]}=$chr."\t".$valid_site[0]."\t".$site_type[1]."/".($heter{$valid_site[0]}{$site_type[1]}/2)."\t".$site_type[0]."/".($heter{$valid_site[0]}{$site_type[1]}/2)."\n";
				
		}
=cut		
	}elsif(@valid_site>=2){
		my @read_group=keys(%group);
		my %mode=();
		my $pattern="";
		my $mode_depth=0;
		for my $A(@read_group){
			$pattern="";
			for my $site(@valid_site){
				if(exists($group{$A}{$site})){
					$pattern.=$group{$A}{$site};
				}else{
					$pattern.=".";
				}
			}
			if(exists($mode{$pattern})){
				$mode{$pattern}++;
			}else{
				$mode{$pattern}=1;
			}
	#		print $pattern."\t$A\n";
		}
		
		my @kmode=sort{$mode{$b}<=>$mode{$a}}keys(%mode);
		my %neighboring_pattern=();
		for my $pattern(@kmode){
			my $sub_patternA=substr($pattern,0,2);
			if($sub_patternA=~/\./){
				#..,X.,.X should be excluded
			}else{
				$mode_depth+=$mode{$pattern};
				if(exists($neighboring_pattern{$sub_patternA})){
					$neighboring_pattern{$sub_patternA}+=$mode{$pattern};
				}else{
					$neighboring_pattern{$sub_patternA}=$mode{$pattern};
				}
			}
		}
		
		my @N_pattern=sort{$neighboring_pattern{$b}<=>$neighboring_pattern{$a}}keys(%neighboring_pattern);
		if(@N_pattern>=2){
			my @Two_site1=split //,$N_pattern[0];
			my @Two_site2=split //,$N_pattern[1];
			if($neighboring_pattern{$N_pattern[1]}>=$sec_to_max_perc*$neighboring_pattern{$N_pattern[0]}){
				$heter_sites{0}=$chr."\t".$valid_site[0]."\t".$Two_site1[0]."/".$neighboring_pattern{$N_pattern[0]}."\t".$Two_site2[0]."/".$neighboring_pattern{$N_pattern[1]}."\n";
				$heter_sites{1}=$chr."\t".$valid_site[1]."\t".$Two_site1[1]."/".$neighboring_pattern{$N_pattern[0]}."\t".$Two_site2[1]."/".$neighboring_pattern{$N_pattern[1]}."\n";
			}elsif($neighboring_pattern{$N_pattern[0]}>$homo_mode_perc*$mode_depth){
				$homo{$valid_site[0]}=$chr."\t".$valid_site[0]."\t".$Two_site1[0]."/".($neighboring_pattern{$N_pattern[0]}/2)."\t".$Two_site1[0]."/".($neighboring_pattern{$N_pattern[0]}/2)."\n";
				$homo{$valid_site[1]}=$chr."\t".$valid_site[1]."\t".$Two_site1[1]."/".($neighboring_pattern{$N_pattern[0]}/2)."\t".$Two_site1[1]."/".($neighboring_pattern{$N_pattern[0]}/2)."\n";
			}
		}elsif(@N_pattern==1){
			my @Two_site1=split //,$N_pattern[0];
			my @Two_site2=(1-$Two_site1[0],1-$Two_site1[1]);
			if($neighboring_pattern{$N_pattern[0]}>$homo_mode_perc*$heter{$valid_site[0]}{$Two_site1[0]} or $neighboring_pattern{$N_pattern[0]}>$homo_mode_perc*$heter{$valid_site[0]}{$Two_site1[1]}){
				$heter_sites{0}=$chr."\t".$valid_site[0]."\t".$Two_site1[0]."/".$neighboring_pattern{$N_pattern[0]}."\t".$Two_site2[0]."/0\n";
				$heter_sites{1}=$chr."\t".$valid_site[1]."\t".$Two_site1[1]."/".$neighboring_pattern{$N_pattern[0]}."\t".$Two_site2[1]."/0\n";
			}
		}
	
		my $former_site1="";
		my $former_site2="";
		my $break=1;
		my $former_real_rate=0;
		for(my $i=2;$i<@valid_site;$i++){
		#	print $valid_site[$i-1]."\t".$valid_site[$i]."\n";
			
			my $mode_depth=0;
			my %neighboring_pattern=();
			for my $pattern(@kmode){
				my $sub_patternA=substr($pattern,$i-1,2);
				if($sub_patternA=~/\./){
					#..,X.,.X should be excluded
				}else{
					$mode_depth+=$mode{$pattern};
					if(exists($neighboring_pattern{$sub_patternA})){
						$neighboring_pattern{$sub_patternA}+=$mode{$pattern};
					}else{
						$neighboring_pattern{$sub_patternA}=$mode{$pattern};
					}
				}
			}
			my @N_pattern=sort{$neighboring_pattern{$b}<=>$neighboring_pattern{$a}}keys(%neighboring_pattern);
			my $min_depth=($ploy_site{$valid_site[$i]}<$ploy_site{$valid_site[$i-1]}?$ploy_site{$valid_site[$i]}:$ploy_site{$valid_site[$i-1]});
			my $Expect_co_depth=($min_depth*(1-(($valid_site[$i]-$valid_site[$i-1])/150)));
			if($Expect_co_depth==0){
				$Expect_co_depth=1;
			}
			my $real_rate=$mode_depth/$Expect_co_depth;

			my $max_p="";		
			my $sec_max_p="";
			my @former_site=();
			my @count_1=();
			my @count_2=();
			if(exists($homo{$valid_site[$i-1]})){
				chomp $homo{$valid_site[$i-1]};
				@former_site=split /\t/,$homo{$valid_site[$i-1]};
				$homo{$valid_site[$i-1]}.="\n";
				@count_1=split /\//,$former_site[2];
				@count_2=split /\//,$former_site[3];
				$former_site1=$count_1[0];
				$former_site2=$count_2[0];
			}elsif(exists($heter_sites{$i-1})){
				chomp $heter_sites{$i-1};
				@former_site=split /\t/,$heter_sites{$i-1};
				$heter_sites{$i-1}.="\n";
				@count_1=split /\//,$former_site[2];
				@count_2=split /\//,$former_site[3];
				$former_site1=$count_1[0];
				$former_site2=$count_2[0];
			}else{
				$former_site1="";
				$former_site2="";
			}
			
		
			if(@N_pattern>=3 && $neighboring_pattern{$N_pattern[1]}==$neighboring_pattern{$N_pattern[2]}){
				my $max_1=0;
				for my $p(@N_pattern){
					if($neighboring_pattern{$p}>=2){
						my $first_s=substr($p,0,1);
						if( $max_p eq ""){
							$max_p=$p;
							if($former_site1 eq $first_s ){
								$max_1=1;
							}elsif($former_site2 eq $first_s){
								$max_1=2;
							}
						}elsif(($max_1==1 && $former_site2 eq $first_s) or ($max_1==2 && $former_site1 eq $first_s) && $sec_max_p eq ""){
							$sec_max_p=$p;
						}
						if($max_p ne "" && $sec_max_p ne ""){
							last;
						}
					}
				}
				if($max_p eq "" or $sec_max_p eq ""){
					if($neighboring_pattern{$N_pattern[0]}>=2){
						$max_p=$N_pattern[0];
					}
					if($neighboring_pattern{$N_pattern[1]}>=2){
						$sec_max_p=$N_pattern[1];
					}
				}
			}elsif(@N_pattern>=2){
				if($neighboring_pattern{$N_pattern[0]}>=2){
					$max_p=$N_pattern[0];
				}
				if($neighboring_pattern{$N_pattern[1]}>=2){
					$sec_max_p=$N_pattern[1];
				}
			}elsif(@N_pattern==1){
				if($neighboring_pattern{$N_pattern[0]}>=2){
					$max_p=$N_pattern[0];
				}
			}
		
			if($max_p ne "" && $sec_max_p ne ""){
				if($break==1){
					my @Two_site1=split //,$max_p;
					my @Two_site2=split //,$sec_max_p;
					#if($neighboring_pattern{$sec_max_p}>=$sec_to_max_perc*$neighboring_pattern{$max_p} && $neighboring_pattern{$max_p}+$neighboring_pattern{$sec_max_p}>= $heter_mode_perc*$mode_depth){
					
					if($neighboring_pattern{$sec_max_p}>=$sec_to_max_perc*$neighboring_pattern{$max_p}){
						if($Two_site1[0] ne $Two_site2[0]){
							$heter_sites{$i-1}=$chr."\t".$valid_site[$i-1]."\t".$Two_site1[0]."/".$neighboring_pattern{$max_p}."\t".$Two_site2[0]."/".$neighboring_pattern{$sec_max_p}."\n";
						
						}else{
							$homo{$valid_site[$i-1]}=$chr."\t".$valid_site[$i-1]."\t".$Two_site1[0]."/".$neighboring_pattern{$max_p}."\t".$Two_site2[0]."/".$neighboring_pattern{$sec_max_p}."\n";
						}
						if($Two_site1[1] ne $Two_site2[1]){
							$heter_sites{$i}=$chr."\t".$valid_site[$i]."\t".$Two_site1[1]."/".$neighboring_pattern{$max_p}."\t".$Two_site2[1]."/".$neighboring_pattern{$sec_max_p}."\n";
						
						}else{
							$homo{$valid_site[$i]}=$chr."\t".$valid_site[$i]."\t".$Two_site1[1]."/".$neighboring_pattern{$max_p}."\t".$Two_site2[1]."/".$neighboring_pattern{$sec_max_p}."\n";
						} 
						$break=0;
			#			print "1.1 ".$neighboring_pattern{$sec_max_p}."\t".$sec_to_max_perc."\t".$neighboring_pattern{$max_p}." $chr\t$valid_site[$i-1]\t$valid_site[$i]\t@N_pattern\n";
					}elsif($neighboring_pattern{$max_p}>$homo_mode_perc*$mode_depth){
						$homo{$valid_site[$i-1]}=$chr."\t".$valid_site[$i-1]."\t".$Two_site1[0]."/".($neighboring_pattern{$max_p}/2)."\t".$Two_site1[0]."/".($neighboring_pattern{$max_p}/2)."\n";
						$homo{$valid_site[$i]}=$chr."\t".$valid_site[$i]."\t".$Two_site1[1]."/".($neighboring_pattern{$max_p}/2)."\t".$Two_site1[1]."/".($neighboring_pattern{$max_p}/2)."\n";
						$break=0;
			#			print "1.2 $chr\t$valid_site[$i-1]\t$valid_site[$i]\t@N_pattern\n";
					}else{
			#			print "skip 1\n";
						$break=1;
					}
					
				}else{
					my $flag1=0;
					my $flag2=0;
					my $ts1=$max_p;
					my $ts2=$sec_max_p;
					
					my @Two_site1=split //, $max_p;
					my @Two_site2=split //, $sec_max_p;
			#		print "$former_site1 eq $Two_site1[0] && $former_site2 eq $Two_site2[0]\n";
					if($former_site1 eq $Two_site1[0] && $former_site2 eq $Two_site2[0]){
						$ts1=$max_p;
						$ts2=$sec_max_p;
						$flag1=1;
						$flag2=1;
					}elsif($former_site1 eq $Two_site2[0] && $former_site2 eq $Two_site1[0]){
						$ts2=$max_p;
						$ts1=$sec_max_p;
						$flag1=1;
						$flag2=1;
					}
					
					
					@Two_site1=();
					@Two_site2=();
					@Two_site1=split //,$ts1;
					@Two_site2=split //,$ts2;
					if($flag1==1 && $flag2==1){
						if($former_site1 eq $former_site2){
							#if($neighboring_pattern{$sec_max_p}>=$sec_to_max_perc*$neighboring_pattern{$max_p} && $neighboring_pattern{$max_p}+$neighboring_pattern{$sec_max_p}>= $heter_mode_perc*$mode_depth){
					
							if($neighboring_pattern{$sec_max_p}>=$sec_to_max_perc*$neighboring_pattern{$max_p}){
					
							#$haplotype{$valid_site[$i]}=$chr."\t".$valid_site[$i]."\t".$Two_site1[1]."/".$neighboring_pattern{$N_pattern[0]}."\t".$Two_site2[1]."/".$neighboring_pattern{$N_pattern[1]}."\n";
		#						print "2.1 break=1\n";
					
								if($Two_site1[1] ne $Two_site2[1]){
									$heter_sites{$i}=$chr."\t".$valid_site[$i]."\t".$Two_site1[1]."/".$neighboring_pattern{$ts1}."\t".$Two_site2[1]."/".$neighboring_pattern{$ts2}."\n";
								
								}
								$break=0;
							}elsif($neighboring_pattern{$max_p}>$homo_mode_perc*$mode_depth){
								$homo{$valid_site[$i]}=$chr."\t".$valid_site[$i]."\t".$Two_site1[1]."/".($neighboring_pattern{$max_p}/2)."\t".$Two_site1[1]."/".($neighboring_pattern{$max_p}/2)."\n";
		#						print "2.2 homo_continue\n";
								$break=0;
							}else{
				#				print "skip 2\n";
								$break=1;
							}
						}else{
							#if($neighboring_pattern{$sec_max_p}>=$sec_to_max_perc*$neighboring_pattern{$max_p} && $neighboring_pattern{$max_p}+$neighboring_pattern{$sec_max_p}>= $heter_mode_perc*$mode_depth){
					
							if($neighboring_pattern{$sec_max_p}>=$sec_to_max_perc*$neighboring_pattern{$max_p}){
								if($Two_site1[1] ne $Two_site2[1]){
									$heter_sites{$i}=$chr."\t".$valid_site[$i]."\t".$Two_site1[1]."/".$neighboring_pattern{$ts1}."\t".$Two_site2[1]."/".$neighboring_pattern{$ts2}."\n";
								}else{
									$homo{$valid_site[$i]}=$chr."\t".$valid_site[$i]."\t".$Two_site1[1]."/".($neighboring_pattern{$ts1}/2)."\t".$Two_site2[1]."/".($neighboring_pattern{$ts1}/2)."\n";
								}
								
								$break=0;
		#						print "3.1.1 $chr\t$valid_site[$i-1]\t$valid_site[$i]\t@N_pattern\n";
							}elsif($neighboring_pattern{$max_p}>$homo_mode_perc*$mode_depth){
								
							#	if($count_1[1]+$count_2[1]<$neighboring_pattern{$max_p}){
								if($former_real_rate<$real_rate){
									
									if(exists($heter_sites{$i-1})){
				#						print "delete $heter_sites{$i-1}";
										delete $heter_sites{$i-1};
									}elsif($homo{$valid_site[$i-1]}){
				#						print "delete $homo{$valid_site[$i-1]}";
										delete $homo{$valid_site[$i-1]};
									}
									$i--;
									
				#					print "3.2.2 homo break $chr\t$valid_site[$i-1]\t$valid_site[$i]\t$N_pattern[0]"."/".$count_1[1]."\t$N_pattern[1]/".$count_2[1]."\n";
				#					print "3.2.1 $chr\t$valid_site[$i-1]\t$valid_site[$i]\t$N_pattern[0]"."/".$neighboring_pattern{$N_pattern[0]}."\t$N_pattern[1]/".$neighboring_pattern{$N_pattern[1]}."\n";
				#					print "-------".$haplotype{$valid_site[$i-1]};
								}elsif($i==@valid_site-1){
									@Two_site1=split //, $max_p;
									$homo{$valid_site[$i]}=$chr."\t".$valid_site[$i]."\t".$Two_site1[1]."/".($neighboring_pattern{$max_p}/2)."\t".$Two_site1[1]."/".($neighboring_pattern{$max_p}/2)."\n";
								}
								$break=1;
				#				print "back 1 \n";
		#						print "3.1.2 $chr\t$valid_site[$i-1]\t$valid_site[$i]\t@N_pattern\n";
							}else{
				#				print "skip 3\n";
								$break=1;
							}
							
						}
					}else{
						
					#	if(($neighboring_pattern{$sec_max_p}>=$sec_to_max_perc*$neighboring_pattern{$max_p} && $count_1[1]+$count_2[1]<$neighboring_pattern{$sec_max_p}+$neighboring_pattern{$max_p}) or ($neighboring_pattern{$max_p}>$homo_mode_perc*$mode_depth  && $count_1[1]+$count_2[1]<$neighboring_pattern{$max_p})){
						if(($neighboring_pattern{$sec_max_p}>=$sec_to_max_perc*$neighboring_pattern{$max_p} or $neighboring_pattern{$max_p}>$homo_mode_perc*$mode_depth ) && $former_real_rate<$real_rate){
							if(exists($heter_sites{$i-1})){
								delete $heter_sites{$i-1};
							}elsif($homo{$valid_site[$i-1]}){
								delete $homo{$valid_site[$i-1]};
							}
							
							$i--;
		#					print "3.2.2 homo break $chr\t$valid_site[$i-1]\t$valid_site[$i]\t$N_pattern[0]"."/".$count_1[1]."\t$N_pattern[1]/".$count_2[1]."\n";
		#					print "3.2.1 $chr\t$valid_site[$i-1]\t$valid_site[$i]\t$N_pattern[0]"."/".$neighboring_pattern{$N_pattern[0]}."\t$N_pattern[1]/".$neighboring_pattern{$N_pattern[1]}."\n";
		#					print "-------".$haplotype{$valid_site[$i-1]};
				#			print "back 2 \n";
						}elsif($i==@valid_site-1){
							if($neighboring_pattern{$sec_max_p}>=$sec_to_max_perc*$neighboring_pattern{$max_p}){
								if($Two_site1[1] ne $Two_site2[1]){
									$heter_sites{$i}=$chr."\t".$valid_site[$i]."\t".$Two_site1[1]."/".$neighboring_pattern{$ts1}."\t".$Two_site2[1]."/".$neighboring_pattern{$ts2}."\n";
								}else{
									
									$homo{$valid_site[$i]}=$chr."\t".$valid_site[$i]."\t".$Two_site1[1]."/".$neighboring_pattern{$ts1}."\t".$Two_site2[1]."/".$neighboring_pattern{$ts2}."\n";
								}
							}elsif($neighboring_pattern{$max_p}>$homo_mode_perc*$mode_depth){
								@Two_site1=split //, $max_p;
								$homo{$valid_site[$i]}=$chr."\t".$valid_site[$i]."\t".$Two_site1[1]."/".($neighboring_pattern{$max_p}/2)."\t".$Two_site1[1]."/".($neighboring_pattern{$max_p}/2)."\n";
							}
						}else{
				#			print "skip 4: $former_site1 ($flag1)\t$former_site2 ($flag2);\n";
						}
						$break=1;
					}
					
				}
			}elsif($max_p ne ""){
				my @Two_site1=split //,$max_p;
				if($neighboring_pattern{$max_p}>$homo_mode_perc*$mode_depth){
					if($break==1){
						$homo{$valid_site[$i-1]}=$chr."\t".$valid_site[$i-1]."\t".$Two_site1[0]."/".($neighboring_pattern{$max_p}/2)."\t".$Two_site1[0]."/".($neighboring_pattern{$max_p}/2)."\n";
						$homo{$valid_site[$i]}=$chr."\t".$valid_site[$i]."\t".$Two_site1[1]."/".($neighboring_pattern{$max_p}/2)."\t".$Two_site1[1]."/".($neighboring_pattern{$max_p}/2)."\n";
						$break=0;
					}else{
						
						if($former_site1 eq $former_site2 && $former_site1 eq $Two_site1[0]){
							$homo{$valid_site[$i]}=$chr."\t".$valid_site[$i]."\t".$Two_site1[1]."/".($neighboring_pattern{$max_p}/2)."\t".$Two_site1[1]."/".($neighboring_pattern{$max_p}/2)."\n";
							$break=0;
		#					print "3.2.1 homo continue $chr\t$valid_site[$i-1]\t$valid_site[$i]\t$N_pattern[0]"."/".$neighboring_pattern{$N_pattern[0]}."\t$N_pattern[1]/".$neighboring_pattern{$N_pattern[1]}."\n";
								
					#	}elsif(($neighboring_pattern{$sec_max_p}>=$sec_to_max_perc*$neighboring_pattern{$max_p} && $count_1[1]+$count_2[1]<$neighboring_pattern{$sec_max_p}+$neighboring_pattern{$max_p}) or ($neighboring_pattern{$max_p}>$homo_mode_perc*$mode_depth  && $count_1[1]+$count_2[1]<$neighboring_pattern{$max_p})){
						}elsif($former_real_rate<$real_rate){
						
							if(exists($heter_sites{$i-1})){
				#					print "Delete....$heter_sites{$i-1}";
									delete $heter_sites{$i-1};
				#					print "Delete....$heter_sites{$i-1}";
								}elsif($homo{$valid_site[$i-1]}){
				#					print "Delete....$homo{$valid_site[$i-1]}";
									delete $homo{$valid_site[$i-1]};
				#					print "Delete....$homo{$valid_site[$i-1]}";
								}
								
								$i--;
								$break=1;
				#				print "back 3 \n";
		#						print "3.2.2 homo break $chr\t$valid_site[$i-1]\t$valid_site[$i]\t$N_pattern[0]"."/".$count_1[1]."\t$N_pattern[1]/".$count_2[1]."\n";
		#						print "-------".$haplotype{$valid_site[$i-1]};
						}elsif($i==@valid_site-1){
							if($neighboring_pattern{$max_p}>$homo_mode_perc*$mode_depth){
								
								$homo{$valid_site[$i]}=$chr."\t".$valid_site[$i]."\t".$Two_site1[1]."/".($neighboring_pattern{$max_p}/2)."\t".$Two_site1[1]."/".($neighboring_pattern{$max_p}/2)."\n";
							
							}
						}else{
				#			print "ignore:........$former_real_rate > $real_rate\n";
							$break=1;
						}
					}
				}else{
				#	print "skip 5\n";
					$break=1;
				}
			}
		}
		
		my @homo_site=sort{$a<=>$b}keys(%homo);
		my @hap_site=sort{$a<=>$b}keys(%heter_sites);
		
		#print "The number of heter_sites ".@valid_site.": @valid_site\n";		
	#	print "All homo".@homo_site.":@homo_site\n";
		#print "The heter_sites in haplotype: ".@hap_site.":@hap_site\n";
		$break=0;
		if(@hap_site !=@valid_site){
			if(@hap_site==1){
				chomp $heter_sites{$hap_site[0]};
				$haplotype{$valid_site[$hap_site[0]]}=$heter_sites{$hap_site[0]}."\tT1\n";
			}elsif(@hap_site>1){
				for(my $i=1;$i<@hap_site;$i++){
				#	print $valid_site[$hap_site[$i-1]]."\t".$valid_site[$hap_site[$i]]."\n";
					
					my $mode_depth=0;
					my %neighboring_pattern=();
					for my $pattern(@kmode){
						my $sub_patternA=substr($pattern,$hap_site[$i-1],1);
						$sub_patternA.=substr($pattern,$hap_site[$i],1);
						if($sub_patternA=~/\./){
							#..,X.,.X should be excluded
						}else{
							$mode_depth+=$mode{$pattern};
							if(exists($neighboring_pattern{$sub_patternA})){
								$neighboring_pattern{$sub_patternA}+=$mode{$pattern};
							}else{
								$neighboring_pattern{$sub_patternA}=$mode{$pattern};
							}
						}
					}
					my @N_pattern=sort{$neighboring_pattern{$b}<=>$neighboring_pattern{$a}}keys(%neighboring_pattern);

					my @former_site=();
					if(exists($haplotype{$valid_site[$hap_site[$i-1]]})){
						chomp $haplotype{$valid_site[$hap_site[$i-1]]};
						@former_site=split /\t/,$haplotype{$valid_site[$hap_site[$i-1]]};
						$haplotype{$valid_site[$hap_site[$i-1]]}.="\n";
					}else{
						chomp $heter_sites{$hap_site[$i-1]};
						@former_site=split /\t/,$heter_sites{$hap_site[$i-1]};
					}
					
					my @count_11=split /\//,$former_site[2];
					my @count_21=split /\//,$former_site[3];
					my @typeA=($count_11[0],$count_21[0]);
					chomp $heter_sites{$hap_site[$i]};
					my @former_site=split /\t/,$heter_sites{$hap_site[$i]};
					my @count_12=split /\//,$former_site[2];
					my @count_22=split /\//,$former_site[3];
					my @typeB=($count_12[0],$count_22[0]);
					if($mode_depth>1){
						
						
						my $AC=$typeA[0].$typeB[0];
						my $BD=$typeA[1].$typeB[1];
						my $f1=0;
						my $AD=$typeA[0].$typeB[1];
						my $BC=$typeA[1].$typeB[0];
						my $f2=0;
						
						if(exists($neighboring_pattern{$AC})){
							$f1=$neighboring_pattern{$AC};
						}else{
							$neighboring_pattern{$AC}=0;
						}
						if(exists($neighboring_pattern{$BD})){
							$f1+=$neighboring_pattern{$BD};
						}else{
							$neighboring_pattern{$BD}=0;
						}
						if(exists($neighboring_pattern{$AD})){
							$f2=$neighboring_pattern{$AD};
						}else{
							$neighboring_pattern{$AD}=0;
						}
						if(exists($neighboring_pattern{$BC})){
							$f2+=$neighboring_pattern{$BC};
						}else{
							$neighboring_pattern{$BC}=0;
						}
			#			print "@typeA\t@typeB\t".$AC."/".$neighboring_pattern{$AC}."\t".$BD."/".$neighboring_pattern{$BD}."\t".$AD."/".$neighboring_pattern{$AD}."\t".$BC."/".$neighboring_pattern{$BC}."\n";
						if(exists($haplotype{$valid_site[$hap_site[$i-1]]})){
							
							if($f1>$f2){
								$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[0]."/".$neighboring_pattern{$AC}."\t".$typeB[1]."/".$neighboring_pattern{$BD}."\tT21\n";
			#					print "E12\n";
							}elsif($f1<$f2){
								$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[1]."/".$neighboring_pattern{$AD}."\t".$typeB[0]."/".$neighboring_pattern{$BC}."\tT21\n";
			#					print "E21\n";
							}elsif($f1==$f2){
								if($f1>0){
			#						print "E00\n";
									if($N_pattern[0] eq $AC or $N_pattern[0] eq $BD){
										$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[0]."/".$neighboring_pattern{$AC}."\t".$typeB[1]."/".$neighboring_pattern{$BD}."\tT221\n";
									}elsif($N_pattern[0] eq $AD or $N_pattern[0] eq $BC){
										$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[1]."/".$neighboring_pattern{$AD}."\t".$typeB[0]."/".$neighboring_pattern{$BC}."\tT222\n";
									}
								}else{
			#						print "E01\n";
									if($i==@hap_site-1){
										$haplotype{$valid_site[$hap_site[$i]]}=$heter_sites{$hap_site[$i]}."\tT23\n";
									}
									$break=1;
								}
							}
						}else{
							if($f1>$f2){
			#					print "NE12\n";
								if($count_11[1]+$count_21[1]>=$count_12[1]+$count_22[1]){
									if($count_11[1]>=$count_21[1]){
										$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[0]."/".$neighboring_pattern{$AC}."\t".$typeA[1]."/".$neighboring_pattern{$BD}."\tT33\n";
										$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[0]."/".$neighboring_pattern{$AC}."\t".$typeB[1]."/".$neighboring_pattern{$BD}."\tT33\n";
									}else{
										$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[1]."/".$neighboring_pattern{$BD}."\t".$typeA[0]."/".$neighboring_pattern{$AC}."\tT33\n";
										$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[1]."/".$neighboring_pattern{$BD}."\t".$typeB[0]."/".$neighboring_pattern{$AC}."\tT33\n";
									}
								}else{
									if($count_12[1]>=$count_22[1]){
										$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[0]."/".$neighboring_pattern{$AC}."\t".$typeA[1]."/".$neighboring_pattern{$BD}."\tT34\n";
										$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[0]."/".$neighboring_pattern{$AC}."\t".$typeB[1]."/".$neighboring_pattern{$BD}."\tT34\n";
									}else{
										$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[1]."/".$neighboring_pattern{$BD}."\t".$typeA[0]."/".$neighboring_pattern{$AC}."\tT34\n";
										$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[1]."/".$neighboring_pattern{$BD}."\t".$typeB[0]."/".$neighboring_pattern{$AC}."\tT34\n";
							
									}
								}	
							}elsif($f1<$f2){
			#					print "NE21\n";
								if($count_11[1]+$count_21[1]>=$count_12[1]+$count_22[1]){
									if($count_11[1]>=$count_21[1]){
										$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[0]."/".$neighboring_pattern{$AD}."\t".$typeA[1]."/".$neighboring_pattern{$BC}."\tT35\n";
										$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[0]."/".$neighboring_pattern{$AD}."\t".$typeB[1]."/".$neighboring_pattern{$BC}."\tT35\n";
									}else{
										$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[1]."/".$neighboring_pattern{$BC}."\t".$typeA[0]."/".$neighboring_pattern{$AD}."\tT35\n";
										$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[1]."/".$neighboring_pattern{$BC}."\t".$typeB[0]."/".$neighboring_pattern{$AD}."\tT35\n";
							
									}
								}else{
									if($count_12[1]>=$count_22[1]){
										$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[0]."/".$neighboring_pattern{$AD}."\t".$typeA[1]."/".$neighboring_pattern{$BC}."\tT36\n";
										$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[0]."/".$neighboring_pattern{$AD}."\t".$typeB[1]."/".$neighboring_pattern{$BC}."\tT36\n";
									}else{
										$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[1]."/".$neighboring_pattern{$BC}."\t".$typeA[0]."/".$neighboring_pattern{$AD}."\tT36\n";
										$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[1]."/".$neighboring_pattern{$BC}."\t".$typeB[0]."/".$neighboring_pattern{$AD}."\tT36\n";
							
									}
								}
							}elsif($f1==$f2){
				#				print "NE00\n";
								if($f1!=0){
									if($N_pattern[0] eq $AC or $N_pattern[0] eq $BD){
										if($count_11[1]+$count_21[1]>$count_12[1]+$count_22[1]){
											if($count_11[1]>$count_21[1]){
												$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[0]."/".$neighboring_pattern{$AC}."\t".$typeA[1]."/".$neighboring_pattern{$BD}."\tT37\n";
												$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[0]."/".$neighboring_pattern{$AC}."\t".$typeB[1]."/".$neighboring_pattern{$BD}."\n";
											}else{
												$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[1]."/".$neighboring_pattern{$BD}."\t".$typeA[0]."/".$neighboring_pattern{$AC}."\tT37\n";
												$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[1]."/".$neighboring_pattern{$BD}."\t".$typeB[0]."/".$neighboring_pattern{$AC}."\tT37\n";
									
											}
										}else{
											if($count_12[1]>$count_22[1]){
												$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[0]."/".$neighboring_pattern{$AC}."\t".$typeA[1]."/".$neighboring_pattern{$BD}."\tT38\n";
												$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[0]."/".$neighboring_pattern{$AC}."\t".$typeB[1]."/".$neighboring_pattern{$BD}."\tT38\n";
											}else{
												$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[1]."/".$neighboring_pattern{$BD}."\t".$typeA[0]."/".$neighboring_pattern{$AC}."\tT38\n";
												$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[1]."/".$neighboring_pattern{$BD}."\t".$typeB[0]."/".$neighboring_pattern{$AC}."\tT38\n";
									
											}
										}
									}elsif($N_pattern[0] eq $AD or $N_pattern[0] eq $BC){
										if($count_11[1]+$count_21[1]>$count_12[1]+$count_22[1]){
											if($count_11[1]>$count_21[1]){
												$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[0]."/".$neighboring_pattern{$AD}."\t".$typeA[1]."/".$neighboring_pattern{$BC}."\tT39\n";
												$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[0]."/".$neighboring_pattern{$AD}."\t".$typeB[1]."/".$neighboring_pattern{$BC}."\tT39\n";
											}else{
												$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[1]."/".$neighboring_pattern{$BC}."\t".$typeA[0]."/".$neighboring_pattern{$AD}."\tT39\n";
												$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[1]."/".$neighboring_pattern{$BC}."\t".$typeB[0]."/".$neighboring_pattern{$AD}."\tT39\n";
									
											}
										}else{
											if($count_12[1]>$count_22[1]){
												$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[0]."/".$neighboring_pattern{$AD}."\t".$typeA[1]."/".$neighboring_pattern{$BC}."\tT301\n";
												$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[0]."/".$neighboring_pattern{$AD}."\t".$typeB[1]."/".$neighboring_pattern{$BC}."\tT301\n";
											}else{
												$haplotype{$valid_site[$hap_site[$i-1]]}=$chr."\t".$valid_site[$hap_site[$i-1]]."\t".$typeA[1]."/".$neighboring_pattern{$BC}."\t".$typeA[0]."/".$neighboring_pattern{$AD}."\tT301\n";
												$haplotype{$valid_site[$hap_site[$i]]}=$chr."\t".$valid_site[$hap_site[$i]]."\t".$typeB[1]."/".$neighboring_pattern{$BC}."\t".$typeB[0]."/".$neighboring_pattern{$AD}."\tT301\n";
									
											}
										}
									}
								}else{
				#					print "NE01\n";
									$haplotype{$valid_site[$hap_site[$i-1]]}=$heter_sites{$hap_site[$i-1]}."\tT4\n";
									if($i==@hap_site-1){
										$haplotype{$valid_site[$hap_site[$i]]}=$heter_sites{$hap_site[$i]}."\tT4\n";
									}
									$break=1;
								}
							}
						
						}	
					}else{
				#		print "C00\n";
						if(!exists($haplotype{$valid_site[$hap_site[$i-1]]})){
							$haplotype{$valid_site[$hap_site[$i-1]]}=$heter_sites{$hap_site[$i-1]}."\tT5\n";
						}
						if($i==@hap_site-1){
							$haplotype{$valid_site[$hap_site[$i]]}=$heter_sites{$hap_site[$i]}."\tT5\n";
						}
						$break=1;
					}

				}	
			}
		}else{
			for(my $i=1;$i<@hap_site;$i++){
				$haplotype{$valid_site[$hap_site[$i]]}=$heter_sites{$hap_site[$i]}."\tT5\n";
				
			}
			
		}
	}
	
		
	my @hsite=sort{$a<=>$b}keys(%haplotype);
#	print "The heter_sites in haplotype: ".@hsite.":@hsite\n";
	
	my $re_string="";
	my @homo_site=sort{$a<=>$b}keys(%homo);
	if(@hsite>0 or @homo_site >0){
	#if(@hsite>0){
	#print OUT "break\n";
		
		for my $hs(@homo_site){
			if(exists($haplotype{$hs})){
#				print "Error:$haplotype{$hs}";
#				print $homo{$hs};
			}else{
				
				my @former_site=split /\t/,$homo{$hs};
				my @count_11=split /\//,$former_site[2];
				my @count_21=split /\//,$former_site[3];
				$haplotype{$hs}=$former_site[0]."\t".$former_site[1]."\t".$count_11[0]."/".($heter{$hs}{$count_11[0]}/2)."\t".$count_11[0]."/".($heter{$hs}{$count_11[0]}/2)."\n"; 
				delete $homo{$hs}; 
			 }
		}

		@hsite=sort{$a<=>$b}keys(%haplotype);
		my $hap1="";
		my $hap1_c="";
		my $hap2="";
		my $hap2_c="";
		for my $site(@hsite){
#			print "Final_print:\t$site\t".$haplotype{$site};
			chomp $haplotype{$site};
			my @former_site=split /\t/,$haplotype{$site};
			my @count_11=split /\//,$former_site[2];
			my @count_21=split /\//,$former_site[3];
			$hap1.=$count_11[0];
			if($hap1_c ne ""){
				$hap1_c.=",".$heter{$site}{$count_11[0]};
			}else{
				$hap1_c=$heter{$site}{$count_11[0]}
			}
			$hap2.=$count_21[0];
			if($hap2_c ne ""){
				$hap2_c.=",".$heter{$site}{$count_21[0]};
			}else{
				$hap2_c=$heter{$site}{$count_21[0]};
			}
			
		#	$haplotype{$site}=$former_site[0]."\t".$former_site[1]."\t".$count_11[0]."/".$heter{$site}{$count_11[0]}."\t".$count_21[0]."/".$heter{$site}{$count_21[0]}."\n"; 
		#	print OUT $haplotype{$site};
		}
		if(@hsite>=3){
			$re_string=@hsite."\t".$hsite[0]."\t".$hsite[-1]."\t".$hap1."\t".$hap2."\t".$hap1_c."\t".$hap2_c."\n";
		}
	}
#	print OUT "break2\n";
	%haplotype=();
	return($re_string);
}

sub load_cluster{
	my $chr=shift;
	open(INC,$CF_file)or die $!;
	while(my $line=<INC>){
		chomp $line;
		my @temp=split /\t/,$line;
		if($temp[1] eq $chr){
			$CF_Cluster{$temp[2]}=$temp[3];
			$CG_Num{$temp[2]}=$temp[5];
			$CG_Id{$temp[2]}=$temp[0];
		}
	}
	close(INC);
}

sub read_mapping_length{
	my ($cgir,$meth_str,$start,$XG)=@_;
	my $len=0;
	my $ref_meth_st="";
	my @meth_array=split //,$meth_str;
	my $array_start=0;
	if($cgir=~/[I|D|M]/){			
		my @CGIR=split //,$cgir;
		my $num=0;
		for(my $i=0;$i<@CGIR;$i++){
			
			if($CGIR[$i] =~/\d/){
				$num=$CGIR[$i];
				$i++;
				while($CGIR[$i] =~/\d/){
					$num=$num*10+$CGIR[$i];
					$i++;
				}
				$i--;
			}elsif($CGIR[$i] =~/M/){
				$len+=$num;
				for (my $j=$array_start;$j<$array_start+$num;$j++){
					$ref_meth_st=$ref_meth_st.$meth_array[$j];
				}
				$array_start=$array_start+$num;
				$num=0;
			}elsif($CGIR[$i] =~/D/){
				$len+=$num;
				for (my $j=$array_start;$j<$array_start+$num;$j++){
					$ref_meth_st=$ref_meth_st.".";
				}
				$num=0;
			}elsif($CGIR[$i] =~/I/){
				$array_start=$array_start+$num;
				$num=0;
			}
		}
	}
	my @meth_state=split //,$ref_meth_st;
	my $read_meth="";
	my $pos=0;
	for (my $ms=0;$ms<@meth_state;$ms++){
		if($meth_state[$ms] ne "."){ #note: for CpG, only Z and z should be counted.
										#for CA,CC,CT, more details should be considered 
				if($meth_state[$ms]=~/[Z]/){
					if($XG eq "GA"){
						$pos=$ms+$start-1;
					}else{
						$pos=$ms+$start;
					}
					if($read_meth ne ""){
						$read_meth.="\t$pos/1";
					}else{
						$read_meth="$pos/1";
					}
				}elsif($meth_state[$ms]=~/[z]/){
					if($XG eq "GA"){
						$pos=$ms+$start-1;
					}else{
						$pos=$ms+$start;
					}
					if($read_meth ne ""){
						$read_meth.="\t$pos/0";
					}else{
						$read_meth="$pos/0";
					}
				}
		}
	}
	return ($len,$read_meth);
}
close(OUT);
close();




