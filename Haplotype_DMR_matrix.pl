=cut
	1. contruct methylation haplotype matrix
	2. construct statistic metric and calculate for each element region
	3. select a random proportion of regions, save their statistic values
	4. predict a DMR and assign a pvalue by comparing with the sample set of step 3.
=cut


use strict;

my $Hapsample_group_inf=shift;
my $pro=shift;
open(LI,$Hapsample_group_inf)or die $!;

my %group_sample=();
my @file_list=();
my @group_num=();
while(my $line=<LI>){
	chomp $line;
	my @temp=split /\s/,$line;
	$group_sample{$temp[1]}{$temp[0]}=1;
}
close(LI);
my @groups=keys(%group_sample);

for my $gs(@groups){
	my @samples=keys(%{$group_sample{$gs}});
	push @file_list,@samples;
	my $num=@samples;
	push @group_num,$num;
}
my @statistic_value=();
my $Matrix=&Hap_matrix();
my @random_list=();
&random($pro);
&DMR_Ex($Matrix);

sub DMR_Ex{
	my $hap_matrix=shift;
	my %DMR=();
	my %DMR_pvalue=();
	my $id=0;
	open(IN, $hap_matrix)or die $!;
	while(my $line=<IN>){
		chomp $line;
		my @temp=split /\t/,$line;
		my @pattern=();
		$temp[5]=~/^(.*)\/(.*)\/(.*)\/(.*)$/;
		push @pattern,$1;
		$temp[6]=~/^(.*)\/(.*)\/(.*)\/(.*)$/;
		push @pattern,$1;
		if($pattern[0] ne $pattern[1]){
			
			my $pvalue=&calculate_pvalue($temp[-1]);
			#print $pattern[0]."\t".$pattern[1]."\t".$temp[-1]."\t".$pvalue."\n";
			if($pvalue<0.05 && $temp[-1]>0.5){
				$DMR{$id}="$line\t$pvalue\n";
				$DMR_pvalue{$id}=$pvalue;
				$id++;
			}
		}
	}
	close(IN);
	my $out=$pro."_DMR_pvalue.txt";
	open(OUT,">$out")or die $!;
	
	my @DMR_id=sort{$DMR_pvalue{$a}<=>$DMR_pvalue{$b}}keys(%DMR_pvalue);
	print @DMR_id."\n";
	for my $dmr_id(@DMR_id){
		print OUT $DMR{$dmr_id};
	}
	close(OUT);
}

sub calculate_pvalue{
	my $statis_value=shift;
	my $num=0;
	for my $sta(@random_list){
		if($statis_value<$sta){
			$num++;
		}
	}
	my $pvalue=$num/@random_list;
	return $pvalue;
}

sub random{
	my $proportion=shift;
	
	my $num=@statistic_value*$proportion;
	srand();
	for my $sta(@statistic_value){
		my $select=rand(1);
		if($select>0.5){
			push @random_list,$sta;
			$num--;
		}
		if($num<0){
			last;
		}
	}
}


sub Hap_matrix{
	my %cluster=();
	my %cluster_id_inf=();
	my @type_cls=("LL","HH","LN","HN","MN","NN","LM","HM","MM","HL");
	my $out="Hap_matrix.txt";
	my $sample_i=0;			
	foreach my $file(@file_list){
		$sample_i++;
		open(LI,$file)or die $!;
		while(my $line=<LI>){
			chomp $line;
			my @temp=split /\t/, $line;
			my $class="";
			my $reliab1=0;
			my $Lab1="";
			my $Lab2="";
			
			$cluster_id_inf{$temp[0]}{$temp[1]}=$temp[2]."\t".$temp[4]."\t".$temp[5];
				
			if($temp[3]>3){
				$reliab1=$temp[3]/$temp[2];			
				my ($L1,$mr1)=&hap_xor($temp[6],0,$temp[8]);
				my ($L2,$mr2)=&hap_xor($temp[7],0,$temp[9]);
				#$L1=$L1/$temp[3];
				#$L2=$L2/$temp[3];
				
				$L1=$mr1;
				$L2=$mr2;
				if($L1<=0.25){
					$Lab1="L";
					
				}elsif($L1<=0.5){
					$Lab1="N";	
					
				}elsif($L1<=0.75){
					$Lab1="M";	
					
				}else{
					$Lab1="H";	
					
				}
				if($L2<=0.25){
					$Lab2="L";
					
				}elsif($L2<=0.5){
					$Lab2="N";	
					
				}elsif($L2<=0.75){
					$Lab2="M";	
					
				}else{
					$Lab2="H";	
					
				}
				my @Labs=($Lab1,$Lab2);
				@Labs=sort(@Labs);
				my $Two_lab=join("",@Labs);
				$class=$Two_lab."/$reliab1/$L1/$L2";

			}else{
				$class=".";
			
			}
			if(!exists($cluster{$temp[0]}{$temp[1]})){
				
				if($sample_i!=1){
					$cluster{$temp[0]}{$temp[1]}=".";
					for (my $k=1;$k<$sample_i-1;$k++){
						$cluster{$temp[0]}{$temp[1]}.="\t.";
					}
					$cluster{$temp[0]}{$temp[1]}.="\t".$class;
				}else{
					$cluster{$temp[0]}{$temp[1]}=$class;
				}
			}else{
				my @field=split /\t/,$cluster{$temp[0]}{$temp[1]};
				if($sample_i>@field+1){
					for (my $k=@field;$k<$sample_i-1;$k++){
						$cluster{$temp[0]}{$temp[1]}.="\t.";
					}
					$cluster{$temp[0]}{$temp[1]}.="\t".$class;
				}else{
					$cluster{$temp[0]}{$temp[1]}.="\t".$class;
				}
			}
			
		}
		close(LI);
	}

	my @chr_id=keys(%cluster);
	open(OUT,">$out")or die $!;
	print OUT "chr\tclusterid\tcpgs\tstart\tend\t@file_list\n";
	for my $chr(@chr_id){		
		my @cluster_id=sort{$a<=>$b}keys(%{$cluster{$chr}}); 
		for my $id(@cluster_id){
			my @field=split /\t/,$cluster{$chr}{$id};
			if(@file_list==@field){
				if($field[0] ne "." && $field[1] ne "."){
					my $max1=0;					
					my $min1=0;
					my $min2=0;
					my $max2=0;
					my $reb1=0;
					my $reb2=0;
				
					if($field[0]=~/^(.*)\/(.*)\/(.*)\/(.*)$/){
						$reb1=$2;
						if($3<$4){
							$max1=$4;
							$min1=$3;
						}else{
							$max1=$3;
							$min1=$4;
						}
					}
					
					if($field[1]=~/^(.*)\/(.*)\/(.*)\/(.*)$/){
						$reb2=$2;
						if($3<$4){
							$max2=$4;
							$min2=$3;
						}else{
							$max2=$3;
							$min2=$4;
						}
					}
					my $weight=abs($max1-$max2)>abs($min1-$min2)?abs($max1-$max2):abs($min1-$min2);
					#my $reb=($reb1<$reb2)?$reb1:$reb2;
					#print $reb."\n";
					#$weight=$weight*$reb;
					push @statistic_value,$weight;
					print OUT $chr."\t".$id."\t".$cluster_id_inf{$chr}{$id}."\t".$cluster{$chr}{$id}."\t".$weight."\n";
				}
			}
		}
		
	}
	close(OUT);
	return ($out);
}

sub hap_xor{
	my ($hap_str,$binary,$depth)=@_;
	my @temp=split //,$hap_str;
	my @deps=split /,/,$depth;
	@deps=sort{$a<=>$b}(@deps);
	my $mid=0;
	my $max_depth=$deps[-1];
	if(@deps%2==0){
		$mid=($deps[@deps/2-1]+$deps[@deps/2])/2;
	}else{
		$mid=$deps[@deps/2];
	}
	my $num=0;
	my $mr=0;
	my $ur=0;
	my $i=0;
	for my $bit(@temp){
		if($bit != $binary){
			$num++;
			$mr+=$deps[$i];
		}else{
			$ur+=$deps[$i]
		}
		$i++;
	}
	if($mr!=0){
		$mr=$mr/($ur+$mr);
	}
	return ($num,$mr);
}
close(ST);


