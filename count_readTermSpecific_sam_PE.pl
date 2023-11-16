#!/usr/bin/perl -w

## version v3

if ((@ARGV == 0) || ($ARGV[0]=~/-h/)) {
	print "OPTION:
     ## Both
        -b samfile
     ## Fwd & Rev
        -f fwd_samfile (-r Rev_samfile )
      \n";
#	exit;

}

foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};

$adapt= exists $option->{-adapt} ? $option->{-adapt}->[0] : "upm";

$terminfo= {};
if ( exists $option->{-b}) {
    $inform= "both";
    ($terminfo)= &create_terminfo($option->{-b}->[0],$inform,$terminfo);
} elsif (exists $option->{-f} || exists $option->{-r}) {
    $inform= "fwd";
    ($terminfo)= &create_terminfo($option->{-f}->[0],$inform,$terminfo,$adapt);
    $inform= "rev";
    ($terminfo)= &create_terminfo($option->{-r}->[0],$inform,$terminfo,$adapt);

}

print "### FINISHED READ SAM\n";

open (LOG2,">out_countTerm_outlier.txt");


($cno)= &output_count_term($terminfo,$adapt);

print "$cno\n";
print "### FINISHED CAL OUTLIER\n";

##=== === === === === === === === === === === === === === === === ===
sub create_terminfo {
	my ($infile,$inform,$terminfo,$adapt)= @_;
	my ($info,@line,@item);
 
	open (IN,"$infile");

	while (<IN>) {
 		chomp;
 		next if (/^(\#|\s*$)/);
# 		if (/^\@SQ\tSN:(\S+) LN:(\d+)/) {
 		if (/^\@SQ/) {
			$cont= $1 if (/SN\:(\S+)/);
			$clen= $1 if (/LN\:(\d+)/);
		#	($cont,$clen)= ($1,$2) if (/SN\:(.+) LN\:(.+)/);
			$terminfo->{$cont}->{conlen}= $clen;
            #         print "$cont $clen\n"; exit;
            
		} elsif (/^[A-Z]/) {
            
 			@item= split/\t/;
            
            next if ($item[2] eq '*');

            $cont=$item[2];
			$start=$item[3];
			$len= length($item[9]);
## estimate end point
			$cigar=$item[5];
            
            @sclips=();
            if ($cigar=~ /(\d+)[SH]/) {
                @sclips= $cigar=~ /(\d+)S/g;
                next;
            }
            
			@insert= $cigar=~ /(\d+)I/g;
			$minussum=0;
			$minussum+= $_ foreach (@sclips,@insert);
			$end= $start+$len-1-$minussum;
            
            if ($inform eq "both") {
                $terminfo->{$cont}->{posit}->{$start}++;
                $terminfo->{$cont}->{negat}->{$end}++;
            } elsif ($inform eq "fwd") {
                
                if ($adapt eq "upm") {
                    $terminfo->{$cont}->{posit}->{$start}++;
                } elsif ($adapt eq "u2") {
                     $terminfo->{$cont}->{posit}->{$end}++;
                }
            } elsif ($inform eq "rev") {
                
                if ($adapt eq "upm") {
                    $terminfo->{$cont}->{negat}->{$end}++;
                } elsif ($adapt eq "u2") {
                    $terminfo->{$cont}->{negat}->{$start}++;
                }
            }
		}
	}

	return ($terminfo);
}



sub output_count_term {
	my ($terminfo,$adapt)= @_;
    my (@titles);
    
	open (OUT,">out_countTerm_summary.txt");
    
    @titles= qw/adapt cont conlen/;
    push @titles, qw/1st_st_pos 1st_st_cnt 1st_st_pval 2nd_st_pos 2nd_st_cnt 2nd_st_pval 3rd_st_pos 3rd_st_cnt 3rd_st_pval post_num/;
    push @titles, qw/1st_ed_pos 1st_ed_cnt 1st_ed_pval 2nd_ed_pos 2nd_ed_cnt 2nd_ed_pval 3rd_ed_pos 3rd_ed_cnt 3rd_ed_pval post_num/;
	print OUT (join "\t",@titles)."\n";
    
	open (LOG,">out_countTerm_list.txt");
    print LOG (join "\t",qw(adapt cont clen loc count type))."\n";
    
    @conts= (sort keys %{$terminfo});
    $contnum= @conts;
    
	$cno= 0;
	foreach $cont ( sort keys %{$terminfo} ) {
		$cno++;
		$conlen= $terminfo->{$cont}->{conlen};
        @sumone= ($adapt,$cont,$conlen);
        
        foreach $post (sort{$a<=>$b} keys %{$terminfo->{$cont}->{posit}} ) {
            $countone= $terminfo->{$cont}->{posit}->{$post};
            print LOG "$adapt\t$cont\t$conlen\t$post\t$countone\tposit\n";
        }
 
        foreach $post (sort{$a<=>$b} keys %{$terminfo->{$cont}->{negat}} ) {
            $countone= $terminfo->{$cont}->{negat}->{$post};
            print LOG "$adapt\t$cont\t$conlen\t$post\t$countone\tnegat\n";
        }


## NEGATIVE ( UPM: 3'-TERM  U2: 5'-TERM)

        $min_postnum=3;
        
## 5'-TERM (UPM: positive U2: negative)
        my($outlier_st,$stats_st,$cc_st);
        
        $strandx= $adapt eq "upm" ? "posit" : "negat";
        @postx= sort{$a<=>$b} keys %{$terminfo->{$cont}->{$strandx}};
        $postx_num=@postx;
        
        ## SELCT >= 3 COUNT
#        @targetpostx=();
#        @targetpostx= grep{$terminfo->{$cont}->{$strandx}->{$_} > 2 }@postx;

        push @sumone,"nd","nd","nd","nd","nd","nd","nd","nd","nd",$postx_num;
        if ($postx_num >= $min_postnum) {
            ($outlier_st,$stats_st)=&calculate_Smirnov_Grubbs_test($terminfo->{$cont},[@postx],"start",$strandx);
            print LOG2 (join "\t",($adapt,$strandx,$cont,$conlen,@$_))."\n" foreach (@$outlier_st);

            @sumone[3..11]= (@{$$outlier_st[0]},@{$$outlier_st[1]},@{$$outlier_st[2]});

        } else {
            if ($postx_num > 0) {
                @sumone[3..5]= ($postx[0],$terminfo->{$cont}->{$strandx}->{$postx[0]},"na");
                if ($postx_num > 1) {
                    @sumone[6..8]= ($postx[1],$terminfo->{$cont}->{$strandx}->{$postx[1]},"na");
                }
            }
                
        }
## 3'-TERM (UPM: negative U2: positive)
        my($outlier_ed,$stats_ed,$cc_ed);
        $strandy= $adapt eq "upm" ? "negat" : "posit";
        @posty= sort{$b<=>$a} keys %{$terminfo->{$cont}->{$strandy}};
        $posty_num=@posty;
        
        ## SELCT >= 3 COUNT
#        @targetposty=();
#        @targetposty= grep{$terminfo->{$cont}->{$strandy}->{$_} > 2 }@posty;

        push @sumone,"nd","nd","nd","nd","nd","nd","nd","nd","nd",$posty_num;
        if ($posty_num >= $min_postnum) {
 
            ($outlier_ed,$stats_ed)=&calculate_Smirnov_Grubbs_test($terminfo->{$cont},[@posty],"end",$strandy);
            print LOG2 (join "\t",($adapt,$strandy,$cont,$conlen,@$_))."\n" foreach (@$outlier_ed);
            
            @sumone[13..21]= (@{$$outlier_ed[0]},@{$$outlier_ed[1]},@{$$outlier_ed[2]});

        } else {
            if ($posty_num > 0) {
                @sumone[13..15]= ($posty[0],$terminfo->{$cont}->{$strandy}->{$posty[0]},"na");
                if ($posty_num > 1) {
                    @sumone[16..18]= ($posty[1],$terminfo->{$cont}->{$strandy}->{$posty[1]},"na");
                }
            }
        }
 
#        print "#### $adapt $postx_num $posty_num\n";
        
        
		print OUT (join "\t",(@sumone))."\n";
        
        print "$cno / $contnum $cont\n";
        
	}

	return ($cno);

}

sub calculate_Smirnov_Grubbs_test {
    my ($postinfo,$post_arr,$term,$strand)= @_;
    my (@posts,$sum,$post_num,$max_post,$max_count);
    my ($postone,$countone,$ave,$cov,$std);
    my (@stats,@outlier,@counts,$SG_values,$bb,$cc,$summary_SG,$pval);

    if ($term eq "start") {
            @posts= sort{$a<=>$b}@$post_arr;
    } elsif ($term eq "end") {
            @posts= sort{$b<=>$a}@$post_arr;
    }

### CAL STATICS
    ($sum,$post_num)= (0,0);
    ($max_post,$max_count)=(0,0);
    
    foreach $postone (@$post_arr) {
        $post_num++;
        $countone= $postinfo->{$strand}->{$postone};

        $sum+= $countone;
        if ($countone > $max_count) {
            ($max_post,$max_count)=($postone,$countone);
        }
    }
    $ave= sprintf("%.1f",$sum/$post_num);
    $cov= 0;
    foreach $postone (@$post_arr) {
        $countone= $postinfo->{$strand}->{$postone};
        $cov+= ($countone-$ave)**2;
    }
    $std= sprintf("%.2f",sqrt($cov/$post_num));
    
    @stats=($post_num,$ave,$max_post,$max_count);

### CAL OUTLIER
    
    @outlier= ();
    @counts= values %{$postinfo->{$strand}};
    $SG_values= "c(".(join ",",@counts).")";

    
      ###  FIRST POSIT & SIGNIFICANT POSIT
    
    $bb=0;
    $cc=0;
    @effect_post=();
    foreach $postone (@posts) {
        $countone= $postinfo->{$strand}->{$postone};

        next if ($countone < 3);
        push @effect_post,$postone;
        
        $bb++;

        # FIRST POST
        if ($bb <= 2) {
            ($summary_SG)= &Smirnov_Grubbs_test($SG_values,$countone);
            $pval= $$summary_SG[5];
            push @outlier,[$postone,$countone,$pval];
            $cc++;

        } else {
            ($summary_SG)= &Smirnov_Grubbs_test($SG_values,$countone);
            $pval= $$summary_SG[5];
            if ($pval ne "NA") {
            ## OUTLIER POST ( P-val < 0.05)
                if ($pval < 5e-2) {
                    push @outlier,[$postone,$countone,$pval];
                    $cc++;
                    last if ($cc >= 3);
                }
            }
        }
    }
    
#    print "$strand $_ $postinfo->{$strand}->{$_}\n" foreach (@posts);

    if ($cc == 2) {
        if (@effect_post > 2) {
            push @outlier,[$effect_post[2],$postinfo->{$strand}->{$effect_post[2]},"nd"];
        } else {
            push @outlier,["nd","nd","nd"];
        }
    } elsif ($cc == 1) {
        push @outlier,["nd","nd","nd"],["nd","nd","nd"];
    } elsif ($cc == 0) {
        push @outlier,[$posts[0],$postinfo->{$strand}->{$posts[0]},"na"],
        ["nd","nd","nd"],["nd","nd","nd"];
    }
    
    
    return ([@outlier],[@stats]);

}



sub Smirnov_Grubbs_test {
	my ($values,$query)= @_;

	$script = "temp_SG.R";
	open(SCRIPT, ">$script");

#####
	print SCRIPT <<EOF;

	SG <- function(x, y) {
	method <- "Smirnov‐Grubbs test"
	data.name <- paste(c("min(", "max("), deparse(substitute(x)), ") = ", range(x, na.rm=TRUE), sep="")
	x <- x[!is.na(x)]
	n <- length(x)
    xrange<-c(min(x),max(x),y)
	t <- abs(xrange-mean(x))/sd(x)
	p <- n*pt(sqrt((n-2)/((n-1)^2/t^2/n-1)), n-2, lower.tail=FALSE)
	p <- sapply(p, function(p0) min(p0, 1))

	result <- list(method=method, parameter=c(df=n-2))
	result1 <- structure(c(result,  list(data.name=data.name[1], 
          statistic=c(t=t[1]), p.value=p[1])), class="htest")
	result2 <- structure(c(result,  list(data.name=data.name[2], 
          statistic=c(t=t[2]), p.value=p[2])), class="htest")
    result3 <- structure(c(result,  list(data.name="query",
          statistic=c(t=t[3]), p.value=p[3])), class="htest")
	return(structure(list(result1, result2, result3), class="SG"))

	}

	array<-$values
	query<-$query

	SG(array,query)

EOF
	close(SCRIPT);

#####
# $reault= system("R --vanilla --slave < $script");
#($Rstat = $?/256) && die "Aborted in R with status $Rstat.?n";

	@reaults= `R --vanilla --slave < $script`;
	unlink $script;
	
	@resum=();
	foreach (@reaults) {
		chomp;
		if (/min\(array\) \= (\d+)/) {
			$min=$1;
			push @resum,$min;
		} elsif (/p-value .+ (\S+)/) {
			$pval=$1;
			push @resum,$pval;
		} elsif (/max\(array\) \= (\d+)/) {
			$max=$1;
			push @resum,$max;
		} elsif (/query/) {
			push @resum,"query";
		}

	}

	return ([@resum]);

}
