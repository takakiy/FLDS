#!/usr/bin/perl -w

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

$terminfo= {};
if ( exists $option->{-b}) {

    $inform= "both";
    ($terminfo)= &create_terminfo($option->{-b}->[0],$inform,$terminfo);
} elsif (exists $option->{-f} || exists $option->{-r}) {
    $inform= "fwd";
    ($terminfo)= &create_terminfo($option->{-f}->[0],$inform,$terminfo);
    $inform= "rev";
    ($terminfo)= &create_terminfo($option->{-r}->[0],$inform,$terminfo);

}

print "### FINISHED READ SAM\n";

open (LOG2,">out_countTerm_outlier.txt");


($cno)= &output_count_term($terminfo);

print "$cno\n";
print "### FINISHED CAL OUTLIER\n";

##=== === === === === === === === === === === === === === === === ===
sub create_terminfo {
	my ($infile,$inform,$terminfo)= @_;
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
                $terminfo->{$cont}->{start}->{$start}++;
                $terminfo->{$cont}->{end}->{$end}++;
            } elsif ($inform eq "fwd") {
                $terminfo->{$cont}->{start}->{$start}++;
            } elsif ($inform eq "rev") {
                $terminfo->{$cont}->{end}->{$end}++;
            }
		}
	}

	return ($terminfo);
}



sub output_count_term {
	my ($terminfo)= @_;
    my (@titles);
    
	open (OUT,">out_countTerm_summary.txt");
    
    @titles= qw/cont conlen/;
    push @titles, qw/1st_st_pos 1st_st_cnt 1st_st_pval 2nd_st_pos 2nd_st_cnt 2nd_st_pval/;
    push @titles, qw/t_start av_start max_st_pos max_st_cnt/;
    push @titles, qw/1st_ed_pos 1st_ed_cnt 1st_ed_pval 2nd_ed_pos 2nd_ed_cnt 2nd_ed_pval/;
    push @titles, qw/t_end av_end max_end_pos max_end_cnt/;
	print OUT (join "\t",@titles)."\n";
    
	open (LOG,">out_countTerm_list.txt");
    print LOG (join "\t",qw(cont clen loc count type))."\n";
    
    @conts= (sort keys %{$terminfo});
    $contnum= @conts;

	$cno= 0;
	foreach $cont ( sort keys %{$terminfo} ) {
		$cno++;
		$conlen= $terminfo->{$cont}->{conlen};
        @sumone= ($cont,$conlen);
        
## START
        foreach $post (sort{$a<=>$b} keys %{$terminfo->{$cont}->{start}} ) {
            $countone= $terminfo->{$cont}->{start}->{$post};
            print LOG "$cont\t$conlen\t$post\t$countone\tstart\n";
        }

        if (keys %{$terminfo->{$cont}->{start}} > 4) {
            
            ($outlier_st,$stats)=&calculate_Smirnov_Grubbs_test($terminfo->{$cont},"start");
            print LOG2 (join "\t",("start",$cont,$conlen,@$_))."\n" foreach (@$outlier_st);

            if (@$outlier_st > 1) {
                push @sumone,@{$$outlier_st[0]},@{$$outlier_st[1]};
            } else {
                push @sumone,@{$$outlier_st[0]},"nd","nd","nd";
            }
            push @sumone,@$stats;
            
		} else {
			push @sumone,"nd","nd","nd","nd","nd","nd","nd","nd","nd","nd";
		}
        
## END
        foreach $post (sort{$a<=>$b} keys %{$terminfo->{$cont}->{end}} ) {
            $countone= $terminfo->{$cont}->{end}->{$post};
            print LOG "$cont\t$conlen\t$post\t$countone\tend\n";
        }

		if (keys %{$terminfo->{$cont}->{end}} > 4) {

            ($outlier_ed,$stats)=&calculate_Smirnov_Grubbs_test($terminfo->{$cont},"end");
            print LOG2 (join "\t",("end",$cont,$conlen,@$_))."\n" foreach (@$outlier_ed);
            
            if (@$outlier_ed > 1) {
                push @sumone,@{$$outlier_ed[0]},@{$$outlier_ed[1]};
            } else {
                push @sumone,@{$$outlier_ed[0]},"nd","nd","nd";
            }

            push @sumone,@$stats;

		} else {
            push @sumone,"nd","nd","nd","nd","nd","nd","nd","nd","nd","nd";
		}

		print OUT (join "\t",@sumone)."\n";
        
        print "$cno / $contnum $cont\n";
        
	}

	return ($cno);

}

sub calculate_Smirnov_Grubbs_test {
    my ($postinfo,$order)= @_;
    my (@counts,@posts,@outlier,$postone,$SG_values,$summary_SG,$pval);


# print "## $order\n";
    
    if ($order eq "start") {
        @posts= sort{$a<=>$b} keys %{$postinfo->{$order}};
    } elsif ($order eq "end") {
        @posts= sort{$b<=>$a} keys %{$postinfo->{$order}};
    }

### CAL STATICS
    ($sum,$post_num)= (0,0);
    ($max_post,$max_count)=(0,0);
    
    foreach $postone (@posts) {
        $post_num++;
        $countone= $postinfo->{$order}->{$postone};

        $sum+= $countone;
        if ($countone > $max_count) {
            ($max_post,$max_count)=($postone,$countone);
        }
    }
    $ave= sprintf("%.1f",$sum/$post_num);
    $cov= 0;
    foreach $postone (@posts) {
        $countone= $postinfo->{$order}->{$postone};
        $cov+= ($countone-$ave)**2;
    }
    $std= sprintf("%.2f",sqrt($cov/$post_num));
    
    @stats=($post_num,$ave,$max_post,$max_count);

### CAL OUTLIER
    
    @outlier= ();
    @counts= values %{$postinfo->{$order}};
    $SG_values= "c(".(join ",",@counts).")";

    $bb=0;
    $cc=0;
    foreach $postone (@posts) {
        $bb++;
        $countone= $postinfo->{$order}->{$postone};
 #print "$postone $countone $ave $std\n";
        # FIRST POST
        if ($bb == 1) {
            ($summary_SG)= &Smirnov_Grubbs_test($SG_values,$countone);
            $pval= $$summary_SG[5];
            
            push @outlier,[$postone,$countone,$pval];
            if ($pval ne "NA") {
                $cc++ if ($pval < 5e-2);
            }
            next;
        } else {
            
            ### SKIP IN (COUNT-AVERAGE) > STD
            
            $pval= 1;
            if (abs($countone - $ave) > $std ) {
                ($summary_SG)= &Smirnov_Grubbs_test($SG_values,$countone);
                $pval= $$summary_SG[5];

            }
            
            ## OUTLIER POST ( P-val < 0.05)
            if ($pval < 5e-2) {
                push @outlier,[$postone,$countone,$pval];
                $cc++;
                last if ($cc > 2);
            
            }
        }

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
