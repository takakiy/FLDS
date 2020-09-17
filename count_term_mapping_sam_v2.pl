#!/usr/bin/perl -w

if ((@ARGV == 0) || ($ARGV[0]=~/-h/)) {
	print "OPTION -i infile\n";
#	exit;

}

foreach (@ARGV){if (/^\-/){$key= $_} else {push @{$option->{$key}},$_}};

($terminfo)= &create_terminfo($option->{-i}->[0]);

 print "### FINISHED READ SAM\n";

($cno)= &output_count_term($terminfo);

print "$cno\n";
print "### FINISHED CAL OUTLIER\n";

sub create_terminfo {
	my ($infile)= @_;
	my ($info,@line,@item);
 
	$terminfo= {};
	open (IN,"$infile");

	open (LOG2,">out_count_term.logs2.txt");

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
			@sclips= $cigar=~ /(\d+)S/g;
			@insert= $cigar=~ /(\d+)I/g;
			$minussum=0;
			$minussum+= $_ foreach (@sclips,@insert);
			$end= $start+$len-1-$minussum;

			$terminfo->{$cont}->{start}->{$start}++;
			$terminfo->{$cont}->{end}->{$end}++;

		}
	}

	return ($terminfo);
}



sub output_count_term {
	my ($terminfo)= @_;
    my (@titles);
    
	open (OUT,">out_count_term.txt");
    
    @titles= qw/cont conlen/;
    push @titles, qw/1st_st_pos 1st_st_cnt 1st_st_pval 2nd_st_pos 2nd_st_cnt 2nd_st_pval/;
    push @titles, qw/t_start av_start max_st_pos max_st_cnt/;
    push @titles, qw/1st_ed_pos 1st_ed_cnt 1st_ed_pval 2nd_ed_pos 2nd_ed_cnt 2nd_ed_pval/;
    push @titles, qw/t_end av_end max_end_pos max_end_cnt/;
	print OUT (join "\t",@titles)."\n";
    
	open (LOG,">out_count_term.logs.txt");
	$cno= 0;
	foreach $cont ( sort keys %{$terminfo} ) {
		$cno++;
		$conlen= $terminfo->{$cont}->{conlen};
        @sumone= ($cont,$conlen);

## START
		if (keys %{$terminfo->{$cont}->{start}} > 0) {
            
            ($outlier_st)=&calculate_Smirnov_Grubbs_test($terminfo->{$cont},"start");
            
            foreach (@$outlier_st) {
                print LOG2 (join "\t",("start",$cont,$conlen,@$_))."\n";
                
            }
            
            if (@$outlier_st > 1) {
                push @sumone,@{$$outlier_st[0]},@{$$outlier_st[1]};
            } else {
                push @sumone,@{$$outlier_st[0]},"nd","nd","nd";
            }
            
            ($sum_st,$count_st)= (0,0);
            ($max_post_st,$max_count_st)=(0,0);
			foreach $post (sort{$a<=>$b} keys %{$terminfo->{$cont}->{start}} ) {
				$count_st++;
				$sum_st+= $terminfo->{$cont}->{start}->{$post};
				print LOG "$cont\t$conlen\t$post\t$terminfo->{$cont}->{start}->{$post}\tstart\n";
                
                if ($terminfo->{$cont}->{start}->{$post} > $max_count_st) {
                    ($max_post_st,$max_count_st)=($post,$terminfo->{$cont}->{start}->{$post});
                }
			}
			$ave_st= sprintf("%.1f",$sum_st/$count_st);
            
            push @sumone,$count_st,$ave_st,$max_post_st,$max_count_st;

		} else {
			push @sumone,"nd","nd","nd","nd","nd","nd","nd","nd","nd","nd";
		}
## END

		if (keys %{$terminfo->{$cont}->{end}} > 0) {

            ($outlier_ed)=&calculate_Smirnov_Grubbs_test($terminfo->{$cont},"end");
            
            foreach (@$outlier_ed) {
                print LOG2 (join "\t",("end",$cont,$conlen,@$_))."\n";
            }
            
            if (@$outlier_ed > 1) {
                push @sumone,@{$$outlier_ed[0]},@{$$outlier_ed[1]};
            } else {
                push @sumone,@{$$outlier_ed[0]},"nd","nd","nd";
            }

            ($sum_ed,$count_ed)= (0,0);
            ($max_post_ed,$max_count_ed)=(0,0);
            foreach $post (sort{$a<=>$b} keys %{$terminfo->{$cont}->{end}} ) {
				$count_ed++;
				$sum_ed+= $terminfo->{$cont}->{end}->{$post};
				print LOG "$cont\t$conlen\t$post\t$terminfo->{$cont}->{end}->{$post}\tend\n";
                if ($terminfo->{$cont}->{end}->{$post} > $max_count_ed) {
                    ($max_post_ed,$max_count_ed)=($post,$terminfo->{$cont}->{end}->{$post});
                }
			}
			$ave_ed= sprintf("%.1f",$sum_ed/$count_ed);
            push @sumone,$count_ed,$ave_ed,$max_post_ed,$max_count_ed;

		} else {
            push @sumone,"nd","nd","nd","nd","nd","nd","nd","nd","nd","nd";
		}

		print OUT (join "\t",@sumone)."\n";

	}

	return ($cno);

}

sub calculate_Smirnov_Grubbs_test {
    my ($postinfo,$order)= @_;
    my (@counts,@posts,@outlier,$postone,$SG_values,$summary_SG,$pval);

    @counts= values %{$postinfo->{$order}};
    if ($order eq "start") {
        @posts= sort{$a<=>$b} keys %{$postinfo->{$order}};
    } elsif ($order eq "end") {
        @posts= sort{$b<=>$a} keys %{$postinfo->{$order}};
    }
    @outlier= ();
    
    $SG_values= "c(".(join ",",@counts).")";

    $bb=0;
    $cc=0;
    foreach $postone (@posts) {
        $bb++;
        $postvalue= $postinfo->{$order}->{$postone};
        ($summary_SG)= &Smirnov_Grubbs_test($SG_values,$postvalue);
        $pval= $$summary_SG[5];
        
        # FIRST POST
        if ($bb == 1) {
            push @outlier,[$postone,$postvalue,$pval];
            $cc++ if ($pval < 1e-2);
            next;
        }
        
        ## OUTLIER POST
        if ($pval < 1e-2) {
            push @outlier,[$postone,$postvalue,$pval];
            $cc++;
            last if ($cc > 2);
            
        }

    }
    
    return ([@outlier]);

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
