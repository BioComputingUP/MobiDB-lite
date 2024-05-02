#!/usr/bin/perl

use warnings;
use strict;


{
	my %args;
	my $seq=check_args(\%args,@ARGV);
	my $s=get_seq($seq);

	my $path=$0;
	$path=~s/anchor\.pl$//;
	if($path=~/^$/) {$path='.'}
	else {$path=~s/\/$//}
	
	my @anchor;
        my $comm;	
	
	if($args{"DIR"}==0) {
		if(!-e "$path/anchor") {die "Error: $path/anchor not found.\n"}
		if(!-e "$path/anchordata") {die "Error: $path/anchordata not found.\n"}
		$comm="$path/anchor $seq -d $path";
	}
	if($args{"DIR"}==1) {
		if(!-e $args{"DIR_NAME"}."/anchor") {die "Error: ",$args{"DIR_NAME"},"/anchor not found.\n"}
		if(!-e $args{"DIR_NAME"}."/anchordata") {die "Error: ",$args{"DIR_NAME"},"/anchordata not found.\n"}
		$comm='$args{"DIR_NAME"}/anchor $seq -d $args{\"DIR_NAME\"}';
	}

        if ($args{"VERBOSE"}==1) {
            $comm=$comm." -v";
	}

	@anchor=`$comm`;

	if($args{"MOTIF"}==1) {
		
	    my $ln=0;
	    while() {
		print $anchor[$ln];
		$ln++;
		if($anchor[$ln]=~/^# Prediction profile output:/) {last;}
	    }
	    
	    
	    my @motifs;
	    my $motif_count=read_motifs(\@motifs,$args{"MOTIF_FILE"});
	    
	    
	    print "# Motifs\n"; 
	    printf "#\t%-15s\t%8s\t%8s\t%s\n","Name","Start","End","Instance";
	    for(my $i=0;$i<$motif_count;$i++) {
		my (@start,@end);
		check_motif($motifs[$i]{"pattern"});
		matching(\@start,\@end,$s,$motifs[$i]{"pattern"});
		for(my $j=0;$j<@start;$j++) {
		    printf ("#\t%-15s\t%8d\t%8d\t",$motifs[$i]{"name"},$start[$j],$end[$j]);
		    if($start[$j]!=1) {print "..."}
		    print substr($s,$start[$j]-1,$end[$j]-$start[$j]+1);
		    if($end[$j]!=length($seq)) {print "..."}
		    print "\n";
		}
	    }
	    
	    print "#\n";
	    for(my $i=$ln;$i<@anchor;$i++) {print $anchor[$i]}
	    
	}
	else {print @anchor}
	
}

sub get_seq {
    
    my $seq=$_[0];
    my $s;
    
    open(SEQ,$seq) || die "Error: file $seq cannot be opened.\n";
    foreach my $line (<SEQ>) {
	if($line=~/^>/) {next;}
	$s.=$line;
    }
    close(SEQ);
    
    $s=~s/\s//g;
    
    return($s);
}

sub check_motif {

    my $mot=$_[0];
    
    if(($mot=~/\(\s*\?{1,2}\s*\{/)||($mot=~/\$\s*\{/)||($mot=~/@\s*\{/)) {
	die "Error: possible malicious code in $mot\n";
    }
    eval { qr/$mot/};
    if($@) {die "Error: $mot is not a valid regular expression.\n"}
}


sub check_args {
	
    my ($args,@argv)=@_;
    my $seq;
    
    $args->{"MOTIF"}=0;
    $args->{"VERBOSE"}=0;
    $args->{"DIR"}=0;
    $args->{"FILTER"}=0;
    
    if($#argv<0) {
	print STDERR "Usage: $0 (options -m/-d/-v) (sequence)\n";
	exit;
    }
	
    if($#argv==0) {
	if($argv[0]=~/-[a-z]/) {
	    print STDERR "Error: no sequence specified.\n";
	    exit;
	}
	return($argv[0]);
    }
    
    for(my $i=0;$i<@argv;$i++) {
	if(!$argv[0]) {
	    print STDERR "Error: no sequence specified.\n";
	    exit;
	}
	if($argv[$i] eq '-v') {
	    $args->{"VERBOSE"}=1;
	    splice(@argv,$i,1);
	    $i--;
	    next;
	}
	if($argv[$i] eq '-m') {
	    $args->{"MOTIF"}=1;
	    splice(@argv,$i,1);
	    if((!$argv[$i])||($argv[$i]=~/-[a-z]/)) {
		print STDERR "Error: no motif file specified after -m\n";
		exit;
	    }
	    else {$args->{"MOTIF_FILE"}=$argv[$i]}
	    splice(@argv,$i,1);
	    $i--;
	    next;
	}
	if($argv[$i] eq '-d') {
	    $args->{"DIR"}=1;
	    splice(@argv,$i,1);
	    if((!$argv[$i])||($argv[$i]=~/a-z/)) {
		print STDERR "Error: no directory specified after -d\n";
		exit;
	    }
	    else {$args->{"DIR_NAME"}=$argv[$i]}
	    splice(@argv,$i,1);
	    $i--;
	    next;
	}
    }
    
    if(($args->{"MOTIF"})&&(! -e $args->{"MOTIF_FILE"})) {
	print STDERR "Error: cannot open motif file ",$args->{"MOTIF_FILE"},"\n";
	exit;
    }
    
    if(!$argv[-1]) {
	print STDERR "Error: no sequence specified.\n";
	exit;
    }
    return($argv[-1]);
}


sub matching {
    
    my ($start,$end,$target,$pattern)=@_;
    my $found=0;
    my $chopped=0;
    
    if($target!~/$pattern/i) { return ;}
    
    while ($target) {
	if($chopped==0) {
	    $target=~/$pattern/i;
	    $start->[$found]=$-[0]+1;
	    $end->[$found]=$+[0];
	    $found++;
	}
	else {
	    if ($target=~/$pattern/i) {
		if (($-[0]+1+$chopped != $start->[$found-1])&&($+[0]+$chopped != $end->[$found-1])) {
		    $start->[$found]=$-[0]+1+$chopped;
		    $end->[$found]=$+[0]+$chopped;
		    $found++;
		}			
	    }
	}
	$target=substr($target,1);
	$chopped++;
    }
}




sub read_motifs {
    my ($motifs,$file)=@_;
    
    open(MOT,$file) || die "# Cannot open \"$file\" for reading - Motif search aborted.\n";
    
    my $count=0;
    my $lc=0;
    foreach my $line (<MOT>) {
	$lc++;
	if($line=~/^#/) { next }
	if($line=~/^$/) { next }
	chomp $line;
	$line=~s/^\s*//;
	my @tmp=split('\s+',$line);
	
	$motifs->[$count]{"pattern"}=$tmp[0];
	if($tmp[1]) {
	    $motifs->[$count]{"name"}=substr($line,length($motifs->[$count]{"pattern"}));
	    $motifs->[$count]{"name"}=~s/^\s+//;
	}
	else {$motifs->[$count]{"name"}=$tmp[0]}
	
	$count++;
    }
    close(MOT);
    
    return($count);
}

