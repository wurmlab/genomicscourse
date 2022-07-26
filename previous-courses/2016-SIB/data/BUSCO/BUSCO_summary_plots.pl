# Author: Robert Waterhouse
#!/usr/local/bin/perl

use strict;

if(scalar(@ARGV)!=1) {
    print LOG "Error: must specify name of full_table results file (or folder of files)\n";
    print LOG "E.g. BUSCO_summary_plots.pl my_full_table_files\n";
    print LOG "E.g. BUSCO_summary_plots.pl full_table_assembly1\n";
    print LOG "E.g. BUSCO_summary_plots.pl full_table_geneset1\n";
    exit();
}

open(LOG,">BUSCO_PLOTS.log.txt") || die $!;
open(OUT,">BUSCO_R_PLOTS.R.txt") || die $!;
my @dsets=();

if(-d $ARGV[0]) {
    my $directory=$ARGV[0];
    if($directory!~/\/$/) { $directory=$directory.'/'; }
    print LOG "Reading contents of $directory ...\n";
    my @allfiles=split(/\s+/,`ls $directory`);
    foreach my $file (@allfiles) {
	my $fullfile=$directory.$file;
	my ($TYP,$B2S,$B2H)=parse_full_table_file($fullfile);
	my $dset=print_R_script($file,$TYP,$B2S,$B2H);
	push(@dsets,$dset);
    }
}
elsif(-f $ARGV[0]) {
    my ($TYP,$B2S,$B2H)=parse_full_table_file($ARGV[0]);
    my $dset=print_R_script($ARGV[0],$TYP,$B2S,$B2H);
    push(@dsets,$dset);
}

# set your colours here
print OUT "buscocols<-c('#0099ff','#3333ff','#ffcc00','#ff3300')\n";    # Traditional
#print OUT "buscocols<-c('#00ff00','#99cc00','#ffff00','#ff0066')\n";    # Bright

print OUT "pdf(\"BUSCO_R_PLOTS.pdf\")\n";
print OUT "allsets<-matrix(c(" . join(',',sort(@dsets)) . "), ncol=" . scalar(@dsets) . ", nrow=4)\n";
print OUT "allsetsP<-(t(t(allsets)/colSums(allsets)))*100\n";
print OUT "barplot(allsetsP, names.arg=c('" . join("','",sort(@dsets)) . "'), las=2, cex.names=0.6, col=buscocols, space=0.3, main='BUSCO Assessment Results', ylab='% BUSCOs', ylim=c(0,115), cex.axis=0.75)\n";

print OUT "xpos<-barplot(allsetsP, space=0.3, plot=FALSE)\n";
foreach my $ds (sort(@dsets)) { print OUT "$ds\_t<-sum($ds)\n"; }
print OUT "clabs<-c(" . join("_t,",sort(@dsets)) . "_t)\n";
print OUT "text(x=xpos,y=101, labels=clabs, cex=0.5, srt=50, pos=3)\n";
print OUT "legend('top', c('Complete (1)','Complete (>1)','Fragmented','Missing'), fill=buscocols, xjust=0.5, bty='n', horiz=TRUE, cex=0.75)\n";

print OUT "dev.off()\n";
close(LOG);
close(OUT);
my $rcmd="R CMD BATCH BUSCO_R_PLOTS.R.txt";
system($rcmd);


sub parse_full_table_file {

    open(IN,$_[0]) || die "Cannot open $_[0] for reading: $!";
    my @lines=<IN>;
    close(IN);
    
    my $type='unrecognised';
    my $header=shift(@lines);
    chomp($header);
    if($header=~/^#BUSCO_group\s+Status\s+Scaffold\s+Start\s+End\s+Bitscore\s+Length$/) { $type='genome'; }
    elsif($header=~/^#BUSCO_group\s+Status\s+Gene\s+Bitscore\s+Length$/) { $type='geneset'; }
    if($type eq 'unrecognised') {
	print LOG "File type not recognised from the header in $_[0], please check file formats and try again.\n";
	exit();
    }

    undef my %busco2status;
    undef my %busco2hitloc;
    undef my %buscohit2scr;
    my @badlines=();
    my $linecount=1;

    foreach my $line (@lines) {
	chomp($line);
	$linecount++;
	my @datas=split(/\t/,$line);
	if($type eq 'genome' && $datas[1] ne 'Missing' && scalar(@datas)!=7) {
	    print LOG "Looks like genome results, but this line does not comply:\n$line\n";
	    print LOG "Line $linecount from $_[0]\n";
	    push(@badlines,$line);
	}
	elsif($type eq 'geneset' && $datas[1] ne 'Missing' && scalar(@datas)!=5) {
	    print LOG "Looks like geneset results, but this line does not comply:\n$line\n";
	    print LOG "Line $linecount from $_[0]\n";
	    push(@badlines,$line);
	}
	if($datas[1] eq 'Missing' && scalar(@datas)!=2) {
	    print LOG "Looks like $type results, but this line does not comply:\n$line\n";
	    print LOG "Line $linecount from $_[0]\n";
	    push(@badlines,$line);
	}
	if(scalar(@badlines)>5) {
	    print LOG "Exiting, too many non-compliant lines in $_[0], please check file formats and try again.\n";
	    exit();
	}
	$busco2status{$datas[0]}=$datas[1];
	my $score=6;
	if($type eq 'geneset') { $score=4; }
	if(defined($busco2hitloc{$datas[0]})) {
	    if($datas[$score]>$buscohit2scr{$datas[0]}) { $busco2hitloc{$datas[0]}=$datas[2]; $buscohit2scr{$datas[0]}=$datas[$score]; }
	}
	else { 
	    $busco2hitloc{$datas[0]}=$datas[2]; 
	    $buscohit2scr{$datas[0]}=$datas[$score];
	}
    }

    print LOG $_[0] . "\twith " . scalar(@lines) . " lines\tof type $type asessement\n";

    return ($type, \%busco2status, \%busco2hitloc);
}



sub print_R_script {
    my ($file,$type,$stat,$loca) = @_;

    my $name='file';
    if($file=~/full_table_(\S+)/) { $name=$1; }
    $name=~s/\-//g;
    $name=~s/\_//g;
    $name=~s/\.//g;
    $name=~s/\+//g;
    $name=~s/\(//g;
    $name=~s/\)//g;
    $name=~s/\[//g;
    $name=~s/\]//g;
    $name=~s/\{//g;
    $name=~s/\}//g;
    $name=~s/\///g;
    $name=~s/\\//g;
    $name=~s/\*//g;
    $name=~s/\&//g;
    $name=~s/\%//g;
    $name=~s/\$//g;
    $name=~s/\#//g;
    $name=~s/\@//g;
    $name=~s/\!//g;
    $name=~s/\?//g;
    $name=~s/\|//g;

    $name=uc($name);
    
    undef my %status2count;
    $status2count{'Complete'}=0;
    $status2count{'Duplicated'}=0;
    $status2count{'Fragmented'}=0;
    $status2count{'Missing'}=0;
    
    foreach my $busco (sort keys %{$stat}) { $status2count{$stat->{$busco}}++; }
    
    print LOG "\tComplete SC: " . $status2count{'Complete'} . "\t";
    print LOG "Complete MC: " . $status2count{'Duplicated'} . "\t";
    print LOG "Fragmented: " . $status2count{'Fragmented'} . "\t";
    print LOG "Missing: " . $status2count{'Missing'} . "\n";

    print OUT "$name<-c(" . $status2count{'Complete'} . "," . $status2count{'Duplicated'} . "," . $status2count{'Fragmented'} . "," . $status2count{'Missing'} . ")\n";

    return($name);
}
