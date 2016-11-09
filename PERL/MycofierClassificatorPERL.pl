#!/usr/bin/perl -w
# Mycofier: A new Machine Learning-based classifier for fungal ITS sequences
# Delgado, et al. 2016
# PUBLICACION ISSUE
# Corporacion CorpoGen & Universidad de Los Andes

####### Location Of WEKA jar FILE ########
my $model = "dbnosing_sps.model.fasta";
##########################################
use strict;
use Getopt::Long;
use Algorithm::NaiveBayes;
use threads;
use Thread::Queue;
use Fcntl qw/:flock/; 

my $q = Thread::Queue->new(); 
my $nb = Algorithm::NaiveBayes->new;
my $infile;
my $outfile;
my $numthreads;

GetOptions ("infile=s"  => \$infile,    # input file in fasta format
            "outfile=s" => \$outfile,   # output file with classification
            "numthreads=i" => \$numthreads); # number of simultaneous processes 

my $USAGE = "
Mycofier - A Fungal ITS1 Classifier
This script takes as argument a fasta file (multi fasta) and calculates the 5mer frequency of all the sequences
in the file to produce a classification based on Naive Bayes Classifier.


USAGE:
perl MycofierClassificatorPERL.pl -infile <fasta file> -outfile <output classification> -numthreads <number of concurrent processes>";

unless (-f $infile){
    die "$USAGE\n";}

unless ($outfile){
    die "$USAGE\n";}

unless ($numthreads){
    die "$USAGE\n";}


&trainModel();

my @Q = `grep -o '>' $infile`;
my %hofSeq;
$q->enqueue(@Q);
my @GI;
my @SEQ;

$/="\n>";
open (SEQ, $infile);
while (<SEQ>){
    my $chunk = $_;
    $chunk =~ s/>//g;
    my @array = split /\n/, $chunk;
    my $gi = shift @array;
    my $seq = join ('', @array);
    push (@GI, $gi);
    push (@SEQ, $seq);
}
close SEQ;
$/="\n";
#################################
open (OUT, ">$outfile");
print OUT "#SeqID\tGENUS_CLASSIFICATION\tBOOTSTRAP_VALUE (0-100%) (50 iterations)\n";
my $counter=0;
my (@thr) = map {
    my $line = $_;
    flock(OUT, LOCK_EX) || die "Could not flock file\n";
    threads->create(
	sub {
	    while (defined (my $item = $q->dequeue_nb())) {
		my $gi = $GI[$counter];
		my $seq = $SEQ[$counter];
		my %attributes = &counting ($seq);
		my @frequencyTaxa;  
		for (my $k = 1;$k<=50;$k++){
		    my %bootstrapVector = &bootstrap (%attributes);
		    my $result = $nb->predict
			(attributes => \%bootstrapVector);
		    my %hash = %$result;
		    foreach my $k (sort {$hash{$b} <=> $hash{$a}} keys %hash){
			push (@frequencyTaxa, $k);
			last;
		    }
		}
		my ($mostFrequent, $frequency) = &mostFrequent(@frequencyTaxa);
		$frequency = ($frequency*100)/50;
		print OUT "$gi\t$mostFrequent\t$frequency\n";
		$counter++;
	    }
	}
	)
} 1..$numthreads;

# terminate.
foreach my $thread (@thr){
    $thread->join();
}


###################
### Subroutines ###
###################
sub bootstrap {
    my %bootstrap = @_;
    my %returnHash;
    foreach my $key (keys %bootstrap){
        my $randomNumber = int(rand(5))+1;
        if ($randomNumber == 5){
            $returnHash{$key}=$bootstrap{$key};
        }
    }
    return %returnHash;
}


sub mostFrequent {
    my @items = @_;
    my %count;
    $count{$_}++ for @items;
    my ($winner, $winner_count) = each %count;
    while (my ($maybe, $maybe_count) = each %count) {
        if ($maybe_count > $winner_count) {
            $winner = $maybe;
            $winner_count = $maybe_count;
        }
    }
    return ($winner, $winner_count);
}


sub trainModel {
    $/="\n>";
    open (SEQ, $model);
    while (<SEQ>){
	my $chunk = $_;
	$chunk =~ s/>//g;
	my @array = split "\n", $chunk;
	my $gi = shift @array;
	my @gisplit = split /\s+/, $gi;
	my $genus = $gisplit[1]; 
	my $seq = join ('', @array);
	my %attributes = &counting ($seq);
	$nb->add_instance
	    (attributes => \%attributes,
	     label => $genus);
    }    
    $nb->train;
    $/="\n";
    close SEQ;
}


sub counting {
    my %returnHash;
    my $sequence = $_[0];
    my $GCContent = ($sequence =~ tr/GC//);
    $returnHash{GC}=$GCContent;
    my $loop_size = (length $sequence)-5+1;
    for (my $k=0;$k<$loop_size;$k++){
        my $word = substr ($sequence,$k,5);
        if (exists $returnHash{$word}){
            $returnHash{$word} = $returnHash{$word}+1;
        }
	else {$returnHash{$word}=1;}
    }
    return %returnHash;
}
