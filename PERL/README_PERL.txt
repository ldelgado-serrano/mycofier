MYCOFIER 1.0
1th January 2016
PERL implementation of Mycofier.
For it to work you need the following:

File model: dbnosing_sps.model.fasta
File with Fungal ITS1 sequences to be classified (here is where your own sequences go): ExampleITS1.fasta 


USAGE:
perl MycofierClassificatorPERL.pl -infile <fasta file> -outfile <output classification> -numthreads <number of concurrent processes>



PERL Library requirements:
Getopt::Long
Algorithm::NaiveBayes
threads
Thread::Queue

Which can be installed from the command line as follows:

> sudo perl -MCPAN -e shell
cpan> install Getopt::Long
cpan> install Algorithm:NaiveBayes
cpan> install threads
cpan> install Thread::Queue
