MYCOFIER 1.0
1th January 2016
Java implementation of Mycofier.
It runs by running a perl script wrapper around weka.
Once you clone the entire project from git you will end up with the folder "mycofier-master"
Weka jar file, perl script, file model (dbnosing_sps.model) must be in this folder
for it to work
Mycofier requires Java 1.5 or higher and perl.

USAGE:
perl MycofierClassificator.pl -infile <fasta file> -outfile <output classification> -numthreads <number of concurrent processes>

DISCLAIMER:
Unfortunately weka is too slow to make it useful to classify sequence data. We recommend the PERL implementation of the classifier, which is way faster.
