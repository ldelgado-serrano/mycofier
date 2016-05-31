#!/usr/bin/perl -w
# Mycofier: A new Machine Learning-based classifier for fungal ITS sequences
# Delgado, et al. 2016
# PUBLICACION ISSUE
# Corporacion CorpoGen & Universidad de Los Andes

####### Location Of WEKA jar FILE ########
my $weka = "weka.jar";
my $wekaModel = "dbnosing_sps.model";
my $ARFFFile = "ARFFfile.arff";
my $wekaMemory = "4g"; #Weka needs lots of memory, here we are indicating to use 4 Gigs of RAM
##########################################
use Getopt::Long;

my $infile;
my $outfile;

GetOptions ("infile=s"  => \$infile,    # input file in fasta format
	    "outfile=s" => \$outfile);   # output file with classification

my $USAGE = "
This script takes as argument a fasta file (multi fasta) and calculates the 5mer frequency of all the sequences
in the file to produce an classification based on the Naive Bayes Classification module from WEKA. 
IMPORTANT: Fasta identifiers must be unique for each sequence, with NO spaces.

USAGE:
perl MycofierClassificator.pl -infile <fasta file> -outfile <output classification>";

unless ($infile){
    die "$USAGE\n";}

unless ($outfile){
    die "$USAGE\n";}

#######################
##Check fasta file ####
#######################
my %IDS;
my %fullSeq;
my $averageSeqLen;
my $allSeqLengths;
$/="\n>";
open (FASTA, $infile);
while (<FASTA>){
    my $chunk = $_;
    $chunk =~ s/>//;
    my @ARRAY = split /\n/, $chunk;
    my $id = shift @ARRAY;
    my $sequence = join ('', @ARRAY);
    $sequence = uc($sequence);
    $allSeqLengths .= length ($sequence);
    $fullSeq{$id}= $sequence;
    
    if (exists $IDS{$id}){
	$IDS{$id} = $IDS{$id}+1;
    }
    else {$IDS{$id}=1;}
}
my $numSeqs = keys %IDS;
$averageSeqLen = $allSeqLengths / $numSeqs;

foreach my $element (keys %IDS){
    if ($IDS{$element}>1){
	print "\n\nIdentificators for each sequence are not unique\n $USAGE\n"; exit;
    }
}


##########################
## kmer generator 5-mer ##
##########################
my %hash1; #Keep order of k-mers
my %hash2; #Actual word count
my %hash3; #Total count for all kmers for all sequences
&seq_generator;
my @genusClasses;

##############################
###  Generating ARFF file###
##############################
my $numFeatures = keys (%hash1) + 3;
my $printLine = "% ARFF file for the classification of ITS1 fungal sequences
\%
\@relation ITS1_classifying
\%
\@attribute normalized_length numeric\n";

for my $kmerElement (sort hashValueAscendingNum keys %hash1){
    $printLine .= "\@attribute $kmerElement numeric\n";
}
$printLine .= "\@attribute GCcontent numeric\n";
&addLastHeaderElement;

foreach my $seqID (sort keys %fullSeq){
    
    my $line = &counting ($averageSeqLen, $fullSeq{$seqID});
    $printLine .= $line."\n";
}   
open (OUT, ">$ARFFFile");
print OUT "$printLine\n";
close OUT;

#### Run weka ####
my $wekaOutput = `java -classpath $weka -Xmx$wekaMemory weka.classifiers.bayes.NaiveBayes -l $wekaModel -T $ARFFFile -p 1-1027`;

my @resultingClassificationArray = &processWekaOut($wekaOutput);

my $counter = 0;
open (OUTFILE, ">$outfile");
print OUTFILE "#SeqID\tGENUS_CLASSIFICATION\n";
foreach my $seqID (sort keys %fullSeq){
    print OUTFILE "$seqID\t$resultingClassificationArray[$counter]\n";
    $counter++;
} 
close OUTFILE;

###########################################
##Counting 5 mers in generated dataset   ##
###########################################
sub counting {
    my ($averageSeqLen, $sequence) = ($_[0], $_[1]);
    my $normLen = length ($sequence) / $averageSeqLen;
    my $GCContent = (($sequence =~ tr/GC//) / length ($sequence));
    my $loop_size = (length $sequence)-5+1;
    for (my $k=0;$k<$loop_size;$k++){
	my $word = substr ($sequence,$k,5);
	if (exists $hash2{$word}){
	    $hash2{$word} = $hash2{$word}+1;
	}
    }
    my $line = "$normLen, ";
    foreach my $key (sort hashValueAscendingNum (keys (%hash1))){
	my $freq = $hash2{$key}/length($sequence);
	$line .= "$freq, ";
	$hash2{$key}=0;
    }
    $line .= "$GCContent, ?";
    return $line;
}


#################################
#Sequence Generator Sobroutine ##
#################################
sub seq_generator
{
    my $num_of_sequences = 4**5;
    my $from = keys %hash1;

    my $to = $from + $num_of_sequences;
    #print "$from\t$to\n";
    for (my $i=$from;$i<$to;$i++)
    {
        my $NumericBaseData = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        my $OutputValue ="7"; #A number out of base, or it will complain                                                                                                                
        my $OutputBase =4;
        my $DecimalValue = $i-$from;
        while ($DecimalValue > 0)
        {
            my $X = int((($DecimalValue/$OutputBase)-int($DecimalValue/$OutputBase))*$OutputBase+1.5);
            $OutputValue = substr($NumericBaseData,$X-1,1).$OutputValue;
            $DecimalValue = int($DecimalValue/$OutputBase);
        }
        my $result = $OutputValue;
        my $temp_seq_size = 5+1;
        my $num = sprintf "%${temp_seq_size}.0d",$result;
        $num =~ s/7//g; #Get rid of the number                                                                                                                                          
        $num =~ tr/ 0123/AACGT/;
        # print $num,"\n";                                                                                                                                                              
        $hash1{$num}=$i;
        $hash2{$num}=0;
        $hash3{$num}=[];

    }
}

sub processWekaOut {
    my $wekaOut = $_[0];
    my @lines = split /\n/, $wekaOut;
    #print "$wekaOut\n";
    my @returnArray;
    foreach my $line (@lines){
	$line =~ s/^\s+//g;
	if ($line =~ /\d+\s/){
	    my ($instance, $actual, $predicted, $error, $vector) = split /\s+/, $line;
	    $vector =~ s/\)|\(//g;
	    #my @vectorArray = split (/\,/,$vector);
	    #my @vectorOfZeros = grep {$_ == 0} @vectorArray;
	    #my $sizeOfVector = @vectorArray;
	    #my $sizeOfZeros = @vectorOfZeros;
	    #my $vote = $sizeOfVector - $sizeOfZeros;
	    #print "$instance\t$vote\n";
	    #Example predicted: 164:Hypocr
	    my ($predicted1, $predicted2) = split /:/,$predicted;
	    --$predicted1;
	    push (@returnArray, $genusClasses[$predicted1]); 
	}
    }
    return @returnArray;
}




sub hashValueAscendingNum{
    $hash1{$a} <=> $hash1{$b};
}


sub addLastHeaderElement{
    $printLine .= "\@attribute class { Acremonium, Calonectria, Gymnopus, Hymenopellis, Trichothecium, Saccharomyces, Conidiobolus, Bryoria, Trichosporon, Herpotrichia, Terriera, Trichaptum, Pyrenophora, Terfezia, Ambispora, Umbelopsis, Dermocybe, Erysiphe, Pucciniastrum, Eremothecium, Magnaporthe, Muscodor, Amanita, Flammulina, Peyronellaea, Laetiporus, Lasiodiplodia, Ilyonectria, Blumeria, Rhizophagus, Heterodermia, Phlebia, Fomitiporia, Pycnoporus, Spongipellis, Hericium, Scleroderma, Tylospora, Uwebraunia, Peniophorella, Nigrospora, Serpula, Postia, Xylaria, Arthroderma, Septoria, Tranzschelia, Cladosporium, Lewia, Agaricus, Oidiodendron, Leptosphaeria, Eutypa, Astraeus, Hamigera, Actinomucor, Pisolithus, Emmonsia, Golovinomyces, Hydnellum, Stereum, Melanogaster, Monilinia, Roccella, Gliocladium, Microbotryum, Macowanites, Curvularia, Marasmius, Rhizoscyphus, Gyalolechia, Cryptosporiopsis, Cetraria, Ampelomyces, Volvariella, Alectoria, Phloeospora, Peziza, Dioszegia, Discula, Flavoparmelia, Ceratobasidium, Phellodon, Metacordyceps, Clonostachys, Diatrype, Puccinia, Rhizophlyctis, Botryotinia, Pleospora, Sphaeropsis, Thamnolia, Otidea, Cenococcum, Ceratocystiopsis, Ramaria, Rhizopogon, Epichloe, Funneliformis, Auxarthron, Pachyphloeus, Leucophleps, Monascus, Evernia, Parastagonospora, Hypoxylon, Myrothecium, Coprinopsis, Parmelia, Rhizoplaca, Hydnum, Thanatephorus, Gloeophyllum, Hygrocybe, Geosmithia, Endothia, Acarospora, Protoblastenia, Ulocladium, Phaeoacremonium, Verticillium, Pseudevernia, Rusavskia, Phaeophyscia, Ramulispora, Lobaria, Gloeoporus, Lactarius, Allomyces, Pyrenochaeta, Capronia, Gremmeniella, Asterostroma, Neofusicoccum, Lophodermium, Wickerhamomyces, Haloguignardia, Punctelia, Hypotrachyna, Cerrena, Fusicoccum, Collybia, Euoidium, Mycena, Bjerkandera, Cavernularia, Psilocybe, Phialophora, Antrodia, Neocallimastix, Xanthoparmelia, Polyporus, Camillea, Boeremia, Geopora, Dendrographa, Candida, Lyophyllum, Ceratocystis, Schizophyllum, Fonsecaea, Caloplaca, Dictyonema, Hypocrea, Talaromyces, Moniliophthora, Protoglossum, Russula, Taiwanofungus, Physcia, Mycoblastus, Gaeumannomyces, Apiosporina, Daldinia, Gymnopilus, Cylindrocladium, Melanohalea, Mucor, Phellopilus, Panellus, Arthrobotrys, Agrocybe, Umbilicaria, Sporothrix, Parmeliopsis, Lentinus, Fusarium, Chaetosartorya, Biscogniauxia, Trametes, Cetrelia, Lentinula, Thelephora, Rhodotorula, Gaertneriomyces, Lecanora, Meliniomyces, Cystoderma, Syncephalastrum, Sparassis, Lasiosphaeris, Leotia, Galerina, Hypocrella, Derxomyces, Amyloporia, Alternaria, Merimbla, Simplicillium, Roccellina, Pseudallescheria, Therrya, Cryptococcus, Zasmidium, Paraconiothyrium, Aureobasidium, Sporobolomyces, Cercospora, Cystodermella, Paecilomyces, Rhodosporidium, Gymnascella, Liberomyces, Porodaedalea, Grosmannia, Neonectria, Stachybotrys, Hyphoderma, Geomyces, Sydowia, Chaetosphaeria, Hypogymnia, Phanerochaete, Tolypocladium, Tulasnella, Rhizomucor, Myrmecridium, Diversispora, Pyxine, Placopyrenium, Malassezia, Neosartorya, Acaulospora, Heterobasidion, Auricularia, Tuber, Metarhizium, Lipomyces, Typhula, Flavopunctelia, Rinodina, Entoleuca, Zygosaccharomyces, Stereocaulon, Hypoderma, Cytospora, Dactylellina, Phialocephala, Stagonosporopsis, Physconia, Purpureocillium, Pseudocyphellaria, Glomus, Pichia, Diplodia, Cochliobolus, Flavoplaca, Epicoccum, Fellomyces, Akanthomyces, Eutypella, Coprinellus, Sporodictyon, Xanthoria, Calvitimela, Paxillus, Macrophomina, Flavocetraria, Thysanophora, Inonotus, Podaxis, Megacollybia, Oidium, Hannaella, Fuscoporia, Hortaea, Zoophthora, Bipolaris, Chalara, Verrucaria, Ceratorhiza, Conocybe, Pseudocercospora, Clavulina, Nectria, Lecidea, Rhizoctonia, Phyllactinia, Corticium, Laccaria, Kluyveromyces, Nakazawaea, Sporisorium, Dermatocarpon, Setosphaeria, Xerocomus, Melampsora, Cladonia, Hirsutella, Irpex, Suillus, Tomentellopsis, Biatora, Cladophialophora, Neoscytalidium, Hymenoscyphus, Phoma, Paraphoma, Phacidiopycnis, Annulohypoxylon, Xanthomendoza, Oudemansiella, Botryosphaeria, Sphaerulina, Stemphylium, Fusicladium, Neofabraea, Colletotrichum, Amphinema, Rhodocollybia, Rasamsonia, Hebeloma, Pneumocystis, Anaptychia, Plectosphaerella, Uromyces, Phellorinia, Saccharata, Mycosphaerella, Amylomyces, Nephroma, Lepraria, Imshaugia, Starmerella, Cadophora, Sphaerophorus, Phakopsora, Calvatia, Sclerotinia, Orbilia, Armillaria, Tomentella, Alpova, Trapeliopsis, Pseudoplagiostoma, Sawadaea, Penicillium, Choanephora, Grifola, Eurotium, Togninia, Fomes, Tilletia, Nakaseomyces, Leohumicola, Lichtheimia, Parasola, Meyerozyma, Lachancea, Dufourea, Henrica, Coniosporium, Sarocladium, Dasyspora, Cystocoleus, Lecanicillium, Morchella, Ephelis, Peniophora, Coriolopsis, Aspergillus, Corynespora, Piloderma, Paraglomus, Phymatotrichopsis, Vermispora, Fibroporia, Coccomyces, Pleurotus, Coprinus, Wilcoxina, Inocybe, Byssochlamys, Beauveria, Leveillula, Phomopsis, Claviceps, Hymenogaster, Athelia, Glomerella, Xerula, Cyphellophora, Lentinellus, Hemileia, Teratosphaeria, Harknessia, Ophiostoma, Leptographium, Drechslera, Rosellinia, Dothiorella, Paracoccidioides, Resinicium, Ophiocordyceps, Exophiala, Artomyces, Dendryphion, Connopus, Scedosporium, Absidia, Fomitopsis, Polycauliona, Eudarluca, Claroideoglomus, Coniophora, Valsa, Cookeina, Apiognomonia, Ogataea, Melanelixia, Boletus, Cylindrocladiella, Ganoderma, Calogaya, Phellinus, Pleopsidium, Debaryomyces, Diaporthe, Phlebiopsis, Hanseniaspora, Amandinea, Hypholoma, Cryphonectria, Hypomyces, Tephromela, Dissoconium, Lasiosphaeria, Cladia, Isaria, Nemania, Aspicilia, Gymnomyces, Alnicola, Scheffersomyces, Marssonina, Chaenotheca, Cronartium, Scytalidium, Emericella, Ajellomyces, Botrytis, Venturia, Podosphaera, Trichophyton, Waitea, Hyphodermella, Gibellulopsis, Mucidula, Pestalotiopsis, Coccocarpia, Microdochium, Ramalina, Nematoctonus, Gnomonia, Cunninghamella, Hypsizygus, Phillipsia, Parmotrema, Chrysoporthe, Lepiota, Leucoagaricus, Strobilurus, Physciella, Pleurostomophora, Clitopilus, Usnea, Filobasidiella, Chaetomium, Microthia, Blakeslea, Cortinarius, Peltigera, Sarcinomyces, Microsporum, Hyperphyscia, Pilobolus, Rhizopus, Batcheloromyces, Cordyceps, Melanoleuca, Neoerysiphe, Collophora, Mycocalicium, Guignardia, Tremella, Sebacina, Xanthopsoroma, Battarrea, Didymella, Tricholoma }
\%
\@data
\%\n";
    @genusClasses = split (/,\s/, "Acremonium, Calonectria, Gymnopus, Hymenopellis, Trichothecium, Saccharomyces, Conidiobolus, Bryoria, Trichosporon, Herpotrichia, Terriera, Trichaptum, Pyrenophora, Terfezia, Ambispora, Umbelopsis, Dermocybe, Erysiphe, Pucciniastrum, Eremothecium, Magnaporthe, Muscodor, Amanita, Flammulina, Peyronellaea, Laetiporus, Lasiodiplodia, Ilyonectria, Blumeria, Rhizophagus, Heterodermia, Phlebia, Fomitiporia, Pycnoporus, Spongipellis, Hericium, Scleroderma, Tylospora, Uwebraunia, Peniophorella, Nigrospora, Serpula, Postia, Xylaria, Arthroderma, Septoria, Tranzschelia, Cladosporium, Lewia, Agaricus, Oidiodendron, Leptosphaeria, Eutypa, Astraeus, Hamigera, Actinomucor, Pisolithus, Emmonsia, Golovinomyces, Hydnellum, Stereum, Melanogaster, Monilinia, Roccella, Gliocladium, Microbotryum, Macowanites, Curvularia, Marasmius, Rhizoscyphus, Gyalolechia, Cryptosporiopsis, Cetraria, Ampelomyces, Volvariella, Alectoria, Phloeospora, Peziza, Dioszegia, Discula, Flavoparmelia, Ceratobasidium, Phellodon, Metacordyceps, Clonostachys, Diatrype, Puccinia, Rhizophlyctis, Botryotinia, Pleospora, Sphaeropsis, Thamnolia, Otidea, Cenococcum, Ceratocystiopsis, Ramaria, Rhizopogon, Epichloe, Funneliformis, Auxarthron, Pachyphloeus, Leucophleps, Monascus, Evernia, Parastagonospora, Hypoxylon, Myrothecium, Coprinopsis, Parmelia, Rhizoplaca, Hydnum, Thanatephorus, Gloeophyllum, Hygrocybe, Geosmithia, Endothia, Acarospora, Protoblastenia, Ulocladium, Phaeoacremonium, Verticillium, Pseudevernia, Rusavskia, Phaeophyscia, Ramulispora, Lobaria, Gloeoporus, Lactarius, Allomyces, Pyrenochaeta, Capronia, Gremmeniella, Asterostroma, Neofusicoccum, Lophodermium, Wickerhamomyces, Haloguignardia, Punctelia, Hypotrachyna, Cerrena, Fusicoccum, Collybia, Euoidium, Mycena, Bjerkandera, Cavernularia, Psilocybe, Phialophora, Antrodia, Neocallimastix, Xanthoparmelia, Polyporus, Camillea, Boeremia, Geopora, Dendrographa, Candida, Lyophyllum, Ceratocystis, Schizophyllum, Fonsecaea, Caloplaca, Dictyonema, Hypocrea, Talaromyces, Moniliophthora, Protoglossum, Russula, Taiwanofungus, Physcia, Mycoblastus, Gaeumannomyces, Apiosporina, Daldinia, Gymnopilus, Cylindrocladium, Melanohalea, Mucor, Phellopilus, Panellus, Arthrobotrys, Agrocybe, Umbilicaria, Sporothrix, Parmeliopsis, Lentinus, Fusarium, Chaetosartorya, Biscogniauxia, Trametes, Cetrelia, Lentinula, Thelephora, Rhodotorula, Gaertneriomyces, Lecanora, Meliniomyces, Cystoderma, Syncephalastrum, Sparassis, Lasiosphaeris, Leotia, Galerina, Hypocrella, Derxomyces, Amyloporia, Alternaria, Merimbla, Simplicillium, Roccellina, Pseudallescheria, Therrya, Cryptococcus, Zasmidium, Paraconiothyrium, Aureobasidium, Sporobolomyces, Cercospora, Cystodermella, Paecilomyces, Rhodosporidium, Gymnascella, Liberomyces, Porodaedalea, Grosmannia, Neonectria, Stachybotrys, Hyphoderma, Geomyces, Sydowia, Chaetosphaeria, Hypogymnia, Phanerochaete, Tolypocladium, Tulasnella, Rhizomucor, Myrmecridium, Diversispora, Pyxine, Placopyrenium, Malassezia, Neosartorya, Acaulospora, Heterobasidion, Auricularia, Tuber, Metarhizium, Lipomyces, Typhula, Flavopunctelia, Rinodina, Entoleuca, Zygosaccharomyces, Stereocaulon, Hypoderma, Cytospora, Dactylellina, Phialocephala, Stagonosporopsis, Physconia, Purpureocillium, Pseudocyphellaria, Glomus, Pichia, Diplodia, Cochliobolus, Flavoplaca, Epicoccum, Fellomyces, Akanthomyces, Eutypella, Coprinellus, Sporodictyon, Xanthoria, Calvitimela, Paxillus, Macrophomina, Flavocetraria, Thysanophora, Inonotus, Podaxis, Megacollybia, Oidium, Hannaella, Fuscoporia, Hortaea, Zoophthora, Bipolaris, Chalara, Verrucaria, Ceratorhiza, Conocybe, Pseudocercospora, Clavulina, Nectria, Lecidea, Rhizoctonia, Phyllactinia, Corticium, Laccaria, Kluyveromyces, Nakazawaea, Sporisorium, Dermatocarpon, Setosphaeria, Xerocomus, Melampsora, Cladonia, Hirsutella, Irpex, Suillus, Tomentellopsis, Biatora, Cladophialophora, Neoscytalidium, Hymenoscyphus, Phoma, Paraphoma, Phacidiopycnis, Annulohypoxylon, Xanthomendoza, Oudemansiella, Botryosphaeria, Sphaerulina, Stemphylium, Fusicladium, Neofabraea, Colletotrichum, Amphinema, Rhodocollybia, Rasamsonia, Hebeloma, Pneumocystis, Anaptychia, Plectosphaerella, Uromyces, Phellorinia, Saccharata, Mycosphaerella, Amylomyces, Nephroma, Lepraria, Imshaugia, Starmerella, Cadophora, Sphaerophorus, Phakopsora, Calvatia, Sclerotinia, Orbilia, Armillaria, Tomentella, Alpova, Trapeliopsis, Pseudoplagiostoma, Sawadaea, Penicillium, Choanephora, Grifola, Eurotium, Togninia, Fomes, Tilletia, Nakaseomyces, Leohumicola, Lichtheimia, Parasola, Meyerozyma, Lachancea, Dufourea, Henrica, Coniosporium, Sarocladium, Dasyspora, Cystocoleus, Lecanicillium, Morchella, Ephelis, Peniophora, Coriolopsis, Aspergillus, Corynespora, Piloderma, Paraglomus, Phymatotrichopsis, Vermispora, Fibroporia, Coccomyces, Pleurotus, Coprinus, Wilcoxina, Inocybe, Byssochlamys, Beauveria, Leveillula, Phomopsis, Claviceps, Hymenogaster, Athelia, Glomerella, Xerula, Cyphellophora, Lentinellus, Hemileia, Teratosphaeria, Harknessia, Ophiostoma, Leptographium, Drechslera, Rosellinia, Dothiorella, Paracoccidioides, Resinicium, Ophiocordyceps, Exophiala, Artomyces, Dendryphion, Connopus, Scedosporium, Absidia, Fomitopsis, Polycauliona, Eudarluca, Claroideoglomus, Coniophora, Valsa, Cookeina, Apiognomonia, Ogataea, Melanelixia, Boletus, Cylindrocladiella, Ganoderma, Calogaya, Phellinus, Pleopsidium, Debaryomyces, Diaporthe, Phlebiopsis, Hanseniaspora, Amandinea, Hypholoma, Cryphonectria, Hypomyces, Tephromela, Dissoconium, Lasiosphaeria, Cladia, Isaria, Nemania, Aspicilia, Gymnomyces, Alnicola, Scheffersomyces, Marssonina, Chaenotheca, Cronartium, Scytalidium, Emericella, Ajellomyces, Botrytis, Venturia, Podosphaera, Trichophyton, Waitea, Hyphodermella, Gibellulopsis, Mucidula, Pestalotiopsis, Coccocarpia, Microdochium, Ramalina, Nematoctonus, Gnomonia, Cunninghamella, Hypsizygus, Phillipsia, Parmotrema, Chrysoporthe, Lepiota, Leucoagaricus, Strobilurus, Physciella, Pleurostomophora, Clitopilus, Usnea, Filobasidiella, Chaetomium, Microthia, Blakeslea, Cortinarius, Peltigera, Sarcinomyces, Microsporum, Hyperphyscia, Pilobolus, Rhizopus, Batcheloromyces, Cordyceps, Melanoleuca, Neoerysiphe, Collophora, Mycocalicium, Guignardia, Tremella, Sebacina, Xanthopsoroma, Battarrea, Didymella, Tricholoma");
}
