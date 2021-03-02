#! /usr/bin/perl -w

# Uses modified to version of matchtract.pl and blastformatProt.pl  PWS 12 May 2020

$usage = "Protexcluder.pl -f bpsofflankingtoremove(default=50)  blastx/blastpfile fastafile \n";

# to exclude the portion matching protein subject in a nucleotide sequence file

if (@ARGV < 2) {die "$usage";}

use Getopt::Std;

getopts("f:");

$Len  = defined $opt_f ? $opt_f : 50;

`rm -f $ARGV[1].ssi`;

`/home/ps997/scripts/ProtExcluder1.1/matchtract.pl $ARGV[0] > $ARGV[0].mt`;

`/home/fay-wei/bin/ProtExcluder1.1/countaanu.pl $ARGV[0].mt > $ARGV[0].mtca`;

`/home/fay-wei/bin/ProtExcluder1.1/rmlowcomplexitymathc.pl  $ARGV[0].mtca 3 60 >  $ARGV[0].mtca_3_60`;

`/home/ps997/scripts/ProtExcluder1.1/blastformatProt.pl $ARGV[0] > $ARGV[0].f`;

`/home/fay-wei/bin/ProtExcluder1.1/rmlowcomfromBF.pl $ARGV[0].mtca_3_60 $ARGV[0].f > $ARGV[0].fnolow`;

`sort -k 6,6 -k 3,3n $ARGV[0].fnolow > $ARGV[0].fnolows`;

`/home/fay-wei/bin/ProtExcluder1.1/mergequeryBF.pl $ARGV[0].fnolows $Len > $ARGV[0].fnolowm50`;

`/home/fay-wei/bin/ProtExcluder1.1/unmatchedregionBF.pl $ARGV[0].fnolowm50 $Len > $ARGV[0].fnolowm50MSP`;

`/home/fay-wei/bin/ProtExcluder1.1/mspesl-sfetch.pl $ARGV[1] $ARGV[0].fnolowm50MSP 0 $ARGV[0].fnolowm50seq`;

`/home/fay-wei/bin/ProtExcluder1.1/mergeunmatchedregion.pl $ARGV[0].fnolowm50seq > $ARGV[0].fnolowm50seqm`;

`/home/fay-wei/bin/ProtExcluder1.1/GCcontent.pl $ARGV[0].fnolowm50seqm > $ARGV[0].fnolowm50seqmGC`;

`/home/fay-wei/bin/ProtExcluder1.1/rmshortseq_noN.pl $ARGV[0].fnolowm50seqmGC $ARGV[0].fnolowm50seqm 50 > $ARGV[0].fnolowm50seqmns`;

`/home/fay-wei/bin/ProtExcluder1.1/getanycolumnuni.pl $ARGV[0].fnolow 6 > $ARGV[0].fnolowlist`;

`/home/fay-wei/bin/ProtExcluder1.1/rmlistedseq.pl $ARGV[0].fnolowlist $ARGV[1] >  $ARGV[1]nPr`;

`cat $ARGV[1]nPr $ARGV[0].fnolowm50seqmns > temp`;

`/home/fay-wei/bin/ProtExcluder1.1/fasta-reformat.pl temp 50 > $ARGV[1]noProtFinal`;
