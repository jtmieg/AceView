#!bin/tcsh
set run=$1
set dd=$2 # equivalent to sortaling -o $dd
set toto=$dd/sa2ali.ace

echo "Ali m2.$run\nRun2 $run" > $toto
#RawCounts
set ff=$dd/runStats.tsf
cat $ff | gawk -F '\t' "/^$run/{print;}" > $ff.1
cat $ff.1 | gawk -F '\t' -f scripts2/m2.sa2ali.awk >> $toto

set ff=$dd/wiggleCumuls.tsf
if (-e $ff) then
  cat $ff | gawk -F '\t' '/^Any/{cumul=$4;exonic=$5;intronic=$6;intergenic=$7;if (exonic+0>0)printf("S_1_exonic %.3f Mb aligned\n", exonic/1000000);if (intronic+0>0)printf("S_1_intronic %.3f Mb aligned\n", intronic/1000000);if (intergenic+0>0)printf("S_1_intergenic %.3f Mb aligned\n", intergenic/1000000);}' >> $toto
endif

set ff=$dd.err
if (-e $ff) then
  cat $ff | gawk '/TIMING/{if ($4=="U") {x=$5;n=split(x,aa,":");if(n==2)t=60*aa[1]+aa[2];if(n==3)t=3600*aa[1]+60*aa[2]+aa[3];if(t0+0=0 || t<t0)t0=t;}}END{printf("CPU %d seconds\n", t0);}' >> $toto
  cat $ff | gawk '/TIMING/{if ($4=="E") {x=$5;n=split(x,aa,":");if(n==2)t=60*aa[1]+aa[2];if(n==3)t=3600*aa[1]+60*aa[2]+aa[3];if(t0+0=0 || t<t0)t0=t;}}END{printf("Elapsed %d seconds\n", t0);}' >> $toto
  cat $ff | gawk '/TIMING/{if ($10=="M") {m=$11;printf("Max_memory sortAlign %s Mb\n", m);exit;}' >> $toto
endif

echo >> $toto
wc $toto
