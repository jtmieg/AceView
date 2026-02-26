#!bin/tcsh
set run=$1
set dd=tmp/SA/$run 
set toto=$dd/sa2ali.ace

echo "Ali $run\nRun $run" > $toto
echo "-D Letter_profile" >> $toto
echo "Counts\nStrandedness\nAli\nUnicity\nGene_expression\nSponge\nAlignments\nLetter_profile\nATGC_kb\nPair_fate\nErrors\nComputer_ressource" >> $toto
#RawCounts
set ff=$dd/runStats.tsf
cat $ff | gawk -F '\t' "/^$run/{print;}" > $ff.1
cat $ff.1 | gawk -F '\t' -f scripts2/m2.sa2ali.awk >> $toto

set ff=$dd/wiggleCumuls.tsf
if (-e $ff) then
  echo "\nAli $run\nRun $run" >> $toto
  cat $ff | gawk -F '\t' '/^Any/{cumul=$4;exonic=$5;intronic=$6;intergenic=$7;if (exonic+0>0)printf("S_1_exonic %.3f Mb aligned\n", exonic/1000000);if (intronic+0>0)printf("S_1_intronic %.3f Mb aligned\n", intronic/1000000);if (intergenic+0>0)printf("S_1_intergenic %.3f Mb aligned\n", intergenic/1000000);}' >> $toto
endif

set ff=$dd/$run.letterProfile.tsf
if (-e $ff) then
  echo "\nAli $run\nRun $run" >> $toto
  cat $ff | gawk -F '\t' '/^#/{next;}{fr=substr($1,length($1),1);if(fr=="r")fr="f2";else fr="f1"; if($2+0>0)printf("Letter_profile %s %d %d %d %d %d %d 0  %d %d %d %d %d\n",fr,$2,$10,$11,$12,$13,$14,$4,$5,$6,$7,$8)}' >> $toto
endif


set ff=$dd/salign.err
if (-e $ff) then
  echo "\nAli $run\nRun $run" >> $toto
  cat $ff | gawk '/TIMING/{if ($4=="U") {t=$5;}}END{if(t+0>0)printf("CPU sortAlign %d seconds\n", t);}' >> $toto
  cat $ff | gawk '/TIMING/{if ($2=="E") {t=0;x=$3;n=split(x,aa,":");if(n==2)t=60*aa[1]+aa[2];if(n==3)t=3600*aa[1]+60*aa[2]+aa[3];}}END{if(t+0>0)printf("Elapsed sortAlign %d seconds\n", t);}' >> $toto
  cat $ff | gawk '/TIMING/{if ($8=="M") {m=$9;}}END{if(m+0>0)printf("Max_memory sortAlign %d Gb\n", m/1000000);}' >> $toto
  cat $ff | gawk '/TIMING/{if ($10=="P") {p=$11;}}END{if(p+0>0)printf("Multi_threading sortAlign %.2f average_running_threads\n", p/100.0);}' >> $toto
endif

echo >> $toto

set ff=$dd/runStats.tsf
cat $toto | gawk  'BEGIN{s2=100;n=0;}/Stranding/{r=run;s=$3+0;p=$4;m=$6;if(s>50 && p>10000 && s > s1+0){s1=s;n=p;k1=$0;}if(s<50 && m>10000 && s < s2){n=m;s2=s;k2=$0;}}END{if(n>0){if(s1-50 > 50-s2)print k1;else print k2;}}' run=$run | gawk '{n=$4+$6+0;cl=$2;if(n)printf("Run %s\nObserved_strandedness_in_ns_mapping %s %s %d plus %d minus\n\n", run,cl,$3,$4,$6);}' run=$run > $toto.1
cat $toto.1 >> $toto

wc $toto
