#!bin/tcsh -f

set MAGIC=$1
set chrom=$2

set minIntronSupport=`cat MetaDB/$MAGIC/RunList | gawk '{n++}END{printf("%d",int(n/20)+1);}'`
if (-e  tmp/introns/d4.$MAGIC.candidate_introns.ace) then
  set n=`cat tmp/introns/d4.$MAGIC.candidate_introns.ace | gawk '/^Ali/{ok=1;}/^Candidate_introns/{if(ok==1)n=$19+0;last;}END{print n+0}' `
  if ($n >0) set minIntronSupport=$n
endif

set sMin="-minimalIntronSupport $minIntronSupport"
echo "Minimal Intron support sMin=$sMin"
if (! -d  tmp/X.$MAGIC/EHITS/$chrom) mkdir  tmp/X.$MAGIC/EHITS/$chrom
set sDuo=""

goto deUno

set sDuo=" -newDoubleIntrons "
if (! -e tmp/X.$MAGIC/EHITS/$chrom/introns.de_duo) then
  echo "collating de_duo"

    # touch  tmp/X.$MAGIC/EHITS/any.$chrom.limit1
    foreach run (`cat MetaDB/$MAGIC/RunList`)
        # was OR pA pT SL1 SL2 SL3 SL4 SL5 SL6 SL7 SL8 SL9 SL10 SL11 SL12, 2016_04_01 we now prefer to do the XA in phase f3
       foreach dd (OR)
         echo $run $chrom $dd
         set ok=0
         foreach n (1 2 5 10 20 50 100 200 1000)
           if ($ok == 1) continue
           if (-e tmp/$dd/$run/$chrom/$run.$n.0.txt) then
             set ok=1
             cat  tmp/$dd/$run/$chrom/$run.$n.*.txt | gawk -F '\t' '{if($3 != chrom)next;}/^pA/{print}/^pT/{print}/^SL/{print}/^INTRON/{print}'  chrom=$chrom |  sort   >>  tmp/X.$MAGIC/EHITS/$chrom/introns.de_duo
           endif
         end
       end
    end

endif
  
set  iDuo=tmp/X.$MAGIC/EHITS/$chrom/introns.de_duo
if (-e $iDuo) then
  set sDuo="-sxxNewIntronsFileName  $iDuo "
  ls -ls $iDuo
endif


deUno:

if (! -e tmp/X.$MAGIC/EHITS/$chrom/f1.introns.de_uno.gz) then
  echo "collating de_uno"
  if (-e tmp/X.$MAGIC/EHITS/$chrom/f1.introns.de_uno.1) \rm tmp/X.$MAGIC/EHITS/$chrom/f1.introns.de_uno.*
  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    gunzip -c tmp/OR/$run/d4.de_uno.txt.gz | gawk -F '\t' '{if($1==chrom)print;}' chrom=$chrom >> tmp/X.$MAGIC/EHITS/$chrom/f1.introns.de_uno.1
  end

  cat tmp/X.$MAGIC/EHITS/$chrom/f1.introns.de_uno.1 | gawk -F '\t' '{z=$1 "\t" $2 "\t" $3 ; n[z] += $4;if (length($5)>1) t[z]=$5;}END{for (k in n) printf("%s\t%d\t%s\n",k,n[k],t[z]);}' >  tmp/X.$MAGIC/EHITS/$chrom/f1.introns.de_uno.2
  cat  tmp/X.$MAGIC/EHITS/$chrom/f1.introns.de_uno.2 | gawk  -F '\t' '{chrom=$1;a1=$2;a2=$3;nn=$4;type=$5;ln=a2-a1;if(ln<0)ln=-ln;ln+=1;if(a1<a2){s="Forward";b1=a1-35;b2=a1-1;c1=a2+1;c2=a2+35;}else{s="Reverse";b1=a1+35;b2=a1+1;c1=a2-1;c2=a2-35;}printf("INTRON\t%s\t%s\t%d\t%d\t%s\t%d\t%d\t\t%d\t%d\n",s,chrom,b1,b2,chrom,c1,c2,ln,nn);}' >  tmp/X.$MAGIC/EHITS/$chrom/f1.introns.de_uno

  \rm  tmp/X.$MAGIC/EHITS/$chrom/f1.introns.de_uno.[12]
  gzip tmp/X.$MAGIC/EHITS/$chrom/f1.introns.de_uno
endif

set iUno=tmp/X.$MAGIC/EHITS/$chrom/f1.introns.de_uno.gz

set sUno=""

if (-e $iUno) then
  set sUno="-sxxDeUnoIntronsFileName  $iUno"
  ls -ls $iUno
endif

set chr=""
# if ($species == worm)   set chr="CHROMOSOME_"

if (! -e tmp/X.$MAGIC/EHITS/$chrom/f1.txts) then

    echo "running bin/geneelements"
    echo " bin/geneelements -newDoubleIntrons -stranded $sUno $sDuo $sMin -sxxChromosome $chr$chrom -t TARGET/CHROMS/$species.chrom_$chrom.fasta.gz"
           bin/geneelements -newDoubleIntrons -stranded $sUno $sDuo $sMin -sxxChromosome $chr$chrom -t TARGET/CHROMS/$species.chrom_$chrom.fasta.gz  | sort -u  > tmp/X.$MAGIC/EHITS/$chrom/f1.txts
    ls -ls  tmp/X.$MAGIC/EHITS/$chrom/f1.txts

endif


echo "collate ace       "
# intron name $8 $9 was u1 u2   2024_02 20
  cat   tmp/X.$MAGIC/EHITS/$chrom/f1.txts | gawk -F '\t' '/^#/{next}/^\/\//{next}{chrom=$1;gsub(/CHROMOSOME_/,"",chrom);u1=$6;u2=$11;dnaa=$12;dnab=$13;dnac=$14;support=$15;feet=$16;if(length(feet)==5){if(feet == "gt_ag" || feet == "gc_ag" || feet == "at_ac" || feet == "ct_ac"){feet=feet  "\n";}else {feet= "Other "  feet  "\n" ;}}else feet="Julot\n" ;typea=$3;typeb=$4;typec=$5; if(typeb == "Exon")next;xx="XI_" ; col="PALEYELLOW";nam=xx group "_" chrom "__" $8 "_" $9 ;names[nam]++; n = names[nam] ;nam = nam "." n ; printf("Sequence %s\ncDNA_clone %s\nIntMap %s %d %d\nIs_Read\nForward\nComposite %d\nColour %s\n%s",nam, nam, $1,u1,u2,support, col,feet) ;sl=0;if(substr(typea,1,2)=="SL"){x1=length(dnaa);sl=substr(typea,3,1)+0;}if(substr(typeb,1,2)=="SL"){x1=length(dnaa);sl=substr(typeb,3,1)+0;};if(sl>0)printf("Transpliced_to SL%d %d\n",sl,x1);pA=0;if(typeb=="pA"){x1=length(dnaa);}if(typec=="pA"){x1=length(dnaa)+length(dnab);pA=length(dnac);};if(pA>0)printf("PolyA_after_base %d\n",x1);if(typeb=="Intron")printf("Intron %s__%d_%d\n",chrom,$8,$9);if($2=="Reverse")dx=-1;else dx=1;if(typea=="Intron")printf("Intron %s__%d_%d\n",chrom,$7+dx,$8-dx);if(typec=="Intron")printf("Intron %s__%d_%d\n",chrom,$9+dx,$10-dx);printf("\n");}' group=$MAGIC  > tmp/X.$MAGIC/EHITS/$chrom/f1.ace
 
echo "collate fasta"

cat   tmp/X.$MAGIC/EHITS/$chrom/f1.txts  | gawk -F '\t' '/^#/{next}/^\/\//{next}{chrom=$1;gsub(/CHROMOSOME_/,"",chrom);u1=$6;u2=$11;typeb=$4;if(typeb == "Exon") { next;xx="XJ_" ;} else {xx="XI_" ; }nam=xx group "_" chrom "__" u1 "_" u2 ; names[nam]++; n = names[nam] ;nam = nam "." n ; printf(">%s\n%s%s%s\n",nam, $12,$13,$14) ;}' group=$MAGIC  > tmp/X.$MAGIC/EHITS/$chrom/f1.fasta

gzip  tmp/X.$MAGIC/EHITS/$chrom/f1.ace
gzip  tmp/X.$MAGIC/EHITS/$chrom/f1.fasta
gzip  tmp/X.$MAGIC/EHITS/$chrom/f1.txts
touch tmp/X.$MAGIC/EHITS/$chrom/f1.done

exit 0

