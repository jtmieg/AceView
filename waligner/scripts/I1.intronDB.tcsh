#!bin/tcsh -f

set chrom=$1
set phases="$2"

# This loop ends in phaseLoop
foreach phase ($phases)

echo "phase=$phase  project=$MAGIC"

if (0 && $chrom != 20) exit 0
if ($phase == chromDB) goto chromDB
if ($phase == parseGenome) goto parseGenome
if ($phase == parseTargetIntrons) goto parseTargetIntrons
if ($phase == IntronFeet) goto IntronFeet
if ($phase == Donors) goto Donors
if ($phase == deUno) goto deUno
if ($phase == deMrna) goto deMrna

if ($phase == deDuo) goto deDuo
if ($phase == collate) goto collate
if ($phase == donorAcceptor) goto donorAcceptor
if ($phase == diffAcceptor) goto diffAcceptor
if ($phase == capture) goto capture
if ($phase == intronComfirmation) goto  intronComfirmation

echo "unknown phase $phase"
exit 1
goto phaseLoop

###############################################################
chromDB:
echo -n "I1 phase $phase start :"
date

if (! -d tmp/INTRON_DB)   mkdir tmp/INTRON_DB
if (! -d tmp/INTRON_DB/$chrom)   mkdir tmp/INTRON_DB/$chrom
if (! -d tmp/INTRON_DB/$chrom/database) then
    pushd tmp/INTRON_DB/$chrom
    mkdir database
    ln -s ../../../metaData/wspec.aceview_web_site wspec
    tace . <<EOF
y
EOF
    popd
endif
goto phaseLoop

###############################################################
### parse the genome 
parseGenome:
echo -n "I1 phase $phase start :"
date

if (! -e tmp/INTRON_DB/$chrom/I1.parse_genome.done) then
  pushd tmp/INTRON_DB/$chrom
    ../../../bin/tace . <<EOF
      parse ../../../TARGET/CHROMS/$species.chrom_$chrom.fasta.gz
      find sequence $chrom
      edit genomic
      parse ../../../TARGET/CHROMS/$species.mito.fasta.gz
      save
      quit
EOF
    if (-e I1.intron_feet.done) \rm I1.intron_feet.done
    touch I1.parse_genome.done
  popd
endif

goto phaseLoop

###############################################################
### parse the gene, mrna, intron relations  : intron->supports

parseTargetIntrons:
echo -n "I1 phase $phase start :"
date

# parse in INTRON_DB/$chrom the introns/genes/mRNA/chromosomes

if (! -e tmp/INTRON_DB/$chrom/I1.parse_genes.done) then  
  zcat tmp/METADATA/gtf.*.introns.gz | gawk -F '\t' '{type=substr($1,4);c=$6;if(c != chrom)next;i1=$7+0;i2=$8+0;ln=i2-i1;;if(ln<0)ln=-ln; ln=ln+1;printf("Intron %s__%d_%d\nIntMap %s %d %d\nLength %d\n%s\nGene \"%s\"\nIn_mRNA \"%s\" %d %d\n\n",c,i1,i2,c,i1,i2,ln,type,$2,$3,$4,$5);}' chrom=$chrom | gzip >  tmp/INTRON_DB/$chrom/I1.TargetIntrons.ace.gz
  pushd tmp/INTRON_DB/$chrom
    ../../../bin/tace . <<EOF
      read-models
      parse I1.TargetIntrons.ace.gz
      save
      quit
EOF
    touch I1.parseTargetIntrons.done
    if (-e I1.IntronFeet.done) \rm I1.IntronFeet.done
  popd
endif

goto phaseLoop

########################################################################
### check we have all intron length feet intmap

IntronFeet:
echo -n "I1 phase $phase start :"
date

if (-d tmp/INTRON_DB/$chrom/database && ! -e tmp/INTRON_DB/$chrom/I1.IntronFeet.done) then
  bin/altintrons --setFeet --db tmp/INTRON_DB/$chrom -p $MAGIC
  touch tmp/INTRON_DB/$chrom/I1.IntronFeet.done
endif

goto phaseLoop

########################################################################
### create the donor/accptors of all introns

Donors:
echo -n "I1 phase $phase start :"
date

if (-d tmp/INTRON_DB/$chrom/database && ! -e tmp/INTRON_DB/$chrom/I1.Donors.done) then
  bin/altintrons --setDA --db tmp/INTRON_DB/$chrom -p $MAGIC
  touch tmp/INTRON_DB/$chrom/I1.Donors.done
endif

goto phaseLoop

########################################################################
### parse the genomic de_uno counts

deUno:
echo -n "I1 phase $phase start :"
date

set ok=0
set ff=tmp/INTRON_DB/$chrom/I2.$MAGIC.deUno
echo ' ' > $ff

foreach run1 (`cat MetaDB/$MAGIC/RunList`) 
  set run2=$run1
  set run2=`cat MetaDB/$MAGIC/r2sublib | gawk -F '\t' '{if($2==run1){r=$1;last;}}END{if(r)print r;else print run1;}' run1=$run1`
  if (-e tmp/INTRON_DB/$chrom/I2.$run2.deUno.gz) continue
  set ok=1
  if ($USEMAGICBLAST == 1) then
    zcat tmp/MAGICBLAST/$run1/f2.*.introns.tsf.gz | bin/tsf -s $run2 --merge | gawk -F '\t' '{split($1,aa,"__");c=aa[1];if (c != chrom) next;plit(aa[2],bb,"_");a2=bb[1];b1=bb[2];if(a2<b1){strand="Forward";a2--;a1=a2-12;b1++;b2=b1+10;ln=b1-a2-1;}else{strand="Reverse";a2++;a1=a2+10;b1--;b2=b1-10;ln=a2-b1-1;}nr=$4;feet=$7;if(nr>0)printf("INTRON\t%s\t%s\t%09d\t%09d\t%s\t%09d\t%09d\t%s\t%d\t%d\t%d\n",strand,chrom,a1,a2,chrom,b1,b2,feet,ln,nr,nperfect);}' chrom=$chom  | sort | gzip > tmp/INTRON_DB/$chrom/I2.$run2.deUno.gz
  else
    gunzip -c tmp/PHITS_genome/$run1/*.introns.gz  tmp/PHITS_mito/$run1/*.introns.gz  | gawk -F '\t' '/^INTRON/{c=$3;if (c != chrom) next;if ($6 != chrom) next;i1=$5+0;i2=$7+0;if ($2=="Forward"){i1++;i2--;}else {i1--;i2++;}nn=$11; printf("%s__%d_%d\t%s\t%d\t%s\n",c,i1,i2,run,nn,$9);}' run=$run2 chrom=$chrom | gzip > tmp/INTRON_DB/$chrom/I2.$run2.deUno.gz
  endif

  zcat tmp/INTRON_DB/$chrom/I2.$run2.deUno.gz >> $ff
end

if ($ok == 1) then
  date
  cat $ff | sort > $ff.sorted
  date
  /usr/bin/time bin/altintrons --deMrna $ff.sorted --db tmp/INTRON_DB/$chrom -p $MAGIC
  date
  if (-e  tmp/INTRON_DB/$chrom/I1.IntronFeet.done) \rm  tmp/INTRON_DB/$chrom/I1.IntronFeet.done
  \rm $ff $ff.ace
endif

goto phaseLoop

########################################################################
### parse the intron counts from mRNA alignments

deMrna:
echo -n "I1 phase $phase start :"
date

# export the not-merged max counts
if (0) then
  bin/tacembly tmp/INTRON_DB/$chrom << EOF
    query find intron de_duo
    select -o tmp/INTRON_DB/$chrom/i2r2n.txt1 ii, run, n from ii in @, run in ii->de_duo, n in run[1] 
EOF
endif

#gzip tmp/INTRON_DB/$chrom/i2r2n.txt1

set ok=0
set f0=tmp/INTRON_DB/$chrom/I2.deMrna
if (-e $f0) \rm $f0
foreach run (`cat MetaDB/$MAGIC/RunsList`) 
  set f1=tmp/INTRON_DB/$chrom/I2.$run.deMrna
  if (-e $f1.gz) continue

  set ff=tmp/INTRONRUNS/$run/$run.u.intronSupport.counts.gz 
  if (! -e $ff) continue 
  set ok=1
  zcat $ff | gawk  "/$chrom\t/{print}" > $f1
  cat $f1 | gawk '{printf("%s__%d_%d\t%s\t%d\n", $1,$2+0,$3+3, run, $4+0) ;}' run=$run >> $f0
  gzip $f1
end
date

if ($ok == 1) then
  cat $f0 | sort > $f0.sorted
  date
  bin/altintrons --deMrna $f0.sorted --db tmp/INTRON_DB/$chrom -p $MAGIC
  date
  if (-e  tmp/INTRON_DB/$chrom/I1.IntronFeet.done) \rm  tmp/INTRON_DB/$chrom/I1.IntronFeet.done

# export the merged max counts
  bin/tacembly tmp/INTRON_DB/$chrom << EOF
    query find intron de_duo
    select -o tmp/INTRON_DB/$chrom/i2r2n.txt ii, run, n from ii in @, run in ii->de_duo, n in run[1] 
EOF
  gzip tmp/INTRON_DB/$chrom/i2r2n.txt
endif

goto phaseLoop

########################################################################
## 
# altIntrons

altIntrons:
echo -n "I1 phase $phase start :"
date

   echo "pparse MetaDB/$MAGIC/runs.ace" | bin/tace tmp/INTRON_DB/$chrom -no_prompt
   bin/altintrons -db tmp/INTRON_DB/$chrom -p $MAGIC -o tmp/INTRON_DB/$chrom/I5.altIntrons
goto phaseLoop

########################################################################
## 
# count introns export tsf

intronCounts:
echo -n "I1 phase $phase start :"
date

   echo "pparse MetaDB/$MAGIC/runs.ace" | bin/tace tmp/INTRON_DB/$chrom -no_prompt
   bin/altintrons -db tmp/INTRON_DB/$chrom -p $MAGIC -o tmp/INTRON_DB/$chrom/I2.$MAGIC.intronCounts
goto phaseLoop

########################################################################
## deDuoCumul

deDuo:
echo -n "I1 phase $phase start :"
date
set toto = tmp/introns/I1.$MAGIC.de_duo.$chrom

if (! -e tmp/INTRON_DB/$chrom/I1.$MAGIC.deDuo.done) then

  # hack, because of the schema, the supports are set in intron->de_duo rather than de_uno in the final ace file
  set toto = tmp/INTRON_DB/$chrom/I1.$MAGIC.de_duo
  echo ' ' > $toto.1
  echo ' ' > $toto.txt
  foreach run (`cat MetaDB/$MAGIC/RunsList`)
    set ff=tmp/OR/$run/d4.de_uno.txt.gz
    if (! -e $ff) set ff=tmp/OR/$run/I1.de_uno.txt.gz
    if (-e  $ff) then
      set minS=`cat  tmp/introns/I1.$MAGIC.minS | gawk -F '\t' 'BEGIN{n=1;}{if($1==run)n=$2+0;}END{print n}' run=$run`
      zcat $ff | gawk -F '\t' "/^$chrom\t/"'{a1=$2+0 ; a2=$3+0 ; ii=$1 "__" a1 "_" a2 ; if ($4>=minS) printf ( "%s\t%s\t%s\t%s\t%d\n",ii,$5,"any",run,$4 ) ; }' run=$run minS=$minS chrom=$chrom >> $toto.txt
    endif
  end

  cat $toto.txt | sort -V | gzip > $toto.txt.gz
  \rm  $toto.txt  $toto.1 $toto.2.gz
  # hack, because of the schema, the supports are set in intron->de_duo rather than de_uno in the final ace file
  zcat $toto.txt.gz  | gawk -F '\t' '{ii=$1; if($5+0<1)next;if(ii!=old){if (n>0)printf("RNA_seq %d\n",n);n=0;printf("\nIntron \"%s\"\n",ii);}old=ii;printf("de_duo %s %d\n",$4,$5);n+=$5;}END{if (n>0)printf("RNA_seq %d\n",n);printf("\n");}' | gzip > tmp/INTRON_DB/$chrom/I1.$MAGIC.de_duo.ace.gz

  if (-e tmp/INTRON_DB/$chrom/I1.$MAGIC.collate.done) \rm tmp/INTRON_DB/$chrom/I1.$MAGIC.collate.done
  touch tmp/INTRON_DB/$chrom/I1.$MAGIC.deDuo.done
endif
ls -ls tmp/INTRON_DB/$chrom/I1.$MAGIC.de_duo.ace.gz
goto phaseLoop

#############################################################################
## 
# check we have all intron length feet intmap
donorAcceptor:
echo -n "I1 phase $phase start :"
date

if (! -e tmp/INTRON_DB/$chrom/I1.$MAGIC.collate.done) goto phaseLoop

# check donor acceptors
if (! -e tmp/INTRON_DB/$chrom/I1.$MAGIC.introns.DA.done) then
  pushd tmp/INTRON_DB/$chrom
    ../../../bin/tace . <<EOF
      query find intron
      bql -o I1.donor.f.txt select ii,chrom,a1,a2,d1,d2 from ii in @, chrom in ii->intmap, a1 in chrom[1], a2 in chrom[2] where a1<a2, s in OBJECT("Sequence",chrom), d1 in DNA(s,a1-50,a1+49), d2 in DNA(s,a2-49,a2+50)
      undo
      bql -o I1.donor.r.txt select ii,chrom,a1,a2,d1,d2 from ii in @, chrom in ii->intmap, a1 in chrom[1], a2 in chrom[2] where a1>a2, s in OBJECT("Sequence",chrom), d1 in DNA(s,a1+50,a1-49), d2 in DNA(s,a2+49,a2-50)
      quit
EOF
    cat I1.donor.f.txt | gawk -F '\t' '{ii=$1;m=$2;a1=$3;a2=$4;d1=$5;d2=$6;printf("Intron %s\n-D DA\nD %s__%d_f\nA %s__%d_f\nDonor %s\nAcceptor %s\n\n",ii,m,a1-1,m,a2+1,d1,d2);}' | grep -v NULL > I1.DA.ace
    cat I1.donor.r.txt | gawk -F '\t' '{ii=$1;m=$2;a1=$3;a2=$4;d1=$5;d2=$6;printf("Intron %s\n-D DA\nD %s__%d_r\nA %s__%d_r\nDonor %s\nAcceptor %s\n\n",ii,m,a1+1,m,a2-1,d1,d2);}' | grep -v NULL  >> I1.DA.ace
    echo "pparse I1.DA.ace" | ../../../bin/tace . -no_prompt
# check for same donor same acceptor
    cat I1.donor.f.txt | cut -f 2,3,4 | sort -k 1,1 -k 2,2n | gawk -F '\t' '{if($1==m && $2==a1){printf("Intron %s__%d_%d\nSame_donor %s__%d_%d\n\n",m,a1,a2,m,$2,$3);}m=$1;a1=$2;a2=$3;}' > I1.sameDA.ace  
    cat I1.donor.r.txt | cut -f 2,3,4 | sort -k 1,1 -k 3,3n | gawk -F '\t' '{if($1==m && $3==a2){printf("Intron %s__%d_%d\nSame_acceptor %s__%d_%d\n\n",m,a1,a2,m,$2,$3);}m=$1;a1=$2;a2=$3;}' >> I1.sameDA.ace  
    echo "pparse I1.sameDA.ace" | ../../../bin/tace . -no_prompt
    touch I1.introns.$MAGIC.DA.done
  popd
endif

echo   I1.DA.done

# associate donor acceptor to from_gene, meaning known in AceView
if (! -e tmp/INTRON_DB/$chrom/I1.$MAGIC.introns.DA2G2.done) then
  pushd tmp/INTRON_DB/$chrom
  ../../../bin/tace . <<EOF
    query find mrna COUNT locuslink == 1
    select -o I1.DAGG2.txt2 m,ll,g,ii,d,a from m in @,ll in m->locuslink,g in m->gene,ii in m->intron,d in ii->D, a in ii->A
    query find mrna COUNT gene==1 && !locuslink 
    select -o I1.DAGG2.txt1 m,ll,g,ii,d,a from m in @,ll=0,g in m->gene,ii in m->intron,d in ii->D, a in ii->A
    quit
EOF
  cat I1.DAGG2.txt1 | gawk -F '\t' '{ll=$2;gg=$3;ii=$4;d=$5;a=$6;printf("Intron %s\nGene %s\n\nDonor %s\nGene %s\n\nAcceptor %s\nGene %s\n\n",ii,gg,d,gg,a,gg);}' > I1.DAGG2.ace
  cat I1.DAGG2.txt2 | gawk -F '\t' '{ll=$2;g=$3;ii=$4;d=$5;a=$6;gg=ll;if(g!=ll)gg=ll"("g")";printf("Intron %s\nGene %s\n\nDonor %s\nGene %s\n\nAcceptor %s\nGene %s\n\n",ii,gg,d,gg,a,gg);}' >> I1.DAGG2.ace
  ../../../bin/tace . <<EOF
    pparse I1.DAGG2.ace
    query find intron ! gene
    kstore ii
    select -o I1.DAGG2.txt3c  ii,g1 from ii in @, d in ii->D, g1 in d->gene, a in ii->A, g2 in a->gene where g1 && g2 && g1 == g2
    kget ii
    select -o I1.DAGG2.txt3d  ii,g1 from ii in @, d in ii->D, g1 in d->gene, a in ii->A, g2 in a->gene where g1 && !g2
    kget ii
    select -o I1.DAGG2.txt3a  ii,g2 from ii in @, d in ii->D, g1 in d->gene, a in ii->A, g2 in a->gene where g2 && !g1
    save
    quit
EOF 
   cat I1.DAGG2.txt3[cda] | gawk -F '\t' '{ii=$1;g=$2;printf("Intron %s\nGene %s\n\n",ii,g);}' > I1.DAGG2.x.ace
  ../../../bin/tace . <<EOF
    pparse I1.DAGG2.x.ace
    query find intron !gene
    select -o I1.DAGG2.txt4 ii,g1,g2 from ii in @, d in ii->D where COUNT d->gene==1, a in ii->A where COUNT a->gene==1, g1 in d->gene,g2 in a->gene
    save
    quit
EOF
cat I1.DAGG2.txt4 | gawk -F '\t' '{ii=$1;split($2,g1,"(");split($3,g2,"(");printf("Intron %s\nGene %s__%s\nFusion\n\n",ii,g1[1],g2[1]);}' > I1.DAGG2.4.ace
  ../../../bin/tace . <<EOF
    read-models
    pparse I1.DAGG2.4.ace
    save
    quit
EOF

touch I1.$MAGIC.introns.DA2G2.done
  popd
endif

echo   I1.DA2G2.done

# check donor acceptors
if (! -e tmp/INTRON_DB/$chrom/I1.introns.$MAGIC.DA2G.done) then
  pushd tmp/INTRON_DB/$chrom
  foreach pass (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)
    ../../../bin/tace . <<EOF
      query find intron ! gene && same_donor
      bql -o I1.d2g.txt select ii,g from ii in @, i2 in ii->same_donor, g in i2->gene where g
      query find intron ! gene && same_acceptor
      bql -o I1.a2g.txt select ii,g from ii in @, i2 in ii->same_acceptor, g in i2->gene where g
      quit
EOF
    cat I1.d2g.txt I1.a2g.txt | gawk -F '\t' '{ii=$1;g=$2; printf ("Intron %s\nGene %s\n\n", ii,g);}' > I1.DA2G.$pass.ace
    tags I1.DA2G.$pass.ace
    set n=`wc I1.DA2G.$pass.ace | gawk '{print $1;}'`
    if ($n == 0) break
    echo "pparse I1.DA2G.$pass.ace" | ../../../bin/tace . -no_prompt
  end
  touch I1.$MAGIC.introns.DA2G.done
  popd
endif
if (! -e tmp/INTRON_DB/$chrom/I1.$MAGIC.introns.DDA2G.done) then
  pushd tmp/INTRON_DB/$chrom
    ../../../bin/tace . <<EOF
      query find intron from_gene
      bql -o I1.dd2g.txt select ii,g,d,a from ii in @, g in ii->from_gene, d in ii->D, a in ii->A
      query find intron from_gene && ! gene
      select -o I1.i2m2g.txt i,g from i in @, m in i->in_mrna, g in m->from_gene
      quit
EOF
    cat I1.i2m2g.txt | gawk -F '\t' '{i=$1;g=$2; printf ("Intron %s\nGene %s\n\n",i,g);}' > I1.DDA2G.ace
    cat I1.dd2g.txt  | gawk -F '\t' '{g=$2;d=$3;a=$4; printf ("Donor %s\nFrom_gene %s\n\nAcceptor %s\nFrom_gene %s\n\n", d,ii,a,g);}' >> I1.DDA2G.ace
    echo "pparse I1.DDA2G.ace" | ../../../bin/tace . -no_prompt
    ../../../bin/tace . <<EOF
      read-models
      query find intron known_donor 
      edit -D known_donor
      query find intron known_acceptor
      edit -D known_acceptor
      query find intron ; in_mrna
      edit known_donor
      spush
      edit known_acceptor
      query find intron known_donor ; >D ; >intron ; ! known_donor
      edit Known_donor
      query find intron known_acceptor ; >A ;  >intron ; ! known_acceptor
      edit Known_acceptor
      query find intron !gene && de_uno && (D || A)
      save
      quit
EOF
    touch I1.$MAGIC.introns.DDA2G.done
  popd
endif

echo   I1.DDA2G.done

if (! -e tmp/INTRON_DB/$chrom/I1.$MAGIC.introns.DA2Gb.done) then
  pushd tmp/INTRON_DB/$chrom

if (0) then
   cat I1.d2ii.txt | gawk -F '\t' '{n=split($2,aa,";");if (n>1){for(i=1;i<n;i++)for(j=i+1;j<=n;j++)printf("Intron %s\nSame_donor %s\n\n",aa[i],aa[j]);}}' > I1.DA2ii.ace
   cat I1.a2ii.txt | gawk -F '\t' '{n=split($2,aa,";");if (n>1){for(i=1;i<n;i++)for(j=i+1;j<=n;j++)printf("Intron %s\nSame_acceptor %s\n\n",aa[i],aa[j]);}}' >> I1.DA2ii.ace
   echo "pparse I1.DA2ii.ace" | ../../../bin/tace . -no_prompt
endif

  foreach pass (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)
    ../../../bin/tace . <<EOF
      query find intron ! gene && same_donor
      bql -o I1.d2g.txt select ii,g from ii in @, i2 in ii->same_donor, g in i2->gene where g
      query find intron ! gene && same_acceptor
      bql -o I1.a2g.txt select ii,g from ii in @, i2 in ii->same_acceptor, g in i2->gene where g
      quit
EOF

    cat I1.d2g.txt I1.a2g.txt | gawk -F '\t' '{ii=$1;g=$2; printf ("Intron %s\nGene %s\n\n", ii,g);}' > I1.DA2Gb.$pass.ace
    tags I1.DA2Gb.$pass.ace
    set n=`wc I1.DA2Gb.$pass.ace | gawk '{print $1;}'`
    if ($n == 0) break
    echo "pparse I1.DA2Gb.$pass.ace" | ../../../bin/tace . -no_prompt
  end
  touch I1.$MAGIC.introns.DA2Gb.done
  touch I1.$MAGIC.introns.donorAcceptor.done
  popd



  pushd tmp/INTRON_DB/$chrom
    ../../../bin/tace . <<EOF
      query find intron ! type
      bql -o I1.notype.txt select ii,d,a  from ii in @, d in ii->donor, a in ii->acceptor
      quit
EOF
cat I1.notype.txt | gawk -F '\t' '{if(length($2)==100 && length($3)==100) printf("Intron %s\nType %s_%s\n\n",$1, substr($2,51,2),substr($3,49,2));}' | gawk '/gt_ag/{print;next;}/ct_ac/{print;next;}/gc_ag/{print;next;}/^Type/{printf("Other %s\n",$2);next;}{print}' > I1.newtype.ace




    ../../../bin/tace . <<EOF
      pparse I1.newtype.ace
      find intron 
      bql -o I1.stats.txt select ii,type,rvy,av,n from ii in @, rvy in ii#rvy, av in ii#av, n in ii->rna_seq, type in ii->type
      save
      quit
EOF
cat I1.stats.txt | gawk '{ii=$1;split(ii,aa,"__");split(aa[2],bb,"_");a1=bb[1];a2=bb[2];if(a1>a2){a0=a1;a1=a2;a2=a0;}ip=ii;if(bb[1]>bb[2]){ip=aa[1] "__" bb[2] "_" bb[1];ii2[ip]=ii;}iip[ip]=1;}END{for(ii in iip){i2=ii2[ii];if(i2)printf("Intron %s\nHas_echo %s\nIs_echo %s\n\n",ii,i2,i2);}}' > I1.echo.ace
    ../../../bin/tace . <<EOF
      read-models
      pparse I1.echo.ace
      query find intron Is_echo && (gt_ag || gc_ag)
      edit -D Is_echo
      query find intron Is_echo && ! type && ! intmap
      kill
      query Find intron type
      bql -o I1.stats.txt select ii,type,echo,rvy,av,n from ii in @, rvy in ii#rvy, av in ii#av, n in ii->rna_seq, type in ii->type, echo in ii#is_echo
      save
      quit
EOF


  popd
endif



goto phaseLoop

########################################################################
## 
# check we have all intron length feet intmap
diffAcceptor:
echo -n "I1 phase $phase start :"
date

if (! -e tmp/INTRON_DB/$chrom/I1.$MAGIC.introns.donorAcceptor.done) goto phaseLoop
if (! -e tmp/INTRON_DB/$chrom/I1.$MAGIC.diffAcceptor.done2) then

  pushd tmp/INTRON_DB/$chrom
    ../../../bin/altintrons -db . --project $MAGIC -o I1.$MAGIC
    touch I1.$MAGIC.diffAcceptor.done
  popd
endif

goto phaseLoop

###############################################################
## intronComfirmation probably obsolete

intronComfirmation:
echo -n "I1 phase $phase start :"
date

if (-e  RESULTS/Expression/AceFiles/$MAGIC.introns.INTRON.u.ace.gz2) then
  set toto=GeneIndexDB/$MAGIC.intron_confirmation.ace
  echo $toto
  if (! -e $toto) then
    if (-e  RESULTS/Expression/AceFiles/$MAGIC.introns.INTRON.u.ace.gz2) then
      zcat RESULTS/Expression/AceFiles/$MAGIC.introns.INTRON.u.ace.gz | gawk '/^Intron/{z=$0}/_SumOfAllReadsInProject/{if ($4 >0) printf("%s\nRun_U %s %s %s %s %s %s\n\n",z,magic,$3,$4,$5,$6,$7);}' magic=$MAGIC > $toto
     echo "pparse $toto" | bin/tacembly GeneIndexDB -no_prompt
    else
      echo "Missing file  RESULTS/Expression/AceFiles/$MAGIC.introns.INTRON.u.ace.gz"
    endif
  endif
  
goto phaseLoop


###############################################################
## donorAcceptorZZ probably obsolete

donorAcceptorZZ:
goto phaseLoop

# donor acceptor
  bin/tacembly GeneIndexDB << EOF
    find intron
    show -a -f tmp/introns/I1.$MAGIC.any_intron.ace
    quit
EOF

  cat tmp/introns/I1.$MAGIC.any_intron.ace | gawk '/^Intron/{ii=$2;gsub(/\"/,"",$2);split($2,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];a2=bb[2];if(a1<a2){s="f";ds=1;}else {s="r";ds=-1;}printf("\nDonor %s__%d_%s\n",chrom, a1,s);printf("IntMap %s %d %d\n",chrom, a1-ds,a1);printf("Intron %s\n",ii);next;}/^Gene/{printf("%s %s\n",$1,$2);next;}/^From/{printf("%s %s\n",$1,$2);next;}/^In_mRNA/{printf("%s %s\n",$1,$2);next;}/^gt_ag/{printf("gt\n");next;}/^gc_ag/{printf("gc\n");next;}/^ct_ac/{printf("ct\n");next;}/^at_ac/{printf("at\n");next;}/^Other/{printf("Other\n");next;}/^AV/{print;next;}/^NM/{print;next;}/^XM/{print;next;}/^pg/{print;next;}/^RefSeq/{print $1;next;}/de_uno/{printf("Intron %s %s %d\n",ii,$2,$3);next;}' > tmp/introns/I1.$MAGIC.any_donor.ace

  cat tmp/introns/I1.$MAGIC.any_intron.ace | gawk '/^Intron/{ii=$2;gsub(/\"/,"",$2);split($2,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];a2=bb[2];if(a1<a2){s="f";ds=1;}else {s="r";ds=-1;}printf("\nAcceptor %s__%d_%s\n",chrom, a2,s);printf("IntMap %s %d %d\n",chrom, a2,a2+ds);printf("Intron %s\n",ii);next;}/^Gene/{printf("%s %s\n",$1,$2);next;}/^From/{printf("%s %s\n",$1,$2);next;}/^In_mRNA/{printf("%s %s\n",$1,$2);next;}/^gt_ag/{printf("ag\n");next;}/^gc_ag/{printf("ag\n");next;}/^ct_ac/{printf("ac\n");next;}/^at_ac/{printf("ac\n");next;}/^Other/{printf("Other\n");next;}/^AV/{print;next;}/^NM/{print;next;}/^XM/{print;next;}/^pg/{print;next;}/^RefSeq/{print $1;next;}/de_uno/{printf("Intron %s %s %d\n",ii,$2,$3);next;}' > tmp/introns/I1.$MAGIC.any_acceptor.ace

  bin/tacembly GeneIndexDB << EOF
    pparse tmp/introns/I1.$MAGIC.any_donor.ace
    pparse tmp/introns/I1.$MAGIC.any_acceptor.ace
    save
    bql -o tmp/introns/I1.$MAGIC.donors select d,t from d in ?Donor, t in d#gt
    bql -o tmp/introns/I1.$MAGIC.acceptors select a,t from a in ?Acceptor, t in a#ag
    quit
EOF
  cat tmp/introns/I1.$MAGIC.donors | gawk -F '\t' '/_f/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="gt"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1<=aOld+3)printf("Donor %s\nNext_gt %s %d\n\n", $1,dd, aOld-a1);}}' >  tmp/introns/I1.$MAGIC.donors.next_gt.ace
  cat tmp/introns/I1.$MAGIC.donors | gawk -F '\t' '/_r/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="gt"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1<=aOld+3)printf("Donor %s\nNext_gt %s %d\n\n", $1,dd, -aOld+a1);}}' >>  tmp/introns/I1.$MAGIC.donors.next_gt.ace
  cat tmp/introns/I1.$MAGIC.donors | sort -V -r | gawk -F '\t' '/_f/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="gt"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1>=aOld-3)printf("Donor %s\nNext_gt %s %d\n\n", $1,dd, -aOld+a1);}}' >>  tmp/introns/I1.$MAGIC.donors.next_gt.ace
  cat tmp/introns/I1.$MAGIC.donors | sort -V -r | gawk -F '\t' '/_r/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="gt"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1>=aOld-3)printf("Donor %s\nNext_gt %s %d\n\n", $1,dd, aOld-a1);}}' >>  tmp/introns/I1.$MAGIC.donors.next_gt.ace

  cat tmp/introns/I1.$MAGIC.acceptors | gawk -F '\t' '/_r/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="ag"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1<=aOld+3)printf("Acceptor %s\nNext_ag %s %d\n\n", $1,dd, aOld-a1);}}' >  tmp/introns/I1.$MAGIC.acceptors.next_ag.ace
  cat tmp/introns/I1.$MAGIC.acceptors | gawk -F '\t' '/_f/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="ag"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1<=aOld+3)printf("Acceptor %s\nNext_ag %s %d\n\n", $1,dd, -aOld+a1);}}' >>  tmp/introns/I1.$MAGIC.acceptors.next_ag.ace
  cat tmp/introns/I1.$MAGIC.acceptors | sort -V -r | gawk -F '\t' '/_r/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="ag"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1>=aOld-3)printf("Acceptor %s\nNext_ag %s %d\n\n", $1,dd, -aOld+a1);}}' >>  tmp/introns/I1.$MAGIC.acceptors.next_ag.ace
  cat tmp/introns/I1.$MAGIC.acceptors | sort -V -r | gawk -F '\t' '/_f/{dd=$1;split(dd,aa,"__");chrom=aa[1];split(aa[2],bb, "_");a1=bb[1];if($2=="ag"){cc=chrom;aOld=a1;dd=$1;}else{if(cc==chrom && a1>=aOld-3)printf("Acceptor %s\nNext_ag %s %d\n\n", $1,dd, aOld-a1);}}' >>  tmp/introns/I1.$MAGIC.acceptors.next_ag.ace

  if (-e  tmp/introns/I1.$MAGIC.captured_genes.txt) \rm  tmp/introns/I1.$MAGIC.captured_genes.txt

  bin/tacembly GeneIndexDB << EOF
    read-models
    parse tmp/introns/I1.$MAGIC.donors.next_gt.ace
    parse tmp/introns/I1.$MAGIC.acceptors.next_ag.ace 
    save
    quit
EOF
  touch tmp/introns/I1.$MAGIC.donorAcceptor.done

goto phaseLoop

###############################################################
### CAPTURE

capture:
  if ($?CAPTURES && -e tmp/METADATA/$MAGIC.av.captured_genes.ace && ! -e tmp/INTRON_DB/$chrom/I1.$MAGIC.capture.done) then 
    pushd tmp/INTRON_DB/$chrom
      ../../../bin/tace . << EOF
        read-models
        find gene
        spush
        parse  ../../../tmp/METADATA/$MAGIC.av.captured_genes.ace
        find gene
        sxor
        spop
        kill
        save
        bql -o  I1.$MAGIC.captured_introns.txt select ii,c from g in ?Gene, c in g->capture where c, ii in g->Intron  
        quit
EOF

     cat I1.$MAGIC.captured_introns.txt  | gawk -F '\t' '{if($1 != old)printf ("\nIntron %s\n",$1);old=$1;printf("Capture %s\n",$2);}END{printf("\n");}' > I1.$MAGIC.captured_introns.ace
     ../../../bin/tace .  << EOF
       parse I1.$MAGIC.captured_introns.ace
       save
       quit
EOF

      touch tmp/INTRON_DB/$chrom/I1.$MAGIC.capture.done
    popd  
   endif
  touch tmp/INTRON_DB/$chrom/I1.$MAGIC.capture.done
goto phaseLoop


###################################################################
## is that useful or obsolete ?

  set capt=A1
  set chrom2=$chrom'__'
  echo "chrom2=$chrom2"
  pushd tmp/INTRON_DB/$chrom
    ../../../bin/tace . <<EOF
      query find intron Deep 
      edit -D Deep
      query find intron 
      spush
      query IS $chrom2*
      sminus
      spop
      kill
      query find intron de_uno 
      show -a -f I1.introns.final.preace
      query find intron de_uno && capture == $capt
      show -a -f I1.introns.final.$capt.preace
      save
      quit
EOF
    cat I1.introns.final.preace | gawk '/^$/{print}/^Intron/{print}/^de_uno/{print}' > I1.$MAGIC.de_uno.ace
    cat I1.introns.final.preace | gawk '/^de_uno/{next}{print}' > I1.$MAGIC.info.ace
    tags  I1.$MAGIC.de_uno.ace
    tags  I1.$MAGIC.info.ace
    cat I1.introns.final.$capt.preace | gawk '/^$/{print}/^Intron/{print}/^de_uno/{print}' > I1.$MAGIC.de_uno.$capt.ace
    cat I1.introns.final.$capt.preace | gawk '/^de_uno/{next}{print}' > I1.$MAGIC.info.$capt.ace
    tags  I1.$MAGIC.de_uno.$capt.ace
    tags  I1.$MAGIC.info.$capt.ace
  popd

  cat MetaDB/$MAGIC/RunsList  MetaDB/$MAGIC/GroupIntronList > tmp/INTRON_DB/$chrom/$MAGIC.RunList 

  touch tmp/INTRON_DB/$chrom/I1.$MAGIC.capture.done
goto phaseLoop

###############################################################
## phaseLoop

phaseLoop:
  echo -n "I1 phase $phase done :"
  date
end
  echo "phases $phases done"
  echo I1.intronDB.tcsh phase $phase  done
  exit 0


###############################################################
###############################################################
