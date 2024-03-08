#!bin/tcsh -f
set chrom=$1
setenv ici `pwd`

if (! -d tmp/X.$MAGIC/XH$chrom) mkdir  tmp/X.$MAGIC/XH$chrom
if (! -d tmp/X.$MAGIC/XH$chrom/database) then
   pushd tmp/X.$MAGIC/XH$chrom
      mkdir database
      ln -s ../../../metaData/tables
      if (-e TABIX) \rm TABIX
      ln -s ../../TABIX
      mkdir wspec
        \cp  $ici/metaData/wspec.aceview_web_site/*.wrm wspec
        pushd wspec
          \rm models.wrm
          ln -s $ici/metaData/wspec.aceview_web_site/models.wrm
        popd
        set mynam=`whoami`
        set n=`cat wspec/passwd.wrm | gawk '{if($1==mynam)n++}END{printf("%d", 0+n)}' mynam=$mynam`
     if ($n == 0) then
           whoami >> wspec/passwd.wrm
        endif
      echo y |  $ici/bin/tacembly .
   popd
endif

echo "Clone MainClone\nMain_clone\nMainTitle $MAGIC\n\n" > tmp/X.$MAGIC/XH$chrom/Main_clone.ace

set target=`echo $Etargets | gawk '{print $1}'`

if (! -e tmp/X.$MAGIC/XH$chrom/f3.genome.done2) then

    if ($target == EBI && $species == Dmelanogaster && -e  TARGET/GTF/$species.$target.gtf.gz && -d TARGET/GENES && ! -e TARGET/GENES/f3.gtf.$target.gene2FB_symbol.ace) then
       gunzip -c TARGET/GTF/$species.$target.gtf.gz | gawk -F '\t' '{if($3 == "gene"){z=$9;split(z,aa,";");split(aa[1],bb," ");split(aa[3],cc," ");if(bb[1]=="gene_id" && cc[1]=="gene_name")printf("Sequence X__%s\nFB_symbol %s\n\n",bb[2],cc[2]);}}' > TARGET/GENES/f3.gtf.$target.gene2FB_symbol.ace
    endif
    if ($target == FlyBase && -e  TARGET/GTF/$species.$target.gtf.gz && -d TARGET/GENES && ! -e TARGET/GENES/f3.gtf.$target.gene2FB_symbol.ace) then
       gunzip -c TARGET/GTF/$species.$target.gtf.gz | gawk -F '\t' '{if($3 == "gene"){z=$9;gsub(/\"/,"",z);split(z,aa,";");split(aa[1],bb," ");split(aa[2],cc," ");if(bb[1]=="gene_id" && cc[1]=="gene_symbol")printf("Sequence X__%s\nFB_symbol %s\n\n",bb[2],cc[2]);}}' > TARGET/GENES/f3.gtf.$target.gene2FB_symbol.ace
    endif

printf "-R Sequence c_$chrom $chrom\n\n" >  tmp/X.$MAGIC/XH$chrom/f3.rename_chrom.ace

    bin/tacembly tmp/X.$MAGIC/XH$chrom <<EOF
      read-models
      find est XH*
      kill
      pparse TARGET/CHROMS/$species.chrom_$chrom.fasta.gz
      pparse  tmp/f1.strategy.ace
      // pparse  tmp/pA/$chrom/$MAGIC.pA.$chrom.feature.ace 
      pparse  tmp/TABIX/$MAGIC/$chrom.tabix.ace
      pparse MetaDB/$MAGIC/runs.ace
      parse tmp/METADATA/gtf.$target.f.intron.ace
      parse tmp/METADATA/gtf.$target.r.intron.ace
      parse tmp/METADATA/gtf.RefSeq.transcripts.ace.gz
      parse tmp/METADATA/$MAGIC.$target.captured_genes.ace
      parse tmp/METADATA/av.GENE.info.ace
      query find sequence locuslink
      spush
      query intmap == $chrom
      sminus
      spop
      kill
      query find sequence $chrom
      kstore ss
      acem
        make_subseq -dna c t 1000000 10000 // this breaks my 6 contigs into tiles
        quit                    // oct 15 2001, i changed from 400 kb to 600 kb
      kget ss
      Follow DNA
      kill 
      undo
      edit -D DNA
      parse tmp/X.$MAGIC/XH$chrom/f3.rename_chrom.ace
      query find sequence genomic
      follow source
      edit  YBR_contig  // NT_sequences 
 
      save
      quit
EOF

# in Xcds the score is independant of the coverage, it just reflects the length of the CDS, so that all long CDS are shown with all their exons
set ffcds=tmp/METADATA/gtf.$target.ns.cds.spongeZZZZZ
if (! -e $ffcds) then
  set ffcds=tmp/METADATA/$target.ns.cds.sponge
endif
echo AAA $ffcds
if (-e $ffcds) then
  cat $ffcds ZZZZZ $ffcds | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if (chrom != $3)next;}{a1=$4;a2=$5;ln=a2-a1;if(ln<0)ln = -ln;ln++;}{if(zz<1){t2ln[$1]+=ln;next;}}{nnn=t2ln[$1];if(nnn >= 450)printf("Xcds_%s__%d_%d\t1\t%d\t%s\t%d\t%d\t%d\n",chrom,a1,a2,ln,chrom,a1,a2,nnn);}' chrom=$chrom | sort -u | sort -V >  tmp/X.$MAGIC/XH$chrom/f3.$target.cds.txt2
    bin/dna2dna -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow  tmp/X.$MAGIC/XH$chrom/f3.$target.cds.txt2 >  tmp/X.$MAGIC/XH$chrom/f3.$target.cds.fasta
    cat  tmp/X.$MAGIC/XH$chrom/f3.$target.cds.txt2 | gawk -F '\t' '{printf("Sequence %s\ncdna_clone %s\nForward\nColour LIGHTVIOLET\nIntMap %s %d %d\nIs_read\nComposite %d\n\n",$1,$1,$4,$5,$6,$7);}' >  tmp/X.$MAGIC/XH$chrom/f3.$target.cds.ace
    bin/tacembly  tmp/X.$MAGIC/XH$chrom << EOF
      pparse  tmp/X.$MAGIC/XH$chrom/f3.$target.cds.ace
      pparse  tmp/X.$MAGIC/XH$chrom/f3.$target.cds.fasta
      save
      quit
EOF
endif

    bin/gene2chrom2 -any -gs -i tmp/X.$MAGIC/XH$chrom  >! tmp/X.$MAGIC/XH$chrom/f3.g2c.gsi.ace
    bin/gene2chrom2 -any -pg -i tmp/X.$MAGIC/XH$chrom  >! tmp/X.$MAGIC/XH$chrom/f3.g2c.pgi.ace

    bin/tacembly  tmp/X.$MAGIC/XH$chrom << EOF
      pparse  tmp/X.$MAGIC/XH$chrom/f3.g2c.gsi.ace
      pparse  tmp/X.$MAGIC/XH$chrom/f3.g2c.pgi.ace
      save
      quit
EOF

    touch  tmp/X.$MAGIC/XH$chrom/f3.genome.done
endif

if (! -e   tmp/X.$MAGIC/XH$chrom/f3.cosmid2map.$chrom.txt) then
    echo "get cosmid2map $chrom"
    bin/tacembly tmp/X.$MAGIC/XH$chrom <<EOF > /dev/null
      query find sequence genomic  &&  NOT Is_gene_tile 
      bql -a -o tmp/X.$MAGIC/XH$chrom/f3.cosmid2map.$chrom.txt  select s,m,a1,a2 from s in @, m in s->intmap, a1 in m[1], a2 in a1[1]
      quit
EOF
endif

echo "XXXXXXX f3.parse"
if (! -e tmp/X.$MAGIC/XH$chrom/f3.parse.done) then

  touch tmp/X.$MAGIC/XH$chrom/f3.genes.ace
  if (1 && -e TARGET/GENES/av.gene.ace) then
    # cree des gene box
    cat TARGET/GENES/av.gene.ace | gawk '/^IntMap/{m=$2;gsub(/\"/,"",m);if(m!=chrom)next;{ok=1;gg = gg "IntMap " m " " $3 " "  $4 "\n"; cc=cc "Genes " g " " $3 " " $4 "\n";}next;}/^av/{next;}/^Transcribed_gene/{next;}/^Gene /{if(ok==1)print gg "\n";ok=0;gg=$0 "\n"; g = $2;next;}{gg = gg $0 "\n" ;next;}END{if(ok)print gg "\n";print "Sequence " chrom "\n" cc "\n"}' chrom=$chrom > tmp/X.$MAGIC/XH$chrom/f3.genes.ace
  else 
    # cree des gene models de type le premier de Etargets (may be magic, on preferererait av ou RefSeq)
    set target=`echo $Etargets | gawk '{print $1}'`
    set target=av
    if (1 && -e tmp/METADATA/gtf.$target.transcripts.ace.gz) then
      gunzip -c tmp/METADATA/gtf.$target.transcripts.ace.gz > tmp/X.$MAGIC/XH$chrom/f3.genes.ace
    endif
  endif

  set ch=""
  if ($species == worm) set ch=CHROMOSOME_

  pushd tmp/X.$MAGIC/EHITS/$chrom
    $ici/scripts/rrf
    $ici/bin/tacembly  $ici/tmp/X.$MAGIC/XH$chrom < _r
  popd

  bin/tacembly  tmp/X.$MAGIC/XH$chrom << EOF
    pparse  tmp/X.$MAGIC/XH$chrom/f3.genes.ace
    parse metaData/methods.ace
    save
    query find intron
    list -a -f tmp/X.$MAGIC/EHITS/$chrom/f3.intron.list
    quit
EOF

if (! -e tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.ace && -d IntronDB) then
  bin/tacembly IntronDB << EOF
    key tmp/X.$MAGIC/EHITS/$chrom/f3.intron.list
    show -a -f tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.preace Magic_any_any
    query find intron IS $chrom* && In_mrna
    show -a -f tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.other.preace
    query find intron IS $chrom* AND Magic_any_any >= 5
    show -a -f tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.any.preace Magic_any_any
    quit
EOF
  ls -ls  tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.any.preace

  if (-e  tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.preace) then
    cat  tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.preace | gawk '/^$/{print}{gsub(/\"/,"",$0);}/^Intron/{printf("Intron %s\n",$2);next;}/^Magic_any_any/{printf("Group_U Magic_any 0 %d\n",$2);}'  > tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.ace
  endif
endif

  if (-e  tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.any.preace && ! -e tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.other.ace) then

    cat  tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.any.preace | gawk '{gsub(/\"/,"",$0);}/^Intron/{if(mx>=5){if(a1==a2)a2=a1+30*s;if(d1==d2)d1=d2-30*s;}if(d1!=d2 && a1 != a2)printf("\tExonExon %d %d %d %d",d1,d2,a1,a2);split($2,aa,"__");split(aa[2],bb,"_");s=1;if(bb[1]>bb[2])s=-1;printf("\n%s %d %d",aa[1],bb[1],bb[2]);d1=d2=bb[1]-s;a1=a2=bb[2]+s;mx=0;next;}/^Donor /{printf("\t%s %s",$1,$2);}/^Acceptor /{printf("\t%s %s",$1,$2);}/^Magic_any_any/{if($2>=5){mx=$2;printf("\tGroup_U Magic_any 0 %d",$2);next;}}/^In_mRNA/{if ($5=="Acceptor_exon"){x=$6;y=$7;if(x != a1)next;da=s*(y-x);if(da<6)next;if(da>30)da=30;if(a2==a1 || s*(a2-a1)>da)a2=a1+s*da;next;}if ($5=="Donor_exon"){x=$6;y=$7;if(y != d2)next;da=s*(y-x);if(da<6)next;if(da>30)da=30;if(d2==d1 || s*(d2-d1)>da)d1=d1-s*da;next;}next;}END{if(d1!=d2 && a1 != a2)printf("\tExonExon %d %d %d %d",d1,d2,a1,a2);printf("\n");}' | grep  Magic_any | grep ExonExon > tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.other.txts
  
 cat  tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.other.txts |   gawk '/ExonExon/{c=$1;ii1=$2;ii2=$3;s=1;if(ii1>ii2)s=-1;i=index($0,"ExonExon");split(substr($0,i),aa," ");d1=aa[2];d2=aa[3];a1=aa[4];a2=aa[5];dd=s*(d2-d1)+1;da=s*(a2-a1)+1;printf("XK_%s__%d_%d\t1\t%d\t%s\t%d\t%d\n",c,ii1,ii2,dd,c,d1,d2);printf("XK_%s__%d_%d\t%d\t%d\t%s\t%d\t%d\n",c,ii1,ii2,dd+1,dd+da,c,a1,a2);next;}' > tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.other.shadow
  
  bin/dna2dna -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow  tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.other.shadow >  tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.other.fasta

  cat  tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.other.txts |   gawk '/ExonExon/{c=$1;ii1=$2;ii2=$3;s=1;if(ii1>ii2)s=-1;i=index($0,"Magic_any");split(substr($0,i),aa," ");mx=aa[3];i=index($0,"ExonExon");split(substr($0,i),aa," ");d1=aa[2];d2=aa[3];a1=aa[4];a2=aa[5];dd=s*(d2-d1)+1;da=s*(a2-a1)+1;printf("Sequence XK_%s__%d_%d\nIs_read\ncDNA_clone XK_%s__%d_%d\nColour PALEORANGE\nForward\nComposite %s\nIntron %s__%d_%d\nIntMap %s %d %d\n\n",c,ii1,ii2,c,ii1,ii2,mx,c,ii1,ii2,c,d1,a2);next;}' > tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.other.ace

  endif

if (-e  RESULTS/Expression/AceFiles/$MAGIC.introns.INTRON.u.ace.gz ) then
  gunzip -c RESULTS/Expression/AceFiles/$MAGIC.introns.INTRON.u.ace.gz | gawk '/^Intron/{ok=0;split($2,aa,"__");gsub(/\"/,"",aa[1]);if(aa[1]==chrom)ok=1;if(ok==1)print;next;}/Group_U  _SumOfAllReadsInProject/{if (ok==1)printf("RNA_seq %d\n\n",$6);}' chrom=$chrom > tmp/X.$MAGIC/EHITS/$chrom/f3.transcriptsIntronSupport.ace
endif

if (-e tmp/X.$MAGIC/EHITS/f2.allDoubleIntrons.txt.gz) then

    gunzip -c tmp/X.$MAGIC/EHITS/f2.allDoubleIntrons.txt.gz | gawk -F '\t' '{if($3<1)next;split($1,dd,"___");split(dd[1],aa,"__");c=aa[1];if(c != chrom)next;split(aa[2],bb,"_");a1=bb[1];a2=bb[2];split(dd[2],bb,"_");b1=bb[1];b2=bb[2];s=1;if(a1>a2)s=-1;printf("Sequence XW_%s__%s_%s_%s_%s\nColour ORANGE\nForward\nComposite %d\nIntron %s__%d_%d\nIntron %s__%d_%d\nIs_read\ncDNA_clone XW__%s_%s_%s_%s_%s\nIntMap %s %d %d\n\n",c,a1,a2,b1,b2,$3,c,a1,a2,c,b1,b2,c,a1,a2,b1,b2,c,a1-30*s,b2+30*s);}' chrom=$chrom > tmp/X.$MAGIC/XH$chrom/f3.allDoubleIntrons.ace
     gunzip -c tmp/X.$MAGIC/EHITS/f2.allDoubleIntrons.txt.gz | gawk -F '\t' '{if($3<1)next;split($1,dd,"___");split(dd[1],aa,"__");c=aa[1];if(c != chrom)next;split(aa[2],bb,"_");a1=bb[1];a2=bb[2];split(dd[2],bb,"_");b1=bb[1];b2=bb[2];s=1;if(a1>a2)s=-1;if(a1>0 && b2 + 30*s >0 && b1>0 && b2 > 0 && s*(b1-a2)+1+30 > 0 && s*(b1-a2)+60 > 0){printf("XW_%s__%s_%s_%s_%s\t1\t30\t%s\t%d\t%d\n",c,a1,a2,b1,b2,c,a1-30*s,a1-s);printf("XW_%s__%s_%s_%s_%s\t%d\t%d\t%s\t%d\t%d\n",c,a1,a2,b1,b2,31,29+s*(b1-a2),c,a2+s,b1-s);printf("XW_%s__%s_%s_%s_%s\t%d\t%d\t%s\t%d\t%d\n",c,a1,a2,b1,b2,s*(b1-a2)+1+30,s*(b1-a2)+60,c,b2+s,b2+30*s);}}' chrom=$chrom  > tmp/X.$MAGIC/XH$chrom/f3.allDoubleIntrons.shadow
 
    bin/dna2dna -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow  tmp/X.$MAGIC/XH$chrom/f3.allDoubleIntrons.shadow >  tmp/X.$MAGIC/XH$chrom/f3.allDoubleIntrons.fasta

 endif

if (0) then
# hack to debug the faulty XY 2016_03_21
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    gunzip -c tmp/PHITS_genome/$run/*.doubleintrons.gz | grep DOUBLEINTRON | grep 22408 | grep 22429 | grep 22511 | gawk '{printf("%s\t",run);print;}' run=$run >> tatou76
  end
  foreach ff (`ls tmp/PHITS_genome/SRR1509510/*doubleintrons.gz `)
    gunzip -c $ff  | grep DOUBLEINTRON |  grep 22408 | grep 22429 | grep 22511 | gawk '{printf("%s\t",run);print;}' run=$ff
  end
endif

# DOUBLEINTRON de_uno found on the fly by the aligner
echo AAAAA
if (! -e tmp/introns/$MAGIC.allDoubleIntronsGenomic.$chrom.txt) then
  # collate the data 
  if (-e   tmp/introns/$MAGIC.allDoubleIntronsGenomic.$chrom.txt1) \rm  tmp/introns/$MAGIC.allDoubleIntronsGenomic.$chrom.txt1
  foreach run (`cat MetaDB/$MAGIC/RunList`)
    gunzip -c tmp/PHITS_genome/$run/*.doubleintrons.gz | grep DOUBLEINTRON | sort -k 2,2 -k 4,4n -k 5,5n -k 6,6n -k 7,7n | gawk -F '\t' '/^DOUBLEINTRON/{c=$2;if(c != chrom)next;a=$2 "\t" $4 "\t" $5 "\t" $6 "\t" $7; nn[a]+=$9;x=$3;if ($5<$6 && (x1[a]+0==0 || x1[a]>x)){if(x<$4-30)x=$4-30;x1[a]=x;}if ($5>$6 && (x1[a]+0==0 || x1[a]<x)){if(x>$1+30)x=$1+30;x1[a]=x;}x=$8;if ($5<$6 && (x2[a]+0==0 || x2[a]<x)){if(x>$7+30)x=$7+30;x2[a]=x;}if ($5>$6 && (x2[a]+0==0 || x2[a]<x)){if(x<$7-30)x=$7-30;x2[a]=x;}}END{for(a in nn){split(a,aa,"\t");printf("DOUBLEINTRON\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",aa[1],x1[a],aa[2],aa[3],aa[4],aa[5],x2[a],nn[a]);}}' chrom=$chrom |  sort -k 2,2 -k 4,4n -k 5,5n -k 6,6n -k 7,7n  >>  tmp/introns/$MAGIC.allDoubleIntronsGenomic.$chrom.txt1
  end
  cat tmp/introns/$MAGIC.allDoubleIntronsGenomic.$chrom.txt1 |  grep DOUBLEINTRON | sort -k 2,2 -k 4,4n -k 5,5n -k 6,6n -k 7,7n | gawk -F '\t' '/^DOUBLEINTRON/{c=$2;if(c != chrom)next;a=$2 "\t" $4 "\t" $5 "\t" $6 "\t" $7; nn[a]+=$9;x=$3;if ($5<$6 && (x1[a]+0==0 || x1[a]>x))x1[a]=x;if ($5>$6 && (x1[a]+0==0 || x1[a]<x))x1[a]=x;x=$8;if ($5<$6 && (x2[a]+0==0 || x2[a]<x))x2[a]=x;if ($5>$6 && (x2[a]+0==0 || x2[a]<x))x2[a]=x;}END{for(a in nn){split(a,aa,"\t");printf("DOUBLEINTRON\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",aa[1],x1[a],aa[2],aa[3],aa[4],aa[5],x2[a],nn[a]);}}' chrom=$chrom | sort -k 2,2 -k 4,4n -k 5,5n -k 6,6n -k 7,7n  >  tmp/introns/$MAGIC.allDoubleIntronsGenomic.$chrom.txt
  \rm  tmp/introns/$MAGIC.allDoubleIntronsGenomic.$chrom.txt1
endif 
  # export the ace file


if (! -e tmp/X.$MAGIC/XH$chrom/f3.allDoubleIntronsGenomic.ace) then
  cat   tmp/introns/$MAGIC.allDoubleIntronsGenomic.$chrom.txt |  gawk  -F '\t' '/^DOUBLEINTRON/{if($9<2)next;s=1;if($4>$5)next;s=1;dx1=$4-$3+1;dx2=$6-$5+1;dx3=$8-$7+1;di1=$5-$4-1;di2=$7-$6-1;if(di1<0 || di2<0 || dx1<0 || dx2<0 || dx3<0)next;printf ("Sequence XY_%s__%d_%d_%d_%d\n", $2,$4+s,$5-s,$6+s,$7-s); printf ("cDNA_clone XY_%s__%d_%d_%d_%d\n", $2,$4+s,$5-s,$6+s,$7-s); printf("COLOUR CERISE\nForward\nComposite %d\nIs_read\n", $9);printf("IntMap %s %d %d\n",$2,$3,$8);printf("Intron %s__%s_%s\nIntron %s__%d_%d\n\n",$2,$4+s,$5-s,$2,$6+s,$7-s);}' >   tmp/X.$MAGIC/XH$chrom/f3.allDoubleIntronsGenomic.ace
  cat  tmp/introns/$MAGIC.allDoubleIntronsGenomic.$chrom.txt |  gawk  -F '\t' '/^DOUBLEINTRON/{if($9<2)next;s=1;if($4<$5)next;s=-1;dx3=$7-$8+1;dx2=$5-$6+1;dx1=$3-$4+1;di2=$6-$7-1;di1=$4-$5-1;if(di1<0 || di2<0 || dx1<0 || dx2<0 || dx3<0)next;printf ("Sequence XY_%s__%d_%d_%d_%d\n", $2,$4+s,$5-s,$6+s,$7-s); printf ("cDNA_clone XY_%s__%d_%d_%d_%d\n",$2,$4+s,$5-s,$6+s,$7-s); printf("COLOUR CERISE\nForward\nComposite %d\nIs_read\n", $9);printf("IntMap %s %d %d\n",$2,$3,$8);printf("Intron %s__%s_%s\nIntron %s__%d_%d\n\n",$2,$4+s,$5-s,$2,$6+s,$7-s);}' >>   tmp/X.$MAGIC/XH$chrom/f3.allDoubleIntronsGenomic.ace

  # create the shadows

  cat  tmp/introns/$MAGIC.allDoubleIntronsGenomic.$chrom.txt |  gawk  -F '\t' '/^DOUBLEINTRON/{if($9<2)next;c=$2;if($4>$5)next;nam="XY_" c "__" $4+1 "_" $5-1 "_" $6+1 "_" $7-1; dx1=$4-$3+1;dx2=$6-$5+1;dx3=$8-$7+1;di1=$5-$4-1;di2=$7-$6-1;if(di1<0 || di2<0 || dx1<0 || dx2<0 || dx3<0)next;printf("%s\t%d\t%d\t%s\t%d\t%d\n",nam,1,dx1,c,$3,$4);printf("%s\t%d\t%d\t%s\t%d\t%d\n",nam,dx1+1,dx1+dx2,c,$5,$6);printf("%s\t%d\t%d\t%s\t%d\t%d\n",nam,dx1+dx2+1,dx1+dx2+dx3,c,$7,$8);}' >  tmp/X.$MAGIC/XH$chrom/f3.allDoubleIntronsGenomic.shadow
  cat  tmp/introns/$MAGIC.allDoubleIntronsGenomic.$chrom.txt |  gawk  -F '\t' '/^DOUBLEINTRON/{if($9<2)next;c=$2;if($4<$5)next;nam="XY_" c "__" $4-1 "_" $5+1 "_" $6-1 "_" $7+1; ;s=-1;dx3=$7-$8+1;dx2=$5-$6+1;dx1=$3-$4+1;di2=$6-$7-1;di1=$4-$5-1;if(di1<0 || di2<0 || dx1<0 || dx2<0 || dx3<0)next;printf("%s\t%d\t%d\t%s\t%d\t%d\n",nam,1,dx1,c,$3,$4);printf("%s\t%d\t%d\t%s\t%d\t%d\n",nam,dx1+1,dx1+dx2,c,$5,$6);printf("%s\t%d\t%d\t%s\t%d\t%d\n",nam,dx1+dx2+1,dx1+dx2+dx3,c,$7,$8);}' >>  tmp/X.$MAGIC/XH$chrom/f3.allDoubleIntronsGenomic.shadow

  # export the fasta
  dna2dna -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow  tmp/X.$MAGIC/XH$chrom/f3.allDoubleIntronsGenomic.shadow > tmp/X.$MAGIC/XH$chrom/f3.allDoubleIntronsGenomic.fasta
  
endif
echo BBBB

set ggNS="toto"
if (-e MetaDB/$MAGIC/GroupW_new_exonList) then
  set ggNS=` cat  MetaDB/$MAGIC/GroupW_new_exonList | gawk '{printf("%s",$1);exit;}'`
endif


if (-d GeneIndexDB) then
  bin/tacembly GeneIndexDB << EOF
    key tmp/X.$MAGIC/EHITS/$chrom/f3.intron.list
    show -a -f tmp/X.$MAGIC/EHITS/$chrom/f3.introns.preace
    quit
EOF

  if (-e tmp/X.$MAGIC/EHITS/$chrom/f3.introns.preace) then
    cat  tmp/X.$MAGIC/EHITS/$chrom/f3.introns.preace | gawk '/^$/{print}{gsub(/\"/,"",$0);}/^Intron/{printf("Intron %s\n",$2);next;}/^Group/{if($2 == nx)print;}/Validated_u/{if($2 == "any")print}/de_[du][nu]o/{if($2 == "any")print}' nx=$ggNS  > tmp/X.$MAGIC/EHITS/$chrom/f3.introns.ace
    \rm   tmp/X.$MAGIC/EHITS/$chrom/f3.introns.preace 
  endif
endif

hello:
set ggs=toto
set ggNS="toto"
if (-e MetaDB/$MAGIC/GroupW_new_exonList) then
  set ggNS1=` cat  MetaDB/$MAGIC/GroupW_new_exonList | gawk '{printf("%s",$1);exit;}'`  
  if ($ggNS1 != "") set ggNS=$ggNS1
endif
set ggS="toto"
ls -ls  MetaDB/$MAGIC/GroupW_strandedList
if (-e MetaDB/$MAGIC/GroupW_strandedList) then
  set ggS1=` cat  MetaDB/$MAGIC/GroupW_strandedList | gawk '{printf("%s",$1);exit;}'`
  if ($ggS1 != "") set ggS=$ggS1
endif

set WGR="toto"
# restranding the non stranded wiggle using the stranded wiggle
  if (-d  tmp/WIGGLEGROUP/$ggS) set WGR=WIGGLEGROUP
  if (-d  tmp/WIGGLERUN/$ggS) set WGR=WIGGLERUN
echo "ggS=$ggs#"
echo "WGR=$WGR#"
if ($ggS != $ggNS && $ggS != toto && $ggNS != toto && $WGR != toto && ! -e  tmp/$WGR/$ggNS/$chrom/R.chrom.u0.r.BF.gz) then

  echo "restranding the non stranded wiggle using the stranded wiggle"
  mv  tmp/$WGR/$ggNS/$chrom/R.chrom.u.f.BF.gz  tmp/$WGR/$ggNS/$chrom/R.chrom.u0.f.BF.gz
  mv  tmp/$WGR/$ggNS/$chrom/R.chrom.u.r.BF.gz  tmp/$WGR/$ggNS/$chrom/R.chrom.u0.r.BF.gz

  scripts/wiggle2tabix.tcsh $ggs $chrom u $chrom

  bin/wiggle -wiggleRatioDamper 10 -I BF -O BF  -wiggle1 tmp/$WGR/$ggNS/$chrom/R.chrom.u0.f.BF.gz  -wiggle2 tmp/$WGR/$ggNS/$chrom/R.chrom.u0.r.BF.gz  -swiggle1 tmp/$WGR/$ggS/$chrom/R.chrom.u.f.BF.gz  -swiggle2 tmp/$WGR/$ggS/$chrom/R.chrom.u.r.BF.gz -gzo -o    tmp/$WGR/$ggNS/$chrom/R.chrom.u.f
  bin/wiggle -wiggleRatioDamper 10 -I BF -O BF  -wiggle1 tmp/$WGR/$ggNS/$chrom/R.chrom.u0.r.BF.gz  -wiggle2 tmp/$WGR/$ggNS/$chrom/R.chrom.u0.f.BF.gz  -swiggle1 tmp/$WGR/$ggS/$chrom/R.chrom.u.r.BF.gz  -swiggle2 tmp/$WGR/$ggS/$chrom/R.chrom.u.f.BF.gz -gzo -o    tmp/$WGR/$ggNS/$chrom/R.chrom.u.r

endif

###### grab the coverons and the gene-ends using multiPeaks (XG stranded, XH non stranded)
\rm  tmp/X.$MAGIC/XH$chrom/f3.Xends.*
set WGR=toto
echo "grab the XG/XH"
foreach XGH (XG XH)
  if ($XGH == XG && $ggS == toto) continue 
  if ($XGH == XH && $ggNS == toto) continue 
  set ggs=$ggS
  if ($XGH == XH) set ggs=$ggNS
  if ($XGH == XH && $ggNS == $ggS) continue

  if (-d  tmp/WIGGLEGROUP/$ggs) set WGR=WIGGLEGROUP
  if (-d  tmp/WIGGLERUN/$ggs) set WGR=WIGGLERUN
  if (! -d  tmp/$WGR/$ggs) continue
  set ok=0

# exon multipeaks
  echo "Construct the multipeaks  $WGR/$ggs"
  echo "  bin/wiggle  -multiPeaks 4 -minCover $minExonCover -I BF -O COUNT -wiggle1 tmp/$WGR/$ggs/$chrom/R.chrom.u.f.BF.gz  -wiggle2 tmp/$WGR/$ggs/$chrom/R.chrom.u.r.BF.gz -stranding 98 -o  tmp/X.$MAGIC/XH$chrom/f3.$XGH.f"
# grab the exons from the wiggle
  bin/wiggle  -multiPeaks 4 -minCover $minExonCover -I BF -O COUNT -wiggle1 tmp/$WGR/$ggs/$chrom/R.chrom.u.f.BF.gz  -wiggle2 tmp/$WGR/$ggs/$chrom/R.chrom.u.r.BF.gz -stranding 98 -o  tmp/X.$MAGIC/XH$chrom/f3.$XGH.f
  bin/wiggle  -multiPeaks 4 -minCover $minExonCover -I BF -O COUNT -wiggle1 tmp/$WGR/$ggs/$chrom/R.chrom.u.r.BF.gz  -wiggle2 tmp/$WGR/$ggs/$chrom/R.chrom.u.f.BF.gz -stranding 98 -o  tmp/X.$MAGIC/XH$chrom/f3.$XGH.r

  cat tmp/X.$MAGIC/XH$chrom/f3.$XGH.f.multiPeaks ZZZZZ | scripts/tab_sort -k 1,1 -k 2n |  gawk -F '\t' '/^#/{next;}{if(a10<1)a10=$2;if($1=="ZZZZZ" || $2 >= a2+10){if(a2-a1 > 0){color=1;s=score;while (color<7 && s>100){s/=5;color++;}if(score>=1)printf("Sequence %s_%s__%d_%d\nForward\ncDNA_clone %s_%s__%d_%d\nIntMap %s %d %d\nIs_read\nTags %d\nColour Green%d%s\n\n",XGH,c,a10,a2,XGH,c,a10,a2,c,a10,a2,int(score),color,cov);cov="";a10=$2;score=0;}}c=$1; a1=$2;a2=$3;if(a1<1)a1=1;cov= cov "\nComposite " a1 - a10 + 1 " " a2 - a10 " " int($6); if($6>score)score=int($6);} ' XGH=$XGH >  tmp/X.$MAGIC/XH$chrom/f3.$XGH.f.ace
  wc  tmp/X.$MAGIC/XH$chrom/f3.$XGH.f.ace
  cat tmp/X.$MAGIC/XH$chrom/f3.$XGH.f.ace  | gawk '/^Sequence/{split($2,aa,"__");chrom=substr(aa[1],4);split(aa[2],bb,"_");a1=bb[1];a2=bb[2];ln=a2-a1;if(ln<0)ln=-ln;ln++;printf("%s\t1\t%d\t%s\t%d\t%d\n",$2,ln,chrom,a1,a2);}' >  tmp/X.$MAGIC/XH$chrom/f3.$XGH.f.shadow
  bin/dna2dna  -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow   tmp/X.$MAGIC/XH$chrom/f3.$XGH.f.shadow | gawk '/^>/{print;next;}{printf("%s\n",$1);}' >  tmp/X.$MAGIC/XH$chrom/f3.$XGH.f.fasta


   cat tmp/X.$MAGIC/XH$chrom/f3.$XGH.r.multiPeaks ZZZZZ | scripts/tab_sort -k 1,1 -k2nr | gawk -F '\t' '/^#/{next;}{if(a10<1){a1=$3;a10=$3;}if($1=="ZZZZZ" || $3 <= a2-10){if(a2-a1 < 0){color=1;s=score;while (color<7 && s>100){s/=5;color++;}if(score>=1)printf("Sequence %s_%s__%d_%d\nForward\ncDNA_clone %s_%s__%d_%d\nIntMap %s %d %d\nIs_read\nTags %d\nColour Green%d%s\n\n",XGH,c,a10,a2,XGH,c,a10,a2,c,a10,a2,int(score),color,cov);cov="";a10=$3;score=0;}}c=$1; a1=$3;a2=$2;if(a2<1)a2=1;cov= cov "\nComposite " a10 -a1  + 1 " " a10 -a2 " " int($6); if($6>score)score=int($6); }' XGH=$XGH >  tmp/X.$MAGIC/XH$chrom/f3.$XGH.r.ace
  wc  tmp/X.$MAGIC/XH$chrom/f3.$XGH.r.ace
  cat tmp/X.$MAGIC/XH$chrom/f3.$XGH.r.ace  | gawk '/^Sequence/{split($2,aa,"__");chrom=substr(aa[1],4);split(aa[2],bb,"_");a1=bb[1];a2=bb[2];ln=a2-a1;if(ln<0)ln=-ln;ln++;printf("%s\t1\t%d\t%s\t%d\t%d\n",$2,ln,chrom,a1,a2);}' >  tmp/X.$MAGIC/XH$chrom/f3.$XGH.r.shadow
  bin/dna2dna  -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow   tmp/X.$MAGIC/XH$chrom/f3.$XGH.r.shadow | gawk '/^>/{print;next;}{printf("%s\n",$1);}' >  tmp/X.$MAGIC/XH$chrom/f3.$XGH.r.fasta

# transcriptsEnds :  

    echo "Construct the transcriptsEnds  $WGR/$ggs"
    echo "  bin/wiggle  -transcriptsEnds tmp/$WGR/$ggs/$chrom/R.chrom.u -I BF -O COUNT -o tmp/X.$MAGIC/XH$chrom/f3.Xends.$ggs -stranding 98 -minCover 100 -wiggleRatioDamper 5"
                bin/wiggle  -transcriptsEnds tmp/$WGR/$ggs/$chrom/R.chrom.u -I BF -gzi -O COUNT -o tmp/X.$MAGIC/XH$chrom/f3.Xends.$ggs -stranding 98 -minCover 100 -wiggleRatioDamper 5

  foreach fr (ELF ELR ERF ERR)
    cat tmp/X.$MAGIC/XH$chrom/f3.Xends.*.$fr.transcriptsEnds >> tmp/X.$MAGIC/XH$chrom/f3.Xends.$fr.transcriptsEnds
  end
end

foreach fr (ELF ELR ERF ERR)
  set ff=tmp/X.$MAGIC/XH$chrom/f3.Xends.$fr.transcriptsEnds
  if (-e  $ff) then

     if ($fr == ELF) then
       cat $ff | gawk -F '\t' '/^#/{next;}{if($2<1)$2=1;printf("Xends_%s.%s__%d_%d\t1\t%d\t%s\t%d\t%d\n", fr,$1,$2,$3,$4,$1,$2,$3);}' fr=$fr > $ff.shadow
       cat $ff.shadow | gawk -F '\t' '{printf("Sequence %s\ncDNA_clone %s\nIs_read\nIntMap %s %d %d\nColour CYAN\nComposite 100\nForward\n\n", $1, $1,$4,$5,$6);}' > $ff.ace
     else if ($fr == ERF) then
       cat $ff | gawk -F '\t' '/^#/{next;}{if($2<1)$2=1;printf("Xends_%s.%s__%d_%d\t1\t%d\t%s\t%d\t%d\n", fr,$1,$2,$3,$4,$1,$2,$3);}' fr=$fr > $ff.shadow
       cat $ff.shadow | gawk -F '\t' '{printf("Sequence %s\ncDNA_clone %s\nIs_read\nIntMap %s %d %d\nColour YELLOW\nComposite 100\nForward\n\n", $1, $1,$4,$5,$6);}' > $ff.ace
     else if ($fr == ERR) then
       cat $ff | gawk -F '\t' '/^#/{next;}{if($2<2)$1=1;printf("Xends_%s.%s__%d_%d\t1\t%d\t%s\t%d\t%d\n", fr,$1,$3,$2,$4,$1,$3,$2);}' fr=$fr > $ff.shadow
       cat $ff.shadow | gawk -F '\t' '{printf("Sequence %s\ncDNA_clone %s\nIs_read\nIntMap %s %d %d\nColour CYAN\nComposite 100\nForward\n\n", $1, $1,$4,$5,$6);}' > $ff.ace
     else if ($fr == ELR) then
       cat $ff | gawk -F '\t' '/^#/{next;}{if($2<1)$2=1;printf("Xends_%s.%s__%d_%d\t1\t%d\t%s\t%d\t%d\n", fr,$1,$3,$2,$4,$1,$3,$2);}' fr=$fr > $ff.shadow
       cat $ff.shadow | gawk -F '\t' '{printf("Sequence %s\ncDNA_clone %s\nIs_read\nIntMap %s %d %d\nColour YELLOW\nComposite 100\nForward\n\n", $1, $1,$4,$5,$6);}' > $ff.ace
     endif

     bin/dna2dna  -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow $ff.shadow > $ff.fasta
  endif
end


    cat   tmp/X.$MAGIC/XH$chrom/f3.X*.ace >    tmp/X.$MAGIC/XH$chrom/f3.XFC2.any.ace 
    cat   tmp/X.$MAGIC/XH$chrom/f3.X*.shadow >    tmp/X.$MAGIC/XH$chrom/f3.XFC2.any.shadow
    bin/dna2dna  -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow   tmp/X.$MAGIC/XH$chrom/f3.XFC2.any.shadow >    tmp/X.$MAGIC/XH$chrom/f3.XFC2.any.fasta


# grab the pA and the SL
  echo "XXXXXXXXXXXXXXXXXXXXXX Grab the XA_"
ls -ls   tmp/SLpA/$ggNS.SLpA.gz tmp/SLpA/$ggS.SLpA.gz
if ((-e tmp/SLpA/$ggNS.SLpA.gz || -e tmp/SLpA/$ggS.SLpA.gz) && ! -e tmp/X.$MAGIC/XH$chrom/f3.XA.fasta)  then
  echo "Grab the XA_"
  gunzip -c tmp/SLpA/$ggNS.SLpA.gz tmp/SLpA/$ggS.SLpA.gz | gawk -F '\t' '/^pA/{if($2 == chrom && $5 > 10) {a1=$3;if($4=="Forward")a2=a1-30;else a2=a1+30;support=$5;printf("Sequence XA_%s__%s_%d\ncDNA_clone XA_%s__%s_%d\nIntMap %s %d %d\nIs_read\nColour LIGHTORANGE\nComposite %d\nReverse\nmReverse\nPolyA_after_base 9\n\n",chrom,a1,a2,chrom,a1,a2,chrom,a1,a2,support);}}' chrom=$chrom >  tmp/X.$MAGIC/XH$chrom/f3.XA.ace
  cat tmp/X.$MAGIC/XH$chrom/f3.XA.ace | gawk '/^Sequence/{split($2,aa,"__");chrom=substr(aa[1],4);split(aa[2],bb,"_");a1=bb[1];a2=bb[2];ln=a2-a1;if(ln<0)ln=-ln;ln++;printf("%s\t1\t%d\t%s\t%d\t%d\n",$2,ln,chrom,a1,a2);}' >  tmp/X.$MAGIC/XH$chrom/f3.XA.shadow
  bin/dna2dna  -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow   tmp/X.$MAGIC/XH$chrom/f3.XA.shadow | gawk '/^>/{print;next;}{printf("AAAAAAAA%s\n",$1);}' >  tmp/X.$MAGIC/XH$chrom/f3.XA.fasta
ls -ls  tmp/X.$MAGIC/XH$chrom/f3.XA.ace

  echo "Grab the XSL_"
  gunzip -c tmp/SLpA/$ggNS.SLpA.gz tmp/SLpA/$ggS.SLpA.gz | gawk -F '\t' '/^SL/{if($2 == chrom && $5 > 4) {a1=$3;if($4=="Forward")a2=a1+30;else a2=a1-30;support=$5;printf("Sequence X%s_%s__%s_%d\ncDNA_clone X%s_%s__%s_%d\nIntMap %s %d %d\nIs_read\nColour BLACK\nComposite %d\nForward\nmForward\n\n",$1,chrom,a1,a2,$1,chrom,a1,a2,chrom,a1,a2,support);}}' chrom=$chrom >  tmp/X.$MAGIC/XH$chrom/f3.XSL.ace
  cat tmp/X.$MAGIC/XH$chrom/f3.XSL.ace | gawk '/^Sequence/{split($2,aa,"__");i=index(aa[1],"_");chrom=substr(aa[1],i+1);split(aa[2],bb,"_");a1=bb[1];a2=bb[2];ln=a2-a1;if(ln<0)ln=-ln;ln++;printf("%s\t1\t%d\t%s\t%d\t%d\n",$2,ln,chrom,a1,a2);}' >  tmp/X.$MAGIC/XH$chrom/f3.XSL.shadow
  bin/dna2dna  -i TARGET/CHROMS/$species.chrom_$chrom.fasta.gz -shadow   tmp/X.$MAGIC/XH$chrom/f3.XSL.shadow | gawk '/^>/{print;next;}{printf("%s\n",$1);}' >  tmp/X.$MAGIC/XH$chrom/f3.XSL.fasta
else
  touch  tmp/X.$MAGIC/XH$chrom/f3.XSL.ace tmp/X.$MAGIC/XH$chrom/f3.XSL.fasta
endif

laba:

  bin/tacembly  tmp/X.$MAGIC/XH$chrom << EOF
    read-models
    parse tmp/TABIX/$MAGIC/$chrom.tabix.ace
    pparse MetaDB/$MAGIC/runs.ace
    pparse MetaDB/$MAGIC/samples.ace
    pparse MetaDB/$MAGIC/ali.ace
    // pparse my.chrom.alias.ace
    pparse tmp/X.$MAGIC/EHITS/$chrom/f1.ace.gz
    pparse tmp/X.$MAGIC/EHITS/$chrom/f1.fasta.gz
    pparse tmp/X.$MAGIC/XH$chrom/f3.XG.f.ace 
    pparse tmp/X.$MAGIC/XH$chrom/f3.XG.f.fasta
    pparse tmp/X.$MAGIC/XH$chrom/f3.XG.r.ace 
    pparse tmp/X.$MAGIC/XH$chrom/f3.XG.r.fasta
    pparse tmp/X.$MAGIC/XH$chrom/f3.XH.f.ace 
    pparse tmp/X.$MAGIC/XH$chrom/f3.XH.f.fasta
    pparse tmp/X.$MAGIC/XH$chrom/f3.XH.r.ace 
    pparse tmp/X.$MAGIC/XH$chrom/f3.XH.r.fasta
    pparse tmp/X.$MAGIC/XH$chrom/f3.XFC2.any.ace 
    pparse tmp/X.$MAGIC/XH$chrom/f3.XFC2.any.fasta
    pparse tmp/X.$MAGIC/XH$chrom/f3.allDoubleIntrons.ace
    pparse tmp/X.$MAGIC/XH$chrom/f3.allDoubleIntrons.fasta
    pparse tmp/X.$MAGIC/XH$chrom/f3.allDoubleIntronsGenomic.ace
    pparse tmp/X.$MAGIC/XH$chrom/f3.allDoubleIntronsGenomic.fasta
    pparse tmp/X.$MAGIC/XH$chrom/f3.XA.ace
    pparse tmp/X.$MAGIC/XH$chrom/f3.XA.fasta
    pparse tmp/X.$MAGIC/XH$chrom/f3.XSL.ace
    pparse tmp/X.$MAGIC/XH$chrom/f3.XSL.fasta
    // pparse tmp/X.$MAGIC/XH$chrom/f3.XLRFR.ace
    // pparse tmp/X.$MAGIC/XH$chrom/f3.XLRFR.fasta
    pparse tmp/X.$MAGIC/EHITS/$chrom/f3.introns.ace
    pparse tmp/X.$MAGIC/EHITS/$chrom/f3.transcriptsIntronSupport.ace
    query find Run IS $ggs && ! W_colour_plus
    edit W_colour_plus BLUE
    edit W_colour_minus GREEN
    save
    quit
EOF

if (-e   tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.ace) then

  bin/tacembly  tmp/X.$MAGIC/XH$chrom << EOF
    pparse tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.ace
    pparse tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.other.ace
    pparse tmp/X.$MAGIC/EHITS/$chrom/f3.intronDB.other.fasta
    save
    quit
EOF
endif



hello:

# these lines add the RefSeq mapping to this chromosome as additional pseudo composite reads
if (1) then
 
  cat tmp/METADATA/RefSeq.mrna_map_ln_gc_gene_geneid.txt | gawk -F '\t' '{gsub(/\"/,"",$0); split($2,aa,":");if(aa[1] != chrom)next;split(aa[2],bb,"-");m=$1;c[m]=aa[1];a1[m]=bb[1];a2[m]=bb[2];}END{for(m in c)printf("Sequence %s\nIntMap %s %d %d\nForward\nComposite 10\nIs_read\n\n",m,c[m],a1[m],a2[m]);}' chrom=$chrom | gzip > tmp/X.$MAGIC/XH$chrom/f3.RefSeq.intmap.ace.gz

  gunzip -c tmp/X.$MAGIC/XH$chrom/f3.RefSeq.intmap.ace.gz ZZZZZ.gz $ici/TARGET/Targets/$species.RefSeq.fasta.gz |   gawk '/^ZZZZZ/{zz++;next;}/^Sequence/{okk[$2]=1;next;}{if(zz<1)next;}/^>/{s=substr($1,2);i=index(s,"|");if(i>0)s=substr(s,1,i-1);ok=0;if(okk[s]==1){ok=1;print ">" s;}next;}{if(ok==1)print}' | gzip > tmp/X.$MAGIC/XH$chrom/f3.RefSeq.fasta.gz

  bin/tacembly  tmp/X.$MAGIC/XH$chrom << EOF
    pparse tmp/X.$MAGIC/XH$chrom/f3.RefSeq.intmap.ace.gz
    pparse tmp/X.$MAGIC/XH$chrom/f3.RefSeq.fasta.gz
EOF

endif

  bin/tacembly  tmp/X.$MAGIC/XH$chrom << EOF
    query find sequence is_read && ! IS X?_* && ! intmap==$chrom
    edit -D is_read
    find method encode
    edit Show_up_strand
    query find sequence XE_*
    // edit -D Is_read
    query find sequence IS XY_* && (Other  || ct_ac)
    edit -D Is_read
     save
    quit
EOF

  bin/tacembly  tmp/X.$MAGIC/XH$chrom << EOF
    query find est XC_*_10
    edit colour green1
    query find est XC_*_20
    edit colour green2
    query find est XC_*_50
    edit colour green3
    query find est XC_*_100
    edit colour green4
    query find est XC_*_200
    edit colour green5
    query find est XC_*_500
    edit colour green6
    query find est XC_*_1000
    edit colour green7
     query find est XF_*_10
    edit colour blue1
    query find est XF_*_20
    edit colour blue2
    query find est XF_*_50
    edit colour blue3
    query find est XF_*_100
    edit colour blue4
    query find est XF_*_200
    edit colour blue5
    query find est XF_*_500
    edit colour blue6
    query find est XF_*_1000
    edit colour blue7
    save
    find sequence Is_read && ! cdna_clone
    acem
      cdna_80
      quit
    query find sequence XD_* 
    edit -D Is_read
    query find run w_stranded
    edit W_colour
    query find run w_new_exon
    edit W_colour
    save
    quit
EOF

  scripts/f3.kill_ct_ac_introns.tcsh $chrom 1

# transfer the RNA_seq support read in XI->composite back into the Intron class, which is visible in the tg
bin/tacembly tmp/X.$MAGIC/XH$chrom <<EOF
  select -o tmp/X.$MAGIC/XH$chrom/f3.intron_support.txt  xi, ii, n from xi in ?Sequence where xi ~ "XI_*", n in xi->composite, ii in xi->intron where n && ii
  date
EOF
cat tmp/X.$MAGIC/XH$chrom/f3.intron_support.txt | gawk -F '\t' '{printf ("Intron %s\nRNA_seq %d\n\n",$2,$3);}' >  tmp/X.$MAGIC/XH$chrom/f3.intron_support.ace
bin/tacembly tmp/X.$MAGIC/XH$chrom <<EOF
  parse tmp/X.$MAGIC/XH$chrom/f3.intron_support.ace
  save
  quit
EOF


##### parse the RefSeq

foreach target ($RNAtargets)
  if ($target == av) continue
  if ($target == rrna) continue
  if (-e tmp/METADATA/gtf.$target.transcripts.ace.gz) then
    bin/tacembly tmp/X.$MAGIC/XH$chrom <<EOF  
      query find predicted_gene method == $target
      kill
      pparse TARGET/GENES/RvY.LocusLink.ace
      pparse tmp/METADATA/gtf.$target.transcripts.ace.gz 
      query is_predicted_gene
      edit -D method
      edit method $target
      spush
      query intmap == $chrom
      sminus
      spop
      kill
      find method $target
      edit Colour GREEN
      edit Show_up_strand
      edit Bumpable
      edit Right_priority 1.6
      save
      quit
EOF
  endif
end


# les introns mangent 10% de leur soutien sur les zones inclues des XH sur les 2 brins

  touch tmp/X.$MAGIC/XH$chrom/f3.parse.done

endif

exit 0

echo 284702 > titi
echo ZZZZZ >> titi

gunzip -c PUBLIC/SEQC_NB_MAV_G_log2.20121127.txt.gz | grep  284702 >> titi


cat titi | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$1]=1;next;}}{n=split($2,aa,";");for(i=1;i<=n;i++)if(gg[aa[i]]==1)print;}' 


gunzip -c  titi.gz  ZZZZZ.gz PUBLIC/SEQC_NB_MAV_G_log2.20121127.txt.gz


  gunzip -c RESULTS/MAV_lostGeneID.txt.gz ZZZZZ.gz PUBLIC/SEQC_NB_MAV_G_log2.20121127.txt.gz | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$1]=1;next;}}{n=split($2,aa,";");for(i=1;i<=n;i++)if(gg[aa[i]]==1)print;}' >> RESULTS/MAV_lostGeneID.SEQC_NB_MAV_G_log2.20121127.txt



cat tatou.17.noncoding.list ZZZZZ tatou.201.coding.list ZZZZZ tatou.29.spliced.list ZZZZZ tatou.6.pseudogenes.list | gawk '/^ZZZZZ/{zz++;next;}/^Sequence/{if(zz==0)g[$2]=g[$2] " non_coding";if(zz==1)g[$2]=g[$2] " coding";if(zz==2)g[$2]=g[$2] " spliced";if(zz==3)g[$2]=g[$2] " pseudogene";}END{for(gg in g)printf("%s\t%s\n",gg,g[gg]);}' | sed -e 's/\"//g' -e 's/X__//' | gawk '{printf("%s\t",$1);if(index($0,"coding")>0)nam="Coding"; else nam="Non_coding" ; if(index($0,"spliced")>0)nam= nam "_spliced"; else nam = nam "_single_exon" ;printf("%s\n",nam);}' | sort


  
gunzip -c RESULTS/MAV_lostGeneID.txt.gz ZZZZZ.gz > tutu
cat tutu tmp/METADATA/av.metadata.txt | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz<1){gg[$1]=1;next;}}{n=split($3,aa,";");for(i=1;i<=n;i++)if(gg[aa[i]]==1)print;}'


# table transcripts2spliced_RefSeq_coding.def used on ace11
scp transcripts2spliced_RefSeq_coding.def ace11:SERVER

cat transcripts2spliced_RefSeq_coding.txt | gawk -F '\t' '{gsub(/\"/,"",$0);g=$1;gg[g]=1;if($2 != "NULL")spl[g]=1;if($3 != "NULL" && index(rs[g],$3)==0)rs[g]=rs[g] "_" $3;if($4 != "NULL")coding[g]=1;if($5 != "NULL")gene[g]=$5;}END{for (g in gg){if(coding[g]==1) nam="Coding"; else nam="Non_coding" ; if(spl[g]==1)nam= nam "_spliced"; else nam = nam "_single_exon" ; if(length(rs[g])>1)nam= nam "_with_RefSeq" ; printf("%s\t%s\t%s\t%s\n",gene[g],g,nam,substr(rs[g],2));}}' | sort | gzip > transcripts2spliced_RefSeq_coding.txt.1.gz

gunzip -c TARGET/Targets/hs.av.fasta.gz  ZZZZZ.gz transcripts2spliced_RefSeq_coding.txt.1.gz | gawk '/^ZZZZZ/{zz++;next;}/^#/{next;}/^>/{split(substr($1,2),aa,"|");s=aa[1];ss[s]=1;next;}{if(zz<1)next;}{k=0;s=$2;if(ss[s]==1)k="kept";else k="rejected";printf("%s\t",k);print;}' > tutu
cat tutu | sort > RESULTS/transcripts2spliced_RefSeq_coding.kept_rejected.txt

##### expression of a few genes
pushd RESULTS/Expression/unique/RefSeq

echo ZZZZZ > ZZZZZ
foreach ff (`ls $MAGIC.*.txt  $MAGIC.*.txt.gz `)
  gunzip -c -f SELECTED/gene.list ZZZZZ $ff | gawk -F '\t' '/^ZZZZZ/{zz++;next;}{if(zz+0<1){gg[$1]=1;next;}}/^#/{print;next;}{if(gg[$1] + gg[$2]> 0)print;next;}' > SELECTED/$ff
end

  
done:
echo 'done'

exit 0

set ggNS=zozo
foreach chrom ($chromSetAll)
  mv tmp/TABIX/$ggNS/$chrom.u0.f.tabix.gz    tmp/TABIX/$ggNS/$chrom.u.f.tabix.gz
  mv tmp/TABIX/$ggNS/$chrom.u0.r.tabix.gz     tmp/TABIX/$ggNS/$chrom.u.r.tabix.gz
  mv tmp/TABIX/$ggNS/$chrom.u0.f.tabix.gz.tbi tmp/TABIX/$ggNS/$chrom.u.f.tabix.gz.tbi
  mv tmp/TABIX/$ggNS/$chrom.u0.r.tabix.gz.tbi tmp/TABIX/$ggNS/$chrom.u.r.tabix.gz.tbi
end



