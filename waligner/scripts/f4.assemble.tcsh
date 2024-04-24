#!bin/tcsh -f

set chrom=$1
setenv ici `pwd`

echo -n "f4.assemble.tcsh start "
date

if ( -e tmp/X.$MAGIC/XH$chrom/f3.parse.done) then

  pushd tmp/X.$MAGIC/XH$chrom
    if (-e TABIX) \rm TABIX
    ln -s $ici/tmp/TABIX TABIX
    if (! -e tables) ln -s ../../metaData/tables
  popd
  if ($species == worm && ! -e tmp/X.$MAGIC/EHITS/$chrom/genes.ace) then
    tbly ~/yknew <<EOF
      query find gene transcribed_gene && IntMap == $chrom
      show -a -f tmp/X.$MAGIC/XH$chrom/f4.genes2intmap.ace IntMap
EOF
    cat  tmp/X.$MAGIC/XH$chrom/f4.genes2intmap.ace | gawk '/^Gene/{g=$2;next;}/^IntMap/{printf("Sequence %s\nGenes %s %s %s\n\n", $2, g, $3,$4);}' >   tmp/X.$MAGIC/XH$chrom/f4.geneParents.ace
    bin/tacembly tmp/X.$MAGIC/XH$chrom << EOF
      read-models
      parse tmp/X.$MAGIC/XH$chrom/f4.genes2intmap.ace
      parse tmp/X.$MAGIC/XH$chrom/f4.geneParents.ace
      save
      quit
EOF
  
  endif

zcat   tmp/X.$MAGIC/EHITS/$chrom/f1.txts.gz  | gawk -F '\t' '/^#/{next}/^\/\//{next}{chrom=$1;gsub(/CHROMOSOME_/,"",chrom);u1=$8;u2=$9;typeb=$4;if(typeb == "Exon") { next;xx="XJ_" ;} else {xx="XI_" ; }nam=xx group "_" chrom "__" u1 "_" u2 ; names[nam]++; n = names[nam] ;nam = nam "." n ; printf(">%s\n%s%s%s\n",nam, $12,$13,$14) ;}' group=$MAGIC  | gzip > tmp/X.$MAGIC/EHITS/$chrom/f1.fasta.gz

  bin/tacembly tmp/X.$MAGIC/XH$chrom << EOF
    read-models
    parse tmp/X.$MAGIC/EHITS/$chrom/f1.fasta.gz
    query find sequence Is_read && ! cdna_clone
    acembly
      cdna_80
      quit
    acembly
      cdna_1 $chrom
      quit
    save

    table -o tmp/X.$MAGIC/XH$chrom/f4.10.polyAsuspect1.txt  -f tables/10.polyAsuspect1.def
    table -o tmp/X.$MAGIC/XH$chrom/f4.10.polyAsuspect2.txt  -f tables/10.polyAsuspect2.def
    table -o tmp/X.$MAGIC/XH$chrom/f4.10.polyAsuspect4.txt  -f tables/10.polyAsuspect4.def
    quit
EOF

laba:
## remove echo introns
## we need the intron coordinates, some are missing
  echo "find intron NOT IntMap"
  bin/tacembly  tmp/X.$MAGIC/XH$chrom << EOF
     query find intron ! IntMap
     list -a -f  tmp/X.$MAGIC/XH$chrom/f4.intron_nomap.list
     quit
EOF
  touch  tmp/X.$MAGIC/XH$chrom/f4.intron_nomap.list
  ls -ls  tmp/X.$MAGIC/XH$chrom/f4.intron_nomap.list

  cat  tmp/X.$MAGIC/XH$chrom/f4.intron_nomap.list | gawk '/^Intron/{z=$2;gsub(/\"/,"",z);n1=split (z,aa,"__");n2=split(aa[2],bb,"_") ;printf("Intron %s\nIntMap %s %s %s\n\n",$2,aa[1],bb[1],bb[2]);}' > tmp/X.$MAGIC/XH$chrom/f4.intron_nomap.ace
  bin/tacembly  tmp/X.$MAGIC/XH$chrom << EOF
     pparse  tmp/X.$MAGIC/XH$chrom/f4.intron_nomap.ace
     save
     quit
EOF

endif
if (-e tmp/X.$MAGIC/XH$chrom/database/lock.wrm) exit 1
# kill echo ct_ac introns
scripts/f3.kill_ct_ac_introns.tcsh $chrom 2
# find all antisence gene, with non classic or ct_ac intron supported 5 times less than an approximately antisense g[tc]_ag intron
  bin/tacembly  tmp/X.$MAGIC/XH$chrom << EOF
     table -o tmp/X.$MAGIC/XH$chrom/f4.killEchoIntron.out -f tables/f4.killEchoIntron.def
     quit
EOF
# export list of reads and gene
cat tmp/X.$MAGIC/XH$chrom/f4.killEchoIntron.out | gawk -F '\t' '/^"/{printf("Sequence %s\n",$3);printf("Transcribed_gene %s\n",$1);}' > tmp/X.$MAGIC/XH$chrom/f4.killEchoIntron.list
# kill the bad reads, recompute the genes, they will usually vanish
echo "kill tmp/X.$MAGIC/XH$chrom/f4.killEchoIntron.list"
if (-e tmp/X.$MAGIC/XH$chrom/database/lock.wrm) exit 1

  bin/tacembly  tmp/X.$MAGIC/XH$chrom << EOF
     key  tmp/X.$MAGIC/XH$chrom/f4.killEchoIntron.list
     spush
     query find sequence IS XY_* && (Other  || ct_ac)
     sor
     follow from_gene
     sor
     spop
     spush
     query CLASS Sequence
     edit -D Is_read
     spop
     query CLASS Transcribed_gene
     list -a -f  tmp/X.$MAGIC/XH$chrom/f4.killEchoIntron.tg.list
     save
     quit
EOF
if (-e tmp/X.$MAGIC/XH$chrom/database/lock.wrm) exit 1
# split this code out of previous one to monitor eventual errors
  bin/tacembly  tmp/X.$MAGIC/XH$chrom << EOF
     key  tmp/X.$MAGIC/XH$chrom/f4.killEchoIntron.tg.list
     acem
        cdna_73
        quit
     query find mrna DNA:2 > 10000 ; > from_gene
     acem
        cdna_73 // 2017_03_11 a hack to rm stupid extra long empty exons
        quit
     save
     find clone
     list -a -f tmp/X.$MAGIC/XH$chrom/f4.killEchoIntron.done
     quit
EOF
if (-e tmp/X.$MAGIC/XH$chrom/database/lock.wrm) exit 1
echo -n "f4.killEchoIntron done"
date
scripts/f3.kill_ct_ac_introns.tcsh $chrom 3

touch  tmp/X.$MAGIC/XH$chrom/f4.assemble.done
if (-e tmp/X.$MAGIC/XH$chrom/database/lock.wrm) exit 1

# suspect 1 finds reverse read starting inside the ORF
# suspect 2 finds forward polyA read endding inside the ORF
# suspect 4 removes polyA of forward read if a reverse read assembles further down

pushd tmp/X.$MAGIC/XH$chrom
  gawk  -F '\t' '/\"/ {printf ("Sequence %s\nColour PALEGREEN\n\ncDNA_clone %s\nInternal_priming\n\n",$8, $9);}'  f4.10.polyAsuspect1.txt >!  f4.10.polyAsuspect1.ace
  gawk  -F '\t' '/\"/ {printf ("Sequence %s\nColour PALEGREEN\n\ncDNA_clone %s\nInternal_priming\n\n",$8, $9);}'  f4.10.polyAsuspect2.txt >!  f4.10.polyAsuspect2.ace
  gawk  -F '\t' '/\"/ {printf ("Sequence %s\nColour PALEGREEN\n-D polya_after_base\n-D Number_of_terminal_A\n\n",$4);}'  f4.10.polyAsuspect4.txt >!  f4.10.polyAsuspect4.ace

  wc f4.10.polyAsuspect*.txt

  $ici/bin/tacembly . << EOF >! f4.10.reverse.log
    pparse f4.10.polyAsuspect1.ace
    spush 
    pparse f4.10.polyAsuspect2.ace
    sor
    pparse f4.10.polyAsuspect4.ace
    sor
    query find est Number_of_terminal_A && !PolyA_after_base
    edit -D  Number_of_terminal_A // since they were removed by the above scripts
    sor
    query find cdna_clone double_fuzzy
    query follow read ; NOT ref_mrna && NOT ref_seq
    edit -D Is_read
    sor
   
    comment "restore missing flags gt_ag in est"
    query find est Flipped
    edit -D Flipped
    query find tg ct_ac
    acem
      cdna_73  // restore missing flags gt_ag in est
      quit
    comment "FlipAllGenes"
    query find est Is_read && Reverse &&  (IS XE_* || IS XF_*)
    edit -D Is_read  // because they were doubly defined
    query find est  (Reverse && (IS XE_* || IS XF_*)) ; from_gene
    follow from_gene
    comment "realign"
    acem
      cdna_73
      quit
    acembly
      cDNA_FlipAllGenes
      quit
    query find est Is_read && Reverse && (IS XE_* || IS XF_*)
    edit -D Is_read  // because they were doubly defined
    query find est  (Reverse && (IS XE_* || IS XF_*)) ; from_gene
    follow from_gene
    comment "realign"
    acem
      cdna_73
      quit
    query find est Flipped
    sor
    spop
    query follow from_gene
    comment "realign, will kill all the non_best_mrna and resize"
    acem
      cdna_73 // will kill all the non_best_mrna and resize and fuse a first time
      quit 
    find tg
    comment cDNA_Flag_suspected_internal_deletion
    acem
      cDNA_Flag_suspected_internal_deletion -ignore // ignore non informative clones
      quit
    query find tg to_be_fused_with
    comment "fuse_genes"
    acem
      cdna_73 // -locally incompatible avec to_be_fused_with
      quit
    query find mrna NOT from_gene
    spush
    follow dna
    sor
    spop
    kill
    save
    acembly
      cdna_50
      gene_intron
      quit
    save
    quit
EOF

  scripts/f3.kill_ct_ac_introns.tcsh $chrom 4
  scripts/f3.kill_ct_ac_introns.tcsh $chrom 5
  scripts/f3.kill_ct_ac_introns.tcsh $chrom 6

#  /home/mieg/bin/gene2gene.2007_12_17 . // obsolete, replaced by cdna_21 cdna_50

  $ici/bin/tacembly . << EOF
    acembly
      cdna_21
      cdna_50
      quit
    save
    find intron
    list -a -f f4.introns.list
    find clone
    list -a -f f4.assemble.done
    quit
EOF

  echo ' ' >  f4.introns.preace
  echo ' ' >  f4.introns.ace
  if (-d $ici/GeneIndexDB) then
    $ici/bin/tacembly $ici/GeneIndexDB << EOF
      key f4.introns.list
      show -a -f f4.introns.preace
      quit
EOF

    cat  f4.introns.preace | gawk '/^$/{print}{gsub(/\"/,"",$0);}/^Intron/{printf("Intron %s\n",$2);next;}/^Group/{if($2 == nx)print;}/Validated_u/{if($2 == "any")print}/de_[du][nu]o/{if($2 == "any")print}' nx=Rat  > f4.introns.ace
  endif

  $ici/bin/tacembly . << EOF
    pparse f4.introns.ace
    find transcribed_gene
    list -a -f f4.tg.list
    query find transcribed_gene Intron_boundaries
    list -a -f f4.spliced_tg.list
    save
    quit
EOF

popd

# edit away the ct_ac XY and XI
scripts/f4.fix.tcsh $chrom




echo -n "f4.assemble.tcsh end "
date

exit 0

bin/tace ~/yk <<EOF
  date
  table -f $ici/tables/f5.tg2genebox2product.def -o  t1
  date
  table -f $ici/tables/f5.tg2genebox2product.def -bqlo  t2
  date
  table -f $ici/tables/f5.tg2genebox2product.def -o  t3
  date
  table -f $ici/tables/f5.tg2genebox2product.def -bqlo  t4
  date
EOF
wc t[1234]
ls -ls t?
