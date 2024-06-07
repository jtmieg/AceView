#!/bin/tcsh -f
# export the snps/tsnps using snpsummary

set zone=$1

if (! -e tmp/TSNP_DB/$zone/mainClone.ace) then
  cat <<EOF > tmp/TSNP_DB/$zone/mainClone.ace
Clone R
Main_clone
MainTitle  "SNPs $zone"
Species $species

EOF
  bin/tace tmp/TSNP_DB/$zone <<EOF
    read-models
y
    pparse tmp/TSNP_DB/$zone/mainClone.ace
    save
    quit
EOF
endif



if (-d DanLi) then
    bin/snpsummary -db tmp/TSNP_DB/$zone -o RESULTS/SNV/DanLi.$zone --snpType 3 
endif


  bin/snpsummary -db tmp/TSNP_DB/$zone -p $MAGIC -o RESULTS/SNV/$MAGIC.$zone.min$minSnpFrequency --snpType 0 # --minSnpFrequency $minSnpFrequency --minSnpCover $minSnpCover
  bin/snpsummary -db tmp/TSNP_DB/$zone -p $MAGIC -o RESULTS/SNV/$MAGIC.$zone.min80 --snpType 0 --minSnpFrequency 80 # --minSnpCover $minSnpCover 
  touch tmp/TSNP_DB/$zone/$MAGIC.tsnp10.done


