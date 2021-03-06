/* this keyset contains names that if they match will be
 * refused as Blast_pfam_kantor titles see kantor.c 
 * isHorribleTitle
 * a PFAM is then done behind 
 */

static char *addLikeToTitle [] = {
"*disease*" ,
"*leukemia*" ,
"*metasta**" ,
"*carcinoma*" ,
"*lymph**" ,
"* node *" ,
"*T-cell*" ,
"*antenna*" ,
"* eye *" ,
"*kidney*" ,
"*liver *" ,
"*renal*" ,
"*craniofacial*" ,
"*placenta**" ,
"*palate*" ,
"*lung*" ,
"*pancreas*" ,
"*stomach*" ,
"*blood*" ,
"*antigen*" ,
"*induced*" ,
"*binding*" ,
"*pneumocyst**" ,
"*faecalis*" ,
"*streptococ**" ,
"*lavendulae*" ,
"*volvulus*" ,
"*african*" ,
"*frog*" ,
"* pig *" ,
"*mammal*" ,
"*syndrome*" ,
"*yeast*" ,
"*tumor*" ,
"*cancer*" ,
"*drosophila*" ,
"*thaliane*" ,
"*melano*" ,
"*droso*" ,
"*arabidopsis*" ,
"*plant*" ,
"* rat *" ,
"*mouse*" ,
"*human*" ,
"*sapiens*" ,
"*pombe",
"* kD*" ,
"*cerevisiae*",
0
} ;



static char *horribleTitlesToKill [] = {
  "* _ *",
  "receptor ?",
  "*cdna sequence*",
  "nuclear",
  "acid",
  "containing source",
  "*gene near*",
  "*chromosome*",
  "*open*reading frame*",
  "*domain containing*",
"family member" ,
  "*caeel*",
"-like 2" ,
". 10 100" ,
"3-7 gene product" ,
"4 superfamily member 9" ,
"16.7 kDa-like protein" ,
"21D7" ,
"32.0k protein Neurospora crassa" ,
"59.2KD protein IN PFK26-SGA1 INTERGENIC REGION" ,
"a C.elegans protein encoded in cosmid" ,
"abundant gene transcript" ,
"AD025" ,
"AK001373_like" ,
"alternative function" ,
  "*alu subfamily sequence contamination*" ,
  "*warning entr*",
"antigen" ,
"AOF1001" ,
"associated 1" ,
"associated protein" ,
"At1g04290 F19P19_27" ,
"At1g07990 T6D22_5" ,
"At1g24120 F3I6_4" ,
"AT3g20810 MOE17_10" ,
"At4g24270 T22A6_100" ,
"At4g32910 F26P21_30" ,
"AT4g36980 C7A10_380" ,
"AT5g20600 F7C8_190" ,
"AT5g51140 MWD22_8" ,
"AT5g65760 MPA24_11" ,
"B-cell CLL lymphoma 7A" ,
"BC85722_1" ,
"C. elegans protein" ,
"C. Elegans protein F17C8.5" ,
"CeCHK2" ,
"CeCRMP DHP-1" ,
"CeMef-2" ,
"CG5805 gene product" ,
"CG5986 PROTEIN. 5 100" ,
"CG12113 gene product" ,
"CG15432 gene product" ,
  "chromo" ,
"chromosome 1 BAC protein" ,
"chromosome 4 open reading frame 1" ,
"chromosome 11 open reading frame 23" ,
"chromosome 20 open reading frame 77" ,
"chromosome 20 open reading frame 81" ,
"chromosome 20 open reading frame 135" ,
"chromosome 22 open reading frame 5" ,
"chromosome II BAC T27A16 sequence" ,
"chromosome protein" ,
"conserved ORF" ,
"conserved protein" ,
"conserved protein." ,
"cosmid C01G5" ,
  "crystal structure of" ,
"cytosolic factor" ,
"cytosolic protein" ,
"D. melanogaster protein" ,
"D.melanogaster" ,
"data source protein" ,
"DDBJ Accession Number restin CLIP-170" ,
"DKFZP434B168-like protein" ,
"DNA segment Chr 1 Wayne State University 40 source" ,
"DNA segment Chr 9 Wayne State University 149 source" ,
"F23149_1" ,
"family protein" ,
"for expression of" ,
"gene 17.5 protein chicken" ,
"gene product" ,
"I F48E3.1" ,
"intrinsic membrane source" ,
"involved IN plasmid" ,
"involved in plasmid" ,
"KIAA0406-like protein" ,
  "long" ,
"membrane protein" ,
"MEMBRANE protein" ,
"mRNA for protein" ,
"Native 66kDa protein" ,
"not really started" ,
"novel protein" ,
"open reading frame 6" ,
"ORF identified by homology. See FEBS Letters" ,
"Pfam family PF00145" ,
"Pfam PF00001" ,
"Pfam PF00071" ,
"Pfam PF01498" ,
"poly" ,
"Poly" ,
"related protein" ,
"predicted protein" ,
"protein of unknown function" ,
"protein -related" ,
"protein -like" ,
"protein" ,
"protein ?"
"protein 2BE2121 Chinese hamster" ,
"PROTEIN 15E1.2" ,
"protein African clawed frog" ,
"protein B0280.5 III" ,
"protein BB" ,
"protein C04A2.7 II" ,
"Protein CGI-126" ,
"protein from EUROIMAGE 1987170" ,
"protein G01D9.2 Caenorhabditis briggsae" ,
"protein G01D9.4 Caenorhabditis briggsae" ,
"protein G01D9.5 Caenorhabditis briggsae" ,
"protein homolog" ,
"Protein homolog" ,
"protein IN IKI1-ERG9 INTERGENIC REGION" ,
"Protein like" ,
"protein Podospora anserina" ,
"protein R31180_1" ,
"protein R31449_3" ,
"protein T19C3.1 III" ,
"protein T30N20_280" ,
"pseudo", 
"putatative protein" ,
"R27216_1" ,
"R29893_1" ,
"R32611_1" ,
"R32611_2" ,
"reserved" ,
"source" ,
"specific protein" ,
"The ha0919 gene product is novel." ,
"The ha1009 gene product is novel." ,
"This gene is novel." ,
"transmembrane protein" ,
"TRANSMEMBRANE protein" ,
"Uncharacterised protein family" ,
"Uncharacterized ACR, COG2106" ,
"Uncharacterized ACR, YagE family COG1723" ,
"x0001" ,
"sequence receptor delta" ,
"RNA recognition motif" ,
"contains C2H2-type zinc finger" ,
"nematode specific ORF" ,
"transmembranous domain*" ,
"possible trans-splice site at*" ,
"Leucine Rich Repeat" ,
"AT5g20600 F7C8_190" ,
"repeat * protein" ,
"Homeobox domain*" ,
"Protein CGI-126" ,
"protein CGI-*" ,
"Pumilio-family RNA binding domains" ,
"Structure Tsg101 Uev Domain" ,
"complex-associated protein" ,
"contains EST protein" ,
"contains EST*" ,
"antigen" ,
"*ankyrin motif*" ,
"involved IN plasmid" ,
"*ankyrin repeat*" ,
"ankyrin " ,
"conserved protein" ,
"Coiled coil protein" ,
"Protein kinase" ,
"mRNA product" ,
"*Zinc finger protein" ,
"helicase domain*" ,
"protein IN IKI1-ERG9 INTERGENIC REGION" ,
"protein IN * INTERGENIC REGION" ,
"*intergenic region*" ,
"EGF domain*" ,
"EGF-domain*" ,
"gene product" ,
"Coiled coil*" ,
"CG15432 gene product" ,
"CG* gene product" ,
"protein R12B2.5 III" ,
"protein *.* I*" ,
"delta-like 4" ,
"Src homology domain*" ,
"binding protein" ,
"helicase" ,
"gene from NF2 meningioma region of 22q12" ,
"gene from *region*" ,
"TcHSLR3" ,
"coiled-coil*" ,
"protein transport" ,
"member of worm-specific protein family" ,
"*protein kinase domain" ,
"HSPC335" ,
"DNA binding protein" ,
"*tetratricopeptide repeat domain*" ,
"resistance protein" ,
"protein 2BE2121 Chinese hamster" ,
"protein * chinese hamster" ,
"NUCLEIC ACID BINDING*" ,
"SER THR-protein kinase" ,
"protein kinase like" ,
"involved in plasmid*" ,
"*type zinc finger*" ,
"putatative protein" ,
"protein G01D9.4 Caenorhabditis briggsae" ,
"protein * Caenorhabditis briggsae" ,
"dehydrogenase" ,
"bZIP DNA binding domain*" ,
"*DNA binding domain*" ,
"protein R12B2.5 III" ,
"guanine nucleotide binding domain*" ,
"DNAJ protein" ,
"secreted protein" ,
"*expressed gene*" ,
"*elegans*" ,
"*caenorhabditis*" ,
"*briggsae*" ,
"*WD40 repeat*" ,
  "ribosomal protein*",
  "hypothetical protein", 
0
} ;

