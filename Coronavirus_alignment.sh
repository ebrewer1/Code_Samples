#Name: Edmond Brewer
#Date: 4/27/20
#Class: Genomics BMMB 551
#Assignment: Final Project

# Make references folder and setup genome/gff downloads
set -uex
mkdir refs
cd refs

#COVID-19
mkdir covid
cd covid
COV=NC_045512.2
#Get genbank genome and convert to FASTA and GFF, then index genome
efetch -db nucleotide -format gb -id $COV > $COV.gb
cat $COV.gb | seqret -filter -feature -osformat fasta -ofname2 $COV.gff > $COV.fa
rm -R *.gb
samtools faidx $COV.fa
cd ..

#Severe Acute Respiratory Syndrome 1
mkdir sars
cd sars
SAR=NC_004718.3
#Get genbank genome and convert to FASTA and GFF, then index genome
efetch -db nucleotide -format gb -id $SAR > $SAR.gb
cat $SAR.gb | seqret -filter -feature -osformat fasta -ofname2 $SAR.gff > $SAR.fa
rm -R *.gb
samtools faidx $SAR.fa
cd ..

#Middle East Respiratory Syndrome
mkdir mers
cd mers
MER=NC_019843.3
# Get genbank genome and convert to FASTA and GFF, then index genome
efetch -db nucleotide -format gb -id $MER > $MER.gb
cat $MER.gb | seqret -filter -feature -osformat fasta -ofname2 $MER.gff > $MER.fa
rm -R *.gb
samtools faidx $MER.fa
cd ..

#Human Coronavirus HKU1
mkdir hku1
cd hku1
HKU=NC_006577.2
# Get genbank genome and convert to FASTA and GFF, then index genome
efetch -db nucleotide -format gb -id $HKU > $HKU.gb
cat $HKU.gb | seqret -filter -feature -osformat fasta -ofname2 $HKU.gff > $HKU.fa
rm -R *.gb
samtools faidx $HKU.fa
cd ..

#Human Coronavirus 229E
mkdir 229e
cd 229e
two29E=NC_002645.1
# Get genbank genome and convert to FASTA and GFF, then index genome
efetch -db nucleotide -format gb -id $two29E > $two29E.gb
cat $two29E.gb | seqret -filter -feature -osformat fasta -ofname2 $two29E.gff > $two29E.fa
rm -R *.gb
samtools faidx $two29E.fa
cd ..

#Human Coronavirus OC43
mkdir oc43
cd oc43
OC43=NC_006213.1
# Get genbank genome and convert to FASTA and GFF, then index genome
efetch -db nucleotide -format gb -id $OC43 > $OC43.gb
cat $OC43.gb | seqret -filter -feature -osformat fasta -ofname2 $OC43.gff > $OC43.fa
rm -R *.gb
samtools faidx $OC43.fa
cd ..

#Human Coronavirus NL63
mkdir nl63
cd nl63
NL63=NC_005831.2
# Get genbank genome and convert to FASTA and GFF, then index genome
efetch -db nucleotide -format gb -id $NL63 > $NL63.gb
cat $NL63.gb | seqret -filter -feature -osformat fasta -ofname2 $NL63.gff > $NL63.fa
rm -R *.gb
samtools faidx $NL63.fa
cd ..
cd ..

# Multiple Genome Alignment
cat refs/covid/*.fa refs/sars/*.fa refs/mers/*.fa refs/hku1/*.fa refs/229e/*.fa refs/oc43/*.fa refs/nl63/*.fa > combined_genomes.fa
clustalo -i combined_genomes.fa -o clustalo_msa.fa -v

#Double Genome Alignment
stretcher refs/covid/*.fa -bsequence refs/sars/*.fa -outfile covid_vs_sars.txt
stretcher refs/covid/*.fa -bsequence refs/mers/*.fa -outfile covid_vs_mers.txt
stretcher refs/covid/*.fa -bsequence refs/hku1/*.fa -outfile covid_vs_hku1.txt
stretcher refs/covid/*.fa -bsequence refs/229e/*.fa -outfile covid_vs_229e.txt
stretcher refs/covid/*.fa -bsequence refs/oc43/*.fa -outfile covid_vs_oc43.txt
stretcher refs/covid/*.fa -bsequence refs/nl63/*.fa -outfile covid_vs_nl63.txt

#Align select SARS-Cov-2 Proteins to SARS Genome to find alignment statistics
water receptor_binding_domain.txt refs/sars/*.fa -gapopen 10.0 -gapextend 0.5 -outfile rbd_sars.txt
water spike_protein_S1.txt refs/sars/*.fa -gapopen 10.0 -gapextend 0.5 -outfile S1_sars.txt
water spike_protein_S2.txt refs/sars/*.fa -gapopen 10.0 -gapextend 0.5 -outfile S2_sars.txt
water RNA_binding.txt refs/sars/*.fa -gapopen 10.0 -gapextend 0.5 -outfile RNA_binding_sars.txt
water ACE2_binding.txt refs/sars/*.fa -gapopen 10.0 -gapextend 0.5 -outfile ACE2_binding_sars.txt
water N_protein.txt refs/sars/*.fa -gapopen 10.0 -gapextend 0.5 -outfile N_protein_sars.txt
water pp1ab.txt refs/sars/*.fa -gapopen 10.0 -gapextend 0.5 -outfile pp1ab_sars.txt
water nsp2.txt refs/sars/*.fa -gapopen 10.0 -gapextend 0.5 -outfile nsp2_sars.txt
water nsp3.txt refs/sars/*.fa -gapopen 10.0 -gapextend 0.5 -outfile nsp3_sars.txt
water pp1a.txt refs/sars/*.fa -gapopen 10.0 -gapextend 0.5 -outfile pp1a_sars.txt
water M_protein.txt refs/sars/*.fa -gapopen 10.0 -gapextend 0.5 -outfile M_protein_sars.txt
water E_protein.txt refs/sars/*.fa -gapopen 10.0 -gapextend 0.5 -outfile E_protein_sars.txt
