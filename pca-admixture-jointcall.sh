# pca-admixture feb 2022 (joint call)

sbatch extract_GSA_vars.sh

for i in {1..22}; do /storage/atkinson/home/u242335/bin/bcftools index GSA_gnomad.genomes.v3.1.2.genotypes.chr${i}.vcf.gz; done

/storage/atkinson/home/u242335/bin/bcftools concat GSA_gnomad.genomes.v3.1.2.genotypes.chr1.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr2.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr3.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr4.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr5.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr6.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr7.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr8.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr9.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr10.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr11.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr12.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr13.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr14.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr15.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr16.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr17.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr18.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr19.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr20.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr21.vcf.gz GSA_gnomad.genomes.v3.1.2.genotypes.chr22.vcf.gz --threads 10 -Oz -o GSA_gnomadjoint_allchr_1kg_hgdp_hg38.vcf.gz

/storage/atkinson/home/u242335/bin/bcftools index GSA_gnomadjoint_allchr_1kg_hgdp_hg38.vcf.gz

/storage/atkinson/shared_resources/software/plink2 --vcf GSA_gnomadjoint_allchr_1kg_hgdp_hg38.vcf.gz --make-bed --out GSA_gnomadjoint_allchr_1kg_hgdp_hg38

# limpando vars multialelicas
/storage/atkinson/shared_resources/software/plink2/plink2 --bfile GSA_gnomadjoint_allchr_1kg_hgdp_hg38 --max-alleles 2 --snps-only just-acgt --make-bed --out GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO

./plink --bfile INPD_kids_hg38 --bmerge GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO.bed GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO.bim GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO.fam --make-bed --out merge_teste
./plink --bfile INPD_kids_hg38 --exclude merge_teste-merge.missnp --make-bed --out INPD_kids_hg38_tmp
./plink --bfile GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO --exclude merge_teste-merge.missnp --make-bed --out GSA_1kg_hgdp_short_merge_hg38_tmp
./plink --bfile INPD_kids_hg38_tmp --bmerge GSA_1kg_hgdp_short_merge_hg38_tmp --make-bed --out merged_INPD_1KG_HGDP

rm GSA_1kg_hgdp_short_merge_hg38_tmp.*
rm INPD_kids_hg38_tmp.*


# Variant QC - look at removed geno/maf vars and see if we're not losing anything population specific?
#383038 variants loaded from merged_INPD_1KG_HGDP.bim.
#--geno: 10708 variants removed due to missing genotype data.
#--hwe: 91347 variants removed due to Hardy-Weinberg exact test (founders only).
#28841 variants removed due to allele frequency threshold(s)
#(--maf/--max-maf/--mac/--max-mac).
#252142 variants remaining after main filters.

## Feb 10 - decided not to use hwe flag because it would exclude variants based on population structure.

/storage/atkinson/shared_resources/software/plink2/plink2 --bfile merged_INPD_1KG_HGDP --geno 0.1 --maf 0.01 --make-bed --out merged_INPD_1KG_HGDP_hg38_QC
# run --mind 0.1 - passed and removed file
/storage/atkinson/shared_resources/software/plink2/plink2 --bfile merged_INPD_1KG_HGDP_hg38_QC --mind 0.1 --make-bed --out merged_INPD_1KG_HGDP_hg38_QCmind
rm *_QCmind.*

# final merge with 6342 samples and 343472 variants

###### update sex and FID = Population
## R code for making update file, based on .tsv metadata file.

# old FID old IID
awk '{print $1, $2, $5}' merged_INPD_1KG_HGDP_hg38_QC.fam > fam_info_old

# #FID IID SEX - make update sex, did not do yet

# remove CHMI_CHMI3_WGS2 because it didnt belong to any pop
nano remove_CHMI_CHMI3_WGS2.txt 
/storage/atkinson/shared_resources/software/plink2/plink2 --bfile merged_INPD_1KG_HGDP_hg38_QC --remove remove_CHMI_CHMI3_WGS2.txt --make-bed --out merged_INPD_1KG_HGDP_hg38_QC_removeCHMI

#updating fid
#edit file to remove inpd ids with underscores
nano update_fid_merge_inpd_1kg_hgdp.txt
/storage/atkinson/shared_resources/software/plink2/plink2 --bfile merged_INPD_1KG_HGDP_hg38_QC_removeCHMI --update-ids update_fid_merge_inpd_1kg_hgdp.txt --make-bed --out merged_INPD_1KG_HGDP_POP

###### PCA

/storage/atkinson/shared_resources/software/plink2/plink2 --bfile merged_INPD_1KG_HGDP_POP --pca 5 --out inpd_1kg_hgdp

##### ADMIXTURE: Submit batch job run_admixture.sh
# running ADMIXTURE 
for k in {2..10}; do /storage/atkinson/shared_resources/software/admixture_linux-1.3.0/admixture --cv merged_INPD_1KG_HGDP_POP.bed $k -j6 | tee merged_INPD_1KG_HGDP_$k.log; done

### after run
#check cross valídation errors for best K. plot the values and choose a K in the 'elbow'
grep -h CV log*.out

# plotting - pong
tbl=read.table("hapmap3.3.Q")
barplot(t(as.matrix(tbl)), col=rainbow(3), xlab="Individual #", ylab="Ancestry", border=NA)

pong -m pong_filemap -n pop_order_expandednames.txt -i ind2pop.txt


### SELECTING ONLY POPS OF INTEREST
grep BRA_RS merged_INPD_1KG_HGDP_POP.fam > pops_interest.txt
grep BRA_SP merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep ACB merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep ASW merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep ESN merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep YRI merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Yoruba merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep BantuKenya merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep BantuSouthAfrica merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep MbutiPygmy merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep CLM merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Colombian merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Surui merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Karitiana merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep PEL merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep PUR merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep MXL merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Maya merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Pima merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Japanese merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep JPT merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep CEU merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep GBR merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep IBS merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep TSI merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Italian merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Sardinian merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Tuscan merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep French merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Basque merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt

awk '{print $1, $2}' pops_interest.txt > keep_pops_interest.txt

/storage/atkinson/shared_resources/software/plink2/plink2 --bfile merged_INPD_1KG_HGDP_POP --keep keep_pops_interest.txt --make-bed --out merged_INPD_1KG_HGDP_pops_interest

###### PCA

/storage/atkinson/shared_resources/software/plink2/plink2 --bfile merged_INPD_1KG_HGDP_pops_interest --pca 5 --out inpd_1kg_hgdp_interest

##### ADMIXTURE: Submit batch job run_admixture.sh
# running ADMIXTURE = sbatch scripts. 2 iterations for each k.


### after run
#check cross valídation errors for best K. plot the values and choose a K in the 'elbow'
grep -h CV log*.out

# plotting - pong
tbl=read.table("hapmap3.3.Q")
barplot(t(as.matrix(tbl)), col=rainbow(3), xlab="Individual #", ylab="Ancestry", border=NA)

pong -m pops_interest.filemap -i pops_interest.ind2pop -n pop_order_expandednames.tsv

pops_all.ind2pop

### Select only pops of interest x2 = agora vai
grep BRA_RS merged_INPD_1KG_HGDP_POP.fam > pops_interest.txt
grep BRA_SP merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep ESN merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep YRI merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Yoruba merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep BantuKenya merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep BantuSouthAfrica merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Surui merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Karitiana merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Maya merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Pima merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep CEU merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep GBR merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep IBS merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep TSI merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Italian merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Sardinian merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep Tuscan merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt
grep French merged_INPD_1KG_HGDP_POP.fam >> pops_interest.txt

awk '{print $1, $2}' pops_interest.txt > keep_pops_interest.txt

/storage/atkinson/shared_resources/software/plink2/plink2 --bfile merged_INPD_1KG_HGDP_POP --keep keep_pops_interest.txt --make-bed --out merged_INPD_1KG_HGDP_pops_interest_two

###### PCA

/storage/atkinson/shared_resources/software/plink2/plink2 --bfile merged_INPD_1KG_HGDP_pops_interest_two --pca 5 --out inpd_1kg_hgdp_interest_two

mkdir run1
mkdir run2

##### ADMIXTURE: Submit batch job run_admixture.sh
# running ADMIXTURE = sbatch scripts. 2 iterations for each k.


### after run
#check cross valídation errors for best K. plot the values and choose a K in the 'elbow'
grep -h CV log*.out

#IN R:
read.table("cv.txt", h=F) -> cv_run1
c(2,3,4,5,6,7,8,9) -> cv_run1$V5
cv_run1$V4 <- as.numeric(cv_run1$V4)
ggplot(cv_run1, aes(V5, V4)) + geom_line()

# plotting - pong
tbl=read.table("hapmap3.3.Q")
barplot(t(as.matrix(tbl)), col=rainbow(3), xlab="Individual #", ylab="Ancestry", border=NA)

pong -m run1.filemap -i run1.ind2pop -n pop_order_expandednames.txt

pops_all.ind2pop