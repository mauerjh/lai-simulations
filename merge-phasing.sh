#/storage/atkinson/home/u242335/phasing
# Using shapeit2 without reference, with the merged-tailored-dataset of BHRC + Interest pops of 1KG and HGDP


# Merging BHRC and Joint call 1KG-HGDP (all hg38)

##### check merge files for same-position variants: merge by chr-pos #####

## removing '.' IDs from ref with plink2
/storage/atkinson/shared_resources/software/plink2/plink2 --bfile GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO --exclude remove_dot_vars.txt --make-bed --out GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO_semponto

# changing RSID to chr:pos:a1:a2
awk 'BEGIN {FS = "\t"};{print $1,$1":"$4":"$5":"$6,$3,$4,$5,$6}' INPD2020.bim> INPD2020_chrpos.bim
cp INPD2020.bed INPD2020_chrpos.bed
cp INPD2020.fam INPD2020_chrpos.fam

awk 'BEGIN {FS = "\t"};{print $1,$1":"$4":"$5":"$6,$3,$4,$5,$6}' GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO_semponto.bim> GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO_semponto_chrpos.bim
cp GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO_semponto.bed GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO_semponto_chrpos.bed
cp GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO_semponto.fam GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO_semponto_chrpos.fam

#remove duplicates from ref
# 1- list duplicates = no dup on jointcall, 114 dups on bhrc
/storage/atkinson/shared_resources/software/plink2/plink2 --bfile GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO_semponto_chrpos --rm-dup list --out jointcall_dup

/storage/atkinson/shared_resources/software/plink2/plink2 --bfile INPD2020_chrpos --rm-dup list --out bhrc_dup

# 2 - remove dup ids from bhrc. --rm-dup: 122 duplicated IDs, 237 variants removed. 
/storage/atkinson/shared_resources/software/plink2/plink2 --bfile INPD2020_chrpos --rm-dup exclude-mismatch --make-bed --out INPD2020_chrpos_rmdup

#intersection of vars - R 
# extract isec sites = 384655 SNPs.
/storage/atkinson/shared_resources/software/plink2/plink2 --bfile GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO_semponto_chrpos --extract isecSNPs_bhrc_jointcall.txt --make-bed --out GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO_semponto_chrpos_isecSNP

/storage/atkinson/shared_resources/software/plink2/plink2 --bfile INPD2020_chrpos_rmdup --extract isecSNPs_bhrc_jointcall.txt --make-bed --out INPD2020_chrpos_rmdup_isecSNP


# merge
./plink --bfile INPD2020_chrpos_rmdup_isecSNP --bmerge GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO_semponto_chrpos_isecSNP.bed GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO_semponto_chrpos_isecSNP.bim GSA_gnomadjoint_allchr_1kg_hgdp_hg38_LIMPO_semponto_chrpos_isecSNP.fam --make-bed --out merged_INPDall_1KG_HGDP_chrpos_semponto_isecSNP

# Variant QC
/storage/atkinson/shared_resources/software/plink2/plink2 --bfile merged_INPDall_1KG_HGDP_chrpos_semponto_isecSNP --geno 0.1 --maf 0.01 --make-bed --out merged_INPDall_1KG_HGDP_chrpos_semponto_isecSNP_QC
#--geno: 2961 variants removed due to missing genotype data.
#24674 variants removed due to allele frequency threshold(s)--maf

# indiv qc - didnt remove anyone so removed the out file
/storage/atkinson/shared_resources/software/plink2/plink2 --bfile merged_INPDall_1KG_HGDP_chrpos_semponto_isecSNP_QC --mind 0.1 --make-bed --out merged_INPDall_1KG_HGDP_chrpos_semponto_isecSNP_QCmind
rm *_QCmind.*

# remove indiv CHMI_CHMI3_WGS2 because it didnt belong to any pop
/storage/atkinson/shared_resources/software/plink2/plink2 --bfile merged_INPDall_1KG_HGDP_chrpos_semponto_isecSNP_QC --remove ../pca-admixture/remove_CHMI_CHMI3_WGS2.txt --make-bed --out merged_INPDall_1KG_HGDP_hg38_QC_removeCHMI

#edit file to remove bhrc ids with underscores 
nano merged_INPDall_1KG_HGDP_hg38_QC_removeCHMI.fam

# select only pops of interest from the merge file

head -n 2191 merged_INPDall_1KG_HGDP_hg38_QC_removeCHMI.fam | awk '{print $1, $2}' >  inpd-all_1kg_hgdp_keep.txt
awk '{print $1, $2}' merged_INPDall_1KG_HGDP_hg38_QC_removeCHMI.fam | grep 'M1' >>  inpd-all_1kg_hgdp_keep.txt
awk '{print $1, $2}' merged_INPDall_1KG_HGDP_hg38_QC_removeCHMI.fam | grep 'M2' >>  inpd-all_1kg_hgdp_keep.txt
awk '{print $1, $2}' merged_INPDall_1KG_HGDP_hg38_QC_removeCHMI.fam | grep 'P1' >>  inpd-all_1kg_hgdp_keep.txt
awk '{print $1, $2}' merged_INPDall_1KG_HGDP_hg38_QC_removeCHMI.fam | grep 'P2' >>  inpd-all_1kg_hgdp_keep.txt

tail -n 1031 ../pca-admixture/keep_pops_interest.txt | awk '{print 0, $2}' >> inpd-all_1kg_hgdp_keep.txt

/storage/atkinson/shared_resources/software/plink2/plink2 --bfile merged_INPDall_1KG_HGDP_hg38_QC_removeCHMI --keep inpd-all_1kg_hgdp_keep.txt --make-bed --out merged_INPDall_1KG_HGDP_hg38_chrpos_pops_interest

#another var qc
/storage/atkinson/shared_resources/software/plink2/plink2 --bfile merged_INPDall_1KG_HGDP_hg38_chrpos_pops_interest --geno 0.05 --maf 0.01 --make-bed --out merged_INPDall_1KG_HGDP_chrpos_pops_interest_QC

# shapeit2 check - 11 trios with high mendelian error
# dividir por chr
for i in {1..22};
do /storage/atkinson/shared_resources/software/plink2/plink2 \
--bfile merged_INPDall_1KG_HGDP_chrpos_pops_interest_QCmendel \
--chr ${i} \
--make-bed \
--out BHRC_1KG_HGDP_tailored_hg38_chrpos_chr${i};
done

# SHAPEIT2 check
for i in {1..22}; 
do ./shapeit \
-check \
-B BHRC_1KG_HGDP_tailored_hg38_chrpos_chr${i} \
--input-map ./genetic_map/genetic_map_chr${i}_hg38_firstLine.txt \
--output-log BHRC_1KG_HGDP_tailored_hg38_chrpos_chr${i}.mendel; 
done

# check for mendelian errors on plink
./plink --bfile merged_INPDall_1KG_HGDP_chrpos_pops_interest_QC --me 0.05 0.1 --make-bed --out merged_INPDall_1KG_HGDP_chrpos_pops_interest_QCmendel
#83 variants and 33 people excluded.
#352806 variants and 6365 people pass filters and QC.

##### Fixing mendelian errors in BHRC

## Discovering why 33 individuals are problematic.
# mendel report plink. using file with --geno 0.05
./plink --bfile merged_INPDall_1KG_HGDP_chrpos_pops_interest_QC --mendel summaries-only --out merged_mendel_summary

# list of individuals
awk '{if($3>15000)print$2}' <merged_mendel_summary.imendel | awk '{print 0, $1}' > removed_mendel_ids.txt

# calculating ibd
./plink --bfile merged_removedmendel --genome --out ibd_calculation

# some parents are not actually parents. manually updated .fam to adjust this - merged_INPDall_1KG_HGDP_chrpos_pops_interest_QC_updfam.fam

#check mendelian error again - excludes 85 variants but no people.
./plink --bfile merged_INPDall_1KG_HGDP_chrpos_pops_interest_QC_updfam --me 0.05 0.1 --make-bed --out merged_INPDall_1KG_HGDP_chrpos_pops_interest_QCmendel_updfam

#### divide by chr
for i in {1..22};
do /storage/atkinson/shared_resources/software/plink2/plink2 \
--bfile merged_INPDall_1KG_HGDP_chrpos_pops_interest_QCmendel_updfam \
--chr ${i} \
--make-bed \
--out BHRC_1KG_HGDP_tailored_hg38_chrpos_updfam_chr${i};
done

## RUN SHAPEIT2 #made a slurm script for this running 3 chr per job. 

for i in {1..22}; 
do ./shapeit \
-B BHRC_1KG_HGDP_tailored_hg38_chrpos_updfam_chr${i} \
--duohmm \
--input-map ./genetic_map/genetic_map_chr${i}_hg38_firstLine.txt \
--output-max BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz \
BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--thread 10; done
