# separate ref and target from phasing output. #running from /storage/atkinson/home/u242335/lai

awk '{print $1, $2}' /storage/atkinson/home/u242335/phasing/merged_INPD_1KG_HGDP_pops_interest_second.fam > /storage/atkinson/home/u242335/lai/BHRC_KIDS
head -n 2191 BHRC_KIDS | awk '{print $2}' > BHRC_probands

egrep "Surui|Karitiana|Maya|Pima|Colombian" BHRC_KIDS | awk '{print $2}' > NAT_HGDP_REF


grep IBS ./simu-jointcall/1000GP_Phase3/1000GP_Phase3.sample | awk '{print $1}' > IBS_1KG_REF

grep YRI ./simu-jointcall/1000GP_Phase3/1000GP_Phase3.sample | awk '{print $1}' > YRI_1KG_REF
#manually edited YRI_1KG_REF NA18874 to NA18874A (as it is on the jointcall)

cat NAT_HGDP_REF IBS_1KG_REF YRI_1KG_REF > REF_NAT_EUR_1kg_AFR_1kg.ref

cp BHRC_probands BHRC.notref


######## 


for i in {1..22}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps ../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind BHRC_probands \
--output-haps BHRC_chr${i}.haps BHRC_chr${i}.sample ; done


for i in {1..22}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps ../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind NAT_HGDP_REF \
--output-haps REF_NAT_HGDP_chr${i}.haps REF_NAT_HGDP_chr${i}.sample ; done

for i in {1..22}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps ../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind IBS_1KG_REF \
--output-haps REF_IBS_1KG_chr${i}.haps REF_IBS_1KG_chr${i}.sample ; done

for i in {1..22}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps ../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind YRI_1KG_REF \
--output-haps REF_YRI_1KG_chr${i}.haps REF_YRI_1KG_chr${i}.sample ; done

# shapeit output to rfmix
for i in {1..22}; do \
python ./ancestry_pipeline-master/shapeit2rfmix.py \
--shapeit_hap_ref ./REF_NAT_HGDP_chr${i}.haps,./REF_IBS_1KG_chr${i}.haps,./REF_YRI_1KG_chr${i}.haps \
--shapeit_hap_admixed ./BHRC_chr${i}.haps \
--shapeit_sample_ref ./REF_NAT_HGDP_chr${i}.sample,./REF_IBS_1KG_chr${i}.sample,./REF_YRI_1KG_chr${i}.sample \
--shapeit_sample_admixed ./BHRC_chr${i}.sample \
--chr ${i} \
--genetic_map ../phasing/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt \
--ref_keep ./REF_NAT_EUR_1kg_AFR_1kg.ref \
--admixed_keep ./BHRC.notref \
--out BHRC_RefHGDPNAT_1kgIBS_1kgYRI; done

# running RFMIX
cd /storage/atkinson/home/u242335/lai/RFMix_v1.5.4/

python RunRFMix.py \
-e 2 \
-w 0.2 \
-n 5 \
-G 12 \
--num-threads 16 \
--use-reference-panels-in-EM \
--forward-backward \
PopPhased \
/storage/atkinson/home/u242335/lai/BHRC_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.alleles \
/storage/atkinson/home/u242335/lai/BHRC_RefHGDPNAT_1kgIBS_1kgYRI.classes \
/storage/atkinson/home/u242335/lai/BHRC_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.snp_locations \
-o /storage/atkinson/home/u242335/lai/BHRC_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix


#### split to see if it works/ reduces runtime #running from /storage/atkinson/home/u242335/lai/bhrcsplit
head -n 1096 ../BHRC_probands > BHRC_probands1
tail -n 1095 ../BHRC_probands > BHRC_probands2

cp BHRC_probands1 BHRC_1.notref
cp BHRC_probands2 BHRC_2.notref

for i in {1..22}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind BHRC_probands1 \
--output-haps BHRC1_chr${i}.haps BHRC1_chr${i}.sample ; done

for i in {1..22}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind BHRC_probands2 \
--output-haps BHRC2_chr${i}.haps BHRC2_chr${i}.sample ; done


# shapeit output to rfmix
for i in {1..22}; do \
python ../ancestry_pipeline-master/shapeit2rfmix.py \
--shapeit_hap_ref ../REF_NAT_HGDP_chr${i}.haps,../REF_IBS_1KG_chr${i}.haps,../REF_YRI_1KG_chr${i}.haps \
--shapeit_hap_admixed ./BHRC1_chr${i}.haps \
--shapeit_sample_ref ../REF_NAT_HGDP_chr${i}.sample,../REF_IBS_1KG_chr${i}.sample,../REF_YRI_1KG_chr${i}.sample \
--shapeit_sample_admixed ./BHRC1_chr${i}.sample \
--chr ${i} \
--genetic_map ../../phasing/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt \
--ref_keep ../REF_NAT_EUR_1kg_AFR_1kg.ref \
--admixed_keep ./BHRC_1.notref \
--out BHRC1_RefHGDPNAT_1kgIBS_1kgYRI; done

for i in {1..22}; do \
python ../ancestry_pipeline-master/shapeit2rfmix.py \
--shapeit_hap_ref ../REF_NAT_HGDP_chr${i}.haps,../REF_IBS_1KG_chr${i}.haps,../REF_YRI_1KG_chr${i}.haps \
--shapeit_hap_admixed ./BHRC2_chr${i}.haps \
--shapeit_sample_ref ../REF_NAT_HGDP_chr${i}.sample,../REF_IBS_1KG_chr${i}.sample,../REF_YRI_1KG_chr${i}.sample \
--shapeit_sample_admixed ./BHRC2_chr${i}.sample \
--chr ${i} \
--genetic_map ../../phasing/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt \
--ref_keep ../REF_NAT_EUR_1kg_AFR_1kg.ref \
--admixed_keep ./BHRC_2.notref \
--out BHRC2_RefHGDPNAT_1kgIBS_1kgYRI; done


# running RFMIX
cd /storage/atkinson/home/u242335/lai/RFMix_v1.5.4/

python RunRFMix.py \
-e 2 \
-w 0.2 \
-n 5 \
-G 12 \
--use-reference-panels-in-EM \
--forward-backward \
PopPhased \
/storage/atkinson/home/u242335/lai/bhrcsplit/BHRC1_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.alleles \
/storage/atkinson/home/u242335/lai/bhrcsplit/BHRC1_RefHGDPNAT_1kgIBS_1kgYRI.classes \
/storage/atkinson/home/u242335/lai/bhrcsplit/BHRC1_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.snp_locations \
-o /storage/atkinson/home/u242335/lai/bhrcsplit/BHRC1_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix


cd /storage/atkinson/home/u242335/lai/RFMix_v1.5.4/

python RunRFMix.py \
-e 2 \
-w 0.2 \
-n 5 \
-G 12 \
--use-reference-panels-in-EM \
--forward-backward \
PopPhased \
/storage/atkinson/home/u242335/lai/bhrcsplit/BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.alleles \
/storage/atkinson/home/u242335/lai/bhrcsplit/BHRC2_RefHGDPNAT_1kgIBS_1kgYRI.classes \
/storage/atkinson/home/u242335/lai/bhrcsplit/BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.snp_locations \
-o /storage/atkinson/home/u242335/lai/bhrcsplit/BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix

### JOINING TWO HALVES ### ARRUMAR AQUI, O \0 está juntando as colunas


#### Join viterbi
#get only indiv cols from BHRC2 Viterbi and FB. set delimiter. substitute 2 for the column where indiv starts (555, because ref has 554 haplotypes)
#join
for i in {20..22}; do \
cut -d ' ' -f 555- BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt > INDIV_BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt; \
paste -d' ' BHRC1_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt INDIV_BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt > \
BHRC_JOIN_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.2.Viterbi.txt ; done

# CHECK - count columns (should be 2744) 2190 BHRC2 Haplotypes. 554 ref haplotypes.
#should be 2190
# should be 4936 (ref haplotypes + bhrc haplotypes)

for i in {20..22}; do \
awk '{print NF}' BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt  | sort -nu | tail -n 1 ; \
awk '{print NF}' INDIV_BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt  | sort -nu | tail -n 1 ; \
awk '{print NF}' BHRC_JOIN_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.2.Viterbi.txt  | sort -nu | tail -n 1 ; done


#### Join ForwardBackward

for i in {20..22}; do \
cut -d ' ' -f 1663- BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt > INDIV_BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt; \
paste -d' ' BHRC1_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt INDIV_BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt > \
BHRC_JOIN_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.2.ForwardBackward.txt ; done

# CHECK - count columns (should be 8232 = 2744 x 3) 2190 BHRC2 Haplotypes. 554 ref haplotypes. 3 ancestries.
#should be 2190 x 3 = 6570

# should be 4936 x 3 = 14808
for i in {20..22}; do \
awk '{print NF}' BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt | sort -nu | tail -n 1 ; \
awk '{print NF}' INDIV_BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt | sort -nu | tail -n 1 ; \
awk '{print NF}' BHRC_JOIN_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.2.ForwardBackward.txt  | sort -nu | tail -n 1 ; done


#check if snp locations file is the same between both -same so just copy

for i in {20..22}; do \
cp BHRC1_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.snp_locations BHRC_JOIN_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.snp_locations; done

## processing viterbi (putting on chr pos)

for i in {20..22}; do \
sed 's/1/0/g' BHRC_JOIN_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.2.Viterbi.txt | sed 's/2/1/g' | sed 's/3/2/g' > BHRC_JOIN_${i}_calls; \
awk '{print $1, $3}' BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.map | sed 's/:.*//g' | awk '{print $2, $1}' > BHRC_chrpos_${i}; \
paste BHRC_chrpos_${i} BHRC_JOIN_${i}_calls | sed 's/\t/ /g' > BHRC_RFMIX_${i}.Viterbi.txt; done


####################
#### split 3 times for chr 10,9,8, 6-1

#### split to see if it works/ reduces runtime #running from /storage/atkinson/home/u242335/lai/bhrcsplit
head -n 731 ../BHRC_probands > BHRC_1st
head -n 1461 ../BHRC_probands | tail -n 730 > BHRC_2nd
tail -n 730 ../BHRC_probands > BHRC_3rd

cp BHRC_1st BHRC_1st.notref
cp BHRC_2nd BHRC_2nd.notref
cp BHRC_3rd BHRC_3rd.notref

for i in {10..1}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind BHRC_1st \
--output-haps BHRC1st_chr${i}.haps BHRC1st_chr${i}.sample ; done

for i in {10..1}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind BHRC_2nd \
--output-haps BHRC2nd_chr${i}.haps BHRC2nd_chr${i}.sample ; done

for i in {10..1}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz ../../phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind BHRC_3rd \
--output-haps BHRC3rd_chr${i}.haps BHRC3rd_chr${i}.sample ; done

# shapeit output to rfmix
for i in {10..1}; do \
python ../ancestry_pipeline-master/shapeit2rfmix.py \
--shapeit_hap_ref ../REF_NAT_HGDP_chr${i}.haps,../REF_IBS_1KG_chr${i}.haps,../REF_YRI_1KG_chr${i}.haps \
--shapeit_hap_admixed ./BHRC1st_chr${i}.haps \
--shapeit_sample_ref ../REF_NAT_HGDP_chr${i}.sample,../REF_IBS_1KG_chr${i}.sample,../REF_YRI_1KG_chr${i}.sample \
--shapeit_sample_admixed ./BHRC1st_chr${i}.sample \
--chr ${i} \
--genetic_map ../../phasing/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt \
--ref_keep ../REF_NAT_EUR_1kg_AFR_1kg.ref \
--admixed_keep ./BHRC_1st.notref \
--out BHRC1st_RefHGDPNAT_1kgIBS_1kgYRI; done

for i in {10..1}; do \
python ../ancestry_pipeline-master/shapeit2rfmix.py \
--shapeit_hap_ref ../REF_NAT_HGDP_chr${i}.haps,../REF_IBS_1KG_chr${i}.haps,../REF_YRI_1KG_chr${i}.haps \
--shapeit_hap_admixed ./BHRC2nd_chr${i}.haps \
--shapeit_sample_ref ../REF_NAT_HGDP_chr${i}.sample,../REF_IBS_1KG_chr${i}.sample,../REF_YRI_1KG_chr${i}.sample \
--shapeit_sample_admixed ./BHRC2nd_chr${i}.sample \
--chr ${i} \
--genetic_map ../../phasing/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt \
--ref_keep ../REF_NAT_EUR_1kg_AFR_1kg.ref \
--admixed_keep ./BHRC_2nd.notref \
--out BHRC2nd_RefHGDPNAT_1kgIBS_1kgYRI; done

for i in {10..1}; do \
python ../ancestry_pipeline-master/shapeit2rfmix.py \
--shapeit_hap_ref ../REF_NAT_HGDP_chr${i}.haps,../REF_IBS_1KG_chr${i}.haps,../REF_YRI_1KG_chr${i}.haps \
--shapeit_hap_admixed ./BHRC3rd_chr${i}.haps \
--shapeit_sample_ref ../REF_NAT_HGDP_chr${i}.sample,../REF_IBS_1KG_chr${i}.sample,../REF_YRI_1KG_chr${i}.sample \
--shapeit_sample_admixed ./BHRC3rd_chr${i}.sample \
--chr ${i} \
--genetic_map ../../phasing/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt \
--ref_keep ../REF_NAT_EUR_1kg_AFR_1kg.ref \
--admixed_keep ./BHRC_3rd.notref \
--out BHRC3rd_RefHGDPNAT_1kgIBS_1kgYRI; done


### JOINING THREE PARTS


#### Join viterbi
#get only indiv cols from BHRC2nd and BHRC3rd Viterbi and FB. set delimiter. substitute 2 for the column where indiv starts (555, because ref has 554 haplotypes)
#join
for i in {6..1}; do \
cut -d ' ' -f 555- BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt > INDIV_BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt; \
cut -d ' ' -f 555- BHRC3_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt > INDIV_BHRC3_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt; \
paste -d' ' BHRC1_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt INDIV_BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt INDIV_BHRC3_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt > \
BHRC_JOIN_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.2.Viterbi.txt ; done

# CHECK - count columns (should be 2014) 1460 BHRC2 Haplotypes. 554 ref haplotypes.
#should be 1460
#should be 1460
# should be 4936 (ref haplotypes + bhrc haplotypes)

for i in {6..1}; do \
awk '{print NF}' BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt  | sort -nu | tail -n 1 ; \
awk '{print NF}' BHRC3_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt  | sort -nu | tail -n 1 ; \
awk '{print NF}' INDIV_BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt  | sort -nu | tail -n 1 ; \
awk '{print NF}' INDIV_BHRC3_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.Viterbi.txt  | sort -nu | tail -n 1 ; \
awk '{print NF}' BHRC_JOIN_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.2.Viterbi.txt  | sort -nu | tail -n 1 ; done


#### Join ForwardBackward

for i in {6..1}; do \
cut -d ' ' -f 1663- BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt > INDIV_BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt; \
cut -d ' ' -f 1663- BHRC3_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt > INDIV_BHRC3_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt; \
paste -d' ' BHRC1_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt INDIV_BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt INDIV_BHRC3_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt > \
BHRC_JOIN_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.2.ForwardBackward.txt ; done

# CHECK - count columns (should be 6042 = 2014 x 3) 1460 BHRC2 Haplotypes. 554 ref haplotypes. 3 ancestries.
#should be 1460 x 3 = 4380

# should be 4936 x 3 = 14808
for i in {6..1}; do \
awk '{print NF}' BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt | sort -nu | tail -n 1 ; \
awk '{print NF}' BHRC3_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt | sort -nu | tail -n 1 ; \
awk '{print NF}' INDIV_BHRC2_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt | sort -nu | tail -n 1 ; \
awk '{print NF}' INDIV_BHRC3_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.rfmix.2.ForwardBackward.txt | sort -nu | tail -n 1 ; \
awk '{print NF}' BHRC_JOIN_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.2.ForwardBackward.txt  | sort -nu | tail -n 1 ; done


#check if snp locations file is the same between both -same so just copy

for i in {6..1}; do \
cp BHRC1st_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.snp_locations BHRC_JOIN_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.snp_locations; done

## processing viterbi (putting on chr pos)

for i in {6..1}; do \
sed 's/1/0/g' BHRC_JOIN_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.2.Viterbi.txt | sed 's/2/1/g' | sed 's/3/2/g' > BHRC_JOIN_${i}_calls; \
awk '{print $1, $3}' BHRC2nd_RefHGDPNAT_1kgIBS_1kgYRI_chr${i}.map | sed 's/:.*//g' | awk '{print $2, $1}' > BHRC_chrpos_${i}; \
paste BHRC_chrpos_${i} BHRC_JOIN_${i}_calls | sed 's/\t/ /g' > BHRC_RFMIX_${i}.Viterbi.txt; done











###############
# collapse RFMix into bed files. from #running from /storage/atkinson/home/u242335/lai/bedcollapse

while read p; do \
python ../ancestry_pipeline-master/collapse_ancestry.py \
--rfmix ../bhrcsplit/BHRC_JOIN_RefHGDPNAT_1kgIBS_1kgYRI_chr1.2.Viterbi.txt  \
--snp_locations ../bhrcsplit/BHRC_JOIN_RefHGDPNAT_1kgIBS_1kgYRI_chr1.snp_locations   \
--fbk ../bhrcsplit/BHRC_JOIN_RefHGDPNAT_1kgIBS_1kgYRI_chr1.2.ForwardBackward.txt \
--fbk_threshold 0.9 \
--ind "$p" \
--ind_info ../BHRC_RefHGDPNAT_1kgIBS_1kgYRI.sample \
--pop_labels NAT,EUR,AFR \
--out ./"$p"
done < ../BHRC.notref

# plot karyograms
# get centromere files

IND='HG02481'; python plot_karyogram.py \
--bed_a ${IND}_A.bed \
--bed_b ${IND}_B.bed \
--ind ${IND} \
--out ${IND}.png

# estimate global ancestry from LAI (to compare with admixture)
for POP in ACB ASW CLM MXL PEL PUR; do python lai_global.py \
--bed_list bed_list_${POP}.txt \
--ind_list ${POP}.inds \
--pops AFR,EUR,NAT \
--out lai_global_${POP}.txt; done

############################################

### TESTING RFMIX V2 With 100 people ###

############################################

# from lai/rfmix2_test
# making reference VCF
for i in {1..22}; do \
/storage/atkinson/home/u242335/phasing/shapeit -convert --input-haps /storage/atkinson/home/u242335/phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.haps.gz /storage/atkinson/home/u242335/phasing/BHRC_1KG_HGDP_tailored_chrpos_hg38_chr${i}.sample \
--include-ind /storage/atkinson/home/u242335/lai/REF_NAT_EUR_1kg_AFR_1kg.ref \
--output-vcf REF_NAT_HGDP_IBS_1KG_YRI_1KG_chr${i}.vcf; bgzip -c REF_NAT_HGDP_IBS_1KG_YRI_1KG_chr${i}.vcf > REF_NAT_HGDP_IBS_1KG_YRI_1KG_chr${i}.vcf.gz ; done

for i in {1..22}; do /storage/atkinson/home/u242335/bcftools-1.14/bcftools index /storage/atkinson/home/u242335/lai/rfmix2_test/REF_NAT_HGDP_IBS_1KG_YRI_1KG_chr${i}.vcf.gz ; done

# selecting first 100 indivs from BHRC vcf
head -n 100 /storage/atkinson/home/u242335/lai/BHRC.notref > /storage/atkinson/home/u242335/lai/rfmix2_test/BHRC_first100
for i in {1..22}; do /storage/atkinson/home/u242335/bcftools-1.14/bcftools index /storage/atkinson/home/u242335/lai/tractor/BHRC_chr${i}.vcf.gz ; done
for i in {1..22}; do /storage/atkinson/home/u242335/bcftools-1.14/bcftools view /storage/atkinson/home/u242335/lai/tractor/BHRC_chr${i}.vcf.gz -S BHRC_first100 -Oz -o /storage/atkinson/home/u242335/lai/rfmix2_test/BHRC_100_chr${i}.vcf.gz ; done

# making tsv file

printf 'AMR\n%.0s' {1..62} > REF_RFMIX2.ind2pop
printf 'IBS\n%.0s' {1..107} >> REF_RFMIX2.ind2pop 
printf 'YRI\n%.0s' {1..108} >> REF_RFMIX2.ind2pop 

paste REF_NAT_EUR_1kg_AFR_1kg.ref REF_RFMIX2.ind2pop > REF_NAT_EUR_1kg_AFR_1kg_forRFMix2.tsv

#making genetic map acceptable for rfmix 2
for i in {1..22}; do \
awk -v i=${i} '{print i, $1, $3}' /storage/atkinson/home/u242335/phasing/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt | sed 's/ /\t/g' | sed '1 s/./#&/' >  /storage/atkinson/home/u242335/phasing/genetic_map/rfmix2/rfmix2_genetic_map_chr${i}_hg38_firstLine.txt ; done

#running rfmix2
## RUNTIME : 2 hours for all 22 chr (GSA Vars)

for i in {22..1}; do \
/storage/atkinson/shared_resources/software/rfmix/rfmix \
-f /storage/atkinson/home/u242335/lai/rfmix2_test/BHRC_100_chr${i}.vcf.gz  \
-r /storage/atkinson/home/u242335/lai/rfmix2_test/REF_NAT_HGDP_IBS_1KG_YRI_1KG_chr${i}.vcf.gz \
-n 5 \
-G 12 \
--reanalyze-reference -e 2 \
-m /storage/atkinson/home/u242335/lai/rfmix2_test/REF_NAT_EUR_1kg_AFR_1kg_forRFMix2.tsv \
-g /storage/atkinson/home/u242335/phasing/genetic_map/rfmix2/rfmix2_genetic_map_chr${i}_hg38_firstLine.txt \
-o /storage/atkinson/home/u242335/lai/rfmix2_test/BHRC_100_chr${i}.deconvoluted \
--chromosome=${i} ; done


## Run v2 for everyone?? Array job.

/storage/atkinson/shared_resources/software/rfmix/rfmix \
-f /storage/atkinson/home/u242335/lai/tractor/BHRC_chr$SLURM_ARRAY_TASK_ID.vcf.gz \
-r /storage/atkinson/home/u242335/lai/rfmix2_test/REF_NAT_HGDP_IBS_1KG_YRI_1KG_chr$SLURM_ARRAY_TASK_ID.vcf.gz \
-n 5 \
-G 12 \
--reanalyze-reference -e 2 \
-m /storage/atkinson/home/u242335/lai/rfmix2_test/REF_NAT_EUR_1kg_AFR_1kg_forRFMix2.tsv \
-g /storage/atkinson/home/u242335/phasing/genetic_map/rfmix2/rfmix2_genetic_map_chr$SLURM_ARRAY_TASK_ID_hg38_firstLine.txt \
-o /storage/atkinson/home/u242335/lai/rfmix2_test/BHRC_ALL_chr$SLURM_ARRAY_TASK_ID.deconvoluted \
--chromosome=$SLURM_ARRAY_TASK_ID