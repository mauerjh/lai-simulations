#####################
#Prepare the haps files to RFmix (from Alicia Martin):
#############
#23/03/21:  criei a pasta /genetica_1/Marq/teste_ref3way para deixar separado do resto
#23/03/21: Repetir o script com AdmixPOP=  MAM, PEL, exEURadmix, exNATadmix, Baiano (Pardo e Afam não foram simulados ainda)
## Edit here the ancestry orders and which population will you test, remember to manually write the .dat file.
export  anc1=
export  anc2=
export  anc3=
export AdmixPOP=Baiano
# 0/50/50 = Baiano
# 0/70/30 = Pardo
# 0/30/70 = Afam
# 33/33/34 = Mix

cp /genetica_1/Jessica/admixsimu/PEL.noref $AdmixPOP.no
# Jessica cuidado aqui, precisa trocar o baiano
sed 's/PEL/MAM/g' $AdmixPOP.no > $AdmixPOP.noref

###
#30/03/21: two-way ref 2way, prepare samples + rfmix. Rodar da pasta /genetica_1/Marq 
# MAM, PEL, exNAT, exEUR
for i in {1..1}; do python /genetica_1/Jessica/admixsimu/ancestry_pipeline-master/shapeit2rfmix.py --shapeit_hap_ref ./REF_NAT_HGDP_chr${i}_nc.hap,./REF_EUR_1kg_HGDP_chr${i}_nc.hap --shapeit_hap_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.haps --shapeit_sample_ref ./REF_NAT_HGDP.sample,./REF_EUR_1kg.sample --shapeit_sample_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.sample --chr ${i} --genetic_map ../Jessica/admixsimu/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt --ref_keep /genetica_1/Marq/REF_NAT_EUR_1kg --admixed_keep ./$AdmixPOP.noref --out $AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS; done
cd /genetica_1/SGDP/RFMix_v1.5.4/
for i in {1..1}; do nohup python RunRFMix.py -e 2 -w 0.2 -n 5 -G 9 --num-threads 20 --forward-backward --use-reference-panels-in-EM TrioPhased /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_chr1.alleles /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS.classes /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_chr1.snp_locations -o /genetica_1/Marq/$AdmixPOP.RFmixed_RefHGDPNAT_Ref1kgIBS_${i}; done &

#Baiano
for i in {1..1}; do python /genetica_1/Jessica/admixsimu/ancestry_pipeline-master/shapeit2rfmix.py --shapeit_hap_ref REF_AFR_1kg_HGDP_chr${i}_nc.hap,REF_EUR_1kg_HGDP_chr${i}_nc.hap --shapeit_hap_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.haps --shapeit_sample_ref ./REF_AFR_1KG.sample,./REF_EUR_1kg.sample --shapeit_sample_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.sample --chr ${i} --genetic_map ../Jessica/admixsimu/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt --ref_keep ./REF_EUR_1kg_REF_AFR_1KG --admixed_keep ./$AdmixPOP.noref --out $AdmixPOP.admixsimu_RFmixed_Ref1kgYRI_Ref1kgIBS; done

cd /genetica_1/SGDP/RFMix_v1.5.4/
for i in {1..1}; do python RunRFMix.py -e 2 -w 0.2 -n 5 -G 9 --num-threads 20 --forward-backward --use-reference-panels-in-EM TrioPhased /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_Ref1kgYRI_Ref1kgIBS_chr1.alleles /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_Ref1kgYRI_Ref1kgIBS.classes /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_Ref1kgYRI_Ref1kgIBS_chr1.snp_locations -o /genetica_1/Marq/$AdmixPOP.2.RFmixed_Ref1kgYRI_Ref1kgIBS_1; done


###
#23/03/21: Aqui usando como referencia o que era usado no three-way, porém nas simulações 2way

for i in {1..1}; do python /genetica_1/Jessica/admixsimu/ancestry_pipeline-master/shapeit2rfmix.py --shapeit_hap_ref ../REF_PEL_1KG_chr${i}_nc.hap,../REF_EUR_1kg_HGDP_chr${i}_nc.hap,../REF_AFR_1kg_HGDP_chr${i}_nc.hap --shapeit_hap_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.haps --shapeit_sample_ref ../REF_PEL_1KG.sample,../REF_EUR_1kg.sample,../REF_AFR_1KG.sample --shapeit_sample_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.sample --chr ${i} --genetic_map ../../Jessica/admixsimu/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt --ref_keep /genetica_1/Marq/REF_PEL_1KG_REF_EUR_1kg_REF_AFR_1KG --admixed_keep ./$AdmixPOP.noref --out $AdmixPOP.admixsimu_RFmixed_Ref1kgPEL_Ref1kgIBS_Ref1kgYRI; done

for i in {1..1}; do python /genetica_1/Jessica/admixsimu/ancestry_pipeline-master/shapeit2rfmix.py --shapeit_hap_ref ../REF_NAT_HGDP_chr${i}_nc.hap,../REF_EUR_1kg_HGDP_chr${i}_nc.hap,../REF_AFR_1kg_HGDP_chr${i}_nc.hap --shapeit_hap_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.haps --shapeit_sample_ref ../REF_NAT_HGDP.sample,../REF_EUR_1kg.sample,../REF_AFR_1KG.sample --shapeit_sample_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.sample --chr ${i} --genetic_map ../../Jessica/admixsimu/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt --ref_keep /genetica_1/Marq/REF_NAT_EUR_1kg_AFR_1KG --admixed_keep ./$AdmixPOP.noref --out $AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI; done

for i in {1..1}; do python /genetica_1/Jessica/admixsimu/ancestry_pipeline-master/shapeit2rfmix.py --shapeit_hap_ref ../REF_PEL_EAS_1KG_chr${i}_nc.hap,../REF_EUR_1kg_HGDP_chr${i}_nc.hap,../REF_AFR_1kg_HGDP_chr${i}_nc.hap --shapeit_hap_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.haps --shapeit_sample_ref ../REF_PEL_EAS_1KG.sample,../REF_EUR_1kg.sample,../REF_AFR_1KG.sample --shapeit_sample_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.sample --chr ${i} --genetic_map ../../Jessica/admixsimu/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt --ref_keep /genetica_1/Marq/REF_PEL_EAS_1KG_REF_EUR_1kg_AFR_1KG --admixed_keep ./$AdmixPOP.noref --out $AdmixPOP.admixsimu_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_Ref1kgYRI; done

###
# Finally Run RFmix. Tem que rodar da pasta do RFmix. O que a Jessiaca baixou ainda nao esta instalado. Estou rodando no
cd /genetica_1/SGDP/RFMix_v1.5.4/

for i in {1..1}; do nohup python RunRFMix.py -e 2 -w 0.2 -n 5 -G 9 --num-threads 30 --forward-backward --use-reference-panels-in-EM TrioPhased /genetica_1/Marq/teste_ref3way/$AdmixPOP.admixsimu_RFmixed_Ref1kgPEL_Ref1kgIBS_Ref1kgYRI_chr1.alleles /genetica_1/Marq/teste_ref3way/$AdmixPOP.admixsimu_RFmixed_Ref1kgPEL_Ref1kgIBS_Ref1kgYRI.classes /genetica_1/Marq/teste_ref3way/$AdmixPOP.admixsimu_RFmixed_Ref1kgPEL_Ref1kgIBS_Ref1kgYRI_chr1.snp_locations -o /genetica_1/Marq/teste_ref3way/$AdmixPOP.RFmixed_Ref1kgPEL_Ref1kgIBS_Ref1kgYRI${i}; done &

for i in {1..1}; do nohup python RunRFMix.py -e 2 -w 0.2 -n 5 -G 9 --num-threads 15 --forward-backward --use-reference-panels-in-EM TrioPhased /genetica_1/Marq/teste_ref3way/$AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_chr1.alleles /genetica_1/Marq/teste_ref3way/$AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI.classes /genetica_1/Marq/teste_ref3way/$AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_chr1.snp_locations -o /genetica_1/Marq/teste_ref3way/$AdmixPOP.RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_${i}; done &

for i in {1..1}; do nohup python RunRFMix.py -e 2 -w 0.2 -n 5 -G 9 --num-threads 30 --forward-backward --use-reference-panels-in-EM TrioPhased /genetica_1/Marq/teste_ref3way/$AdmixPOP.admixsimu_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_Ref1kgYRI_chr1.alleles /genetica_1/Marq/teste_ref3way/$AdmixPOP.admixsimu_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_Ref1kgYRI.classes /genetica_1/Marq/teste_ref3way/$AdmixPOP.admixsimu_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_Ref1kgYRI_chr1.snp_locations -o /genetica_1/Marq/teste_ref3way/$AdmixPOP.RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_Ref1kgYRI_${i}; done &


# Jessica Parar aqui

### I am here on 3rd July 2020
### restart on 8th July

## Following the post RFmix steps from Admix-simu_130519
# mudei o script da Alicia para incluir só chr 1 e 2
## nao estou usando isso!!
python2.7 /genetica_1/SGDP/ancestry_pipeline-master/collapse_ancestry.py --rfmix PEL_RFmixed_RefHGDP_1.2.Viterbi.txt --snp_locations PELadmixsimu_RFmixed_RefHGDP_chr1.snp_locations --fbk PEL_RFmixed_RefHGDP_1.2.ForwardBackward.txt --fbk_threshold 0.9 --ind PEL1 --ind_info ../PEL.noref --pop_labels NAT,EUR,AFR --out PEL1

sed 's/EUR/1/g' PEL1_A.bed > temp1_PEL1_A.bed
sed 's/NAT/0/g' temp1_PEL1_A.bed > temp2_PEL1_A.bed
sed 's/UNK/2/g' temp2_PEL1_A.bed > NA-PEL1.bed

sed 's/EUR/1/g' PEL1_B.bed > temp1_PEL1_B.bed
sed 's/NAT/0/g' temp1_PEL1_B.bed > temp2_PEL1_B.bed
sed 's/UNK/2/g' temp2_PEL1_B.bed > NB-PEL1.bed

# preparar Simu Real para verificar a acurácia ## 20 Julho

#Cuidado aqui!!! Jessica começ ar aqui. o Exemplo abaixo é do Brasa, precisa alterar os radicais para o Mix e para o Baiano 

for i in {1..1}; do awk '{print $1, $4}' ./chr${i}.pos > chr${i}.posmini; paste chr${i}.posmini ./Baiano${i}.hanc2 | sed 's/\t/ /g' > SimuReal_Baiano_chr${i}.AncCalls; done

# Prepare data for accuracy script 
##### pros eur_nat 2/0 1/1
##### pros afr_nat 3/1 2/0
# HGDPNAT_IBS_YRI - chamei de 1
export AdmixPOP=exNATadmix
for i in {1..1}; do sed 's/2/0/g' $AdmixPOP.RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_${i}.2.Viterbi.txt > $AdmixPOP.1_Lat2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' $AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_chr${i}.map | sed 's/:.*//g' > $AdmixPOP.1_mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.1_mappin_${i}_Lat $AdmixPOP.1_Lat2_${i} | sed 's/\t/ /g' > $AdmixPOP.1_Lat3_${i}; done

export AdmixPOP=exEURadmix
for i in {1..1}; do sed 's/2/0/g' $AdmixPOP.RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_${i}.2.Viterbi.txt > $AdmixPOP.1_Lat2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' $AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_chr${i}.map | sed 's/:.*//g' > $AdmixPOP.1_mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.1_mappin_${i}_Lat $AdmixPOP.1_Lat2_${i} | sed 's/\t/ /g' > $AdmixPOP.1_Lat3_${i}; done

export AdmixPOP=MAM
for i in {1..1}; do sed 's/2/0/g' $AdmixPOP.RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_${i}.2.Viterbi.txt > $AdmixPOP.1_Lat2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' $AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_chr${i}.map | sed 's/:.*//g' > $AdmixPOP.1_mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.1_mappin_${i}_Lat $AdmixPOP.1_Lat2_${i} | sed 's/\t/ /g' > $AdmixPOP.1_Lat3_${i}; done

export AdmixPOP=PEL
for i in {1..1}; do sed 's/2/0/g' $AdmixPOP.RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_${i}.2.Viterbi.txt > $AdmixPOP.1_Lat2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' $AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_chr${i}.map | sed 's/:.*//g' > $AdmixPOP.1_mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.1_mappin_${i}_Lat $AdmixPOP.1_Lat2_${i} | sed 's/\t/ /g' > $AdmixPOP.1_Lat3_${i}; done

export AdmixPOP=Baiano
for i in {1..1}; do sed 's/1/4/g' $AdmixPOP.RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_${i}.2.Viterbi.txt > $AdmixPOP.1_Lat2_${i}; done
for i in {1..1}; do sed 's/2/0/g' $AdmixPOP.1_Lat2_${i} > $AdmixPOP.1_Lat2.1_${i}; done
for i in {1..1}; do sed 's/3/1/g' $AdmixPOP.1_Lat2.1_${i} > $AdmixPOP.1_Lat2.2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' $AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_chr${i}.map | sed 's/:.*//g' > $AdmixPOP.1_mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.1_mappin_${i}_Lat $AdmixPOP.1_Lat2.2_${i} | sed 's/\t/ /g' > $AdmixPOP.1_Lat3_${i}; done

#30/03
export AdmixPOP=Baiano
for i in {1..1}; do sed 's/2/0/g' $AdmixPOP.2.RFmixed_Ref1kgYRI_Ref1kgIBS_${i}.2.Viterbi.txt > $AdmixPOP.1_Lat2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' $AdmixPOP.admixsimu_RFmixed_Ref1kgYRI_Ref1kgIBS_chr${i}.map | sed 's/:.*//g' > $AdmixPOP.mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.mappin_${i}_Lat $AdmixPOP.1_Lat2_${i} | sed 's/\t/ /g' > $AdmixPOP.1_Lat3_${i}; done

export AdmixPOP=PEL #depois MAM, exNATadmix, exEURadmix, Baiano
for i in {1..1}; do sed 's/2/0/g' $AdmixPOP.RFmixed_RefHGDPNAT_Ref1kgIBS_${i}.2.Viterbi.txt > $AdmixPOP.1_Lat2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' $AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_chr${i}.map | sed 's/:.*//g' > $AdmixPOP.1_mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.1_mappin_${i}_Lat $AdmixPOP.1_Lat2_${i} | sed 's/\t/ /g' > $AdmixPOP.1_Lat3_${i}; done
