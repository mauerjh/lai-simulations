###### Making reference files

### 31/03 PERGUNTAR SOBRE GSA VARS
### mudar pasta p nao sobrescrever arquivos
/genetica_1/Marq/teste_pops

# EUR_1KG, Iberic sample only 
egrep "IBS" ../../SGDP/1000GP_Phase3.sample | awk '{print $1}' > EUR_1kg
tail -n 77 EUR_1kg > EUR_1kg2 
egrep "FIN" ../../SGDP/1000GP_Phase3.sample | awk '{print $1}' >> EUR_1kg2
egrep "TSI" ../../SGDP/1000GP_Phase3.sample | awk '{print $1}' >> EUR_1kg2
egrep "GBR" ../../SGDP/1000GP_Phase3.sample | awk '{print $1}' >> EUR_1kg2

for i in {1..1}; do bcftools view ../n1kg_hgdp_short_chr1.vcf -S ./EUR_1kg2 > EUR_1kg_chr1.vcf; bgzip EUR_1kg_chr1.vcf; bcftools convert --hapsample REF_EUR_1KG_chr1 EUR_1kg_chr1.vcf.gz; done

# AFR - YRI alterado em 1810
egrep "YRI" ../../SGDP/1000GP_Phase3.sample | awk '{print $1}' > AFR_1kg
tail -n 77 AFR_1kg > AFR_1kg2 # to remove the other 30+1 that were included in the admix-simu
egrep "LWK" ../../SGDP/1000GP_Phase3.sample | awk '{print $1}' >> AFR_1kg2
egrep "GWD" ../../SGDP/1000GP_Phase3.sample | awk '{print $1}' >> AFR_1kg2
egrep "MSL" ../../SGDP/1000GP_Phase3.sample | awk '{print $1}' >> AFR_1kg2
egrep "ESN" ../../SGDP/1000GP_Phase3.sample | awk '{print $1}' >> AFR_1kg2

for i in {1..1}; do bcftools view ../n1kg_hgdp_short_chr1.vcf -S ./AFR_1kg2 > AFR_1kg_chr1.vcf; bgzip AFR_1kg_chr1.vcf; bcftools convert --hapsample REF_AFR_1KG_chr1 AFR_1kg_chr1.vcf.gz; done

#0810
for i in {1..1}; do awk '{gsub(/^chr/,""); print}' REF_EUR_1KG_chr1.hap > REF_EUR_1kg_chr${i}_nc.hap ; done
for i in {1..1}; do awk '{gsub(/^chr/,""); print}' REF_AFR_1KG_chr1.hap > REF_AFR_1kg_chr${i}_nc.hap ; done

# printar as amostras ref em um mesmo arquivo. #### J ok
# 2 way NAT eEUR

# awk 'FNR>2 {print $1}' REF_NAT_HGDP_chr1.sample > REF_NAT_HGDP 
awk 'FNR>2 {print $1}' REF_EUR_1KG_chr1.sample > REF_EUR_1kg
awk 'FNR>2 {print $1}' REF_AFR_1KG_chr1.sample > REF_AFR_1KG

cat ../REF_NAT_HGDP REF_EUR_1kg REF_AFR_1KG > ./REF_NAT_EUR_1kg_AFR_1KG

# ajustar o .sample como o programa pede
#0910 ?
awk '{print "0",$1,"0","0","0","1","1"}' REF_EUR_1KG_chr1.sample > ./REF_EUR_1kg.sample
awk '{print "0",$1,"0","0","0","1","1"}' REF_AFR_1KG_chr1.sample > ./REF_AFR_1KG.sample

# inserir a informacao do cabecalho na mao #03/04 jess alterei aqui em relação ao script original. estava o nome do arquivo errado
vi ./REF_EUR_1kg.sample
vi ./REF_AFR_1KG.sample


#exemplo do .sample final
ID_1 ID_2 missing father mother sex plink_pheno
0 0 0 D D D B
0 PEL1 0 0 0 1 1
0 PEL2 0 0 0 1 1
0 PEL3 0 0 0 1 1
0 PEL4 0 0 0 1 1

##################
# fazer para Baiano, PEL, MAM, exNATadmix e exEURadmix
export AdmixPOP=exNATadmix
cp /genetica_1/Jessica/admixsimu/PEL.noref $AdmixPOP.no
# Jessica cuidado aqui, precisa trocar o baiano
sed 's/PEL/Brasa/g' $AdmixPOP.no > $AdmixPOP.noref

### Testar na MIX ye BRASA; Se mudar rodar 2way-3way e se mudar rodar 2way-2way
# Prepare samples for RFMix
#23/03/21: Aqui usando como referencia o que era usado no three-way, porém nas simulações 2way

for i in {1..1}; do python /genetica_1/Jessica/admixsimu/ancestry_pipeline-master/shapeit2rfmix.py --shapeit_hap_ref ../REF_NAT_HGDP_chr${i}_nc.hap,./REF_EUR_1kg_chr${i}_nc.hap,./REF_AFR_1kg_chr${i}_nc.hap --shapeit_hap_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.haps --shapeit_sample_ref ../REF_NAT_HGDP.sample,./REF_EUR_1kg.sample,./REF_AFR_1KG.sample --shapeit_sample_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.sample --chr ${i} --genetic_map ../../Jessica/admixsimu/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt --ref_keep ./REF_NAT_EUR_1kg_AFR_1KG --admixed_keep ./$AdmixPOP.noref --out $AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI; done

###
# Finally Run RFmix. Tem que rodar da pasta do RFmix. O que a Jessiaca baixou ainda nao esta instalado. Estou rodando no
cd /genetica_1/SGDP/RFMix_v1.5.4/

for i in {1..1}; do nohup python RunRFMix.py -e 2 -w 0.2 -n 5 -G 9 --num-threads 20 --forward-backward --use-reference-panels-in-EM TrioPhased /genetica_1/Marq/teste_pops/$AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_chr1.alleles /genetica_1/Marq/teste_pops/$AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI.classes /genetica_1/Marq/teste_pops/$AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_chr1.snp_locations -o /genetica_1/Marq/teste_pops/$AdmixPOP.RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_${i}; done &

#Cuidado aqui!!! Jessica começ ar aqui. o Exemplo abaixo é do Brasa, precisa alterar os radicais para o Mix e para o Baiano 
cd /genetica_1/Marq/teste_pops
for i in {1..1}; do awk '{print $1, $4}' ../chr${i}.pos > chr${i}.posmini; paste chr${i}.posmini ../PEL${i}.hanc2 | sed 's/\t/ /g' > SimuReal_PEL_chr${i}.AncCalls; done

# Prepare data for accuracy script 
#RFMIXED: 1 - NAT, 2 - EUR, 3 - AFR
# SimuReal 0 - NAT, 1 -EUR, 2 - AFR
# usar o script do R JM_Accus_2020.R (faz o gsub de 1 por 3 depois)
# brasa apenas 20 indivs
# HGDPNAT_IBS_YRI - chamei de 1
export AdmixPOP=Brasa
for i in {1..1}; do sed 's/1/0/g' $AdmixPOP.RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_${i}.2.Viterbi.txt > $AdmixPOP.1_Lat2_${i}; done
for i in {1..1}; do sed 's/2/1/g' $AdmixPOP.1_Lat2_${i} > $AdmixPOP.1_Lat2.1_${i}; done
for i in {1..1}; do sed 's/3/2/g' $AdmixPOP.1_Lat2.1_${i} > $AdmixPOP.1_Lat2.2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' $AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_chr${i}.map | sed 's/:.*//g' > $AdmixPOP.1_mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.1_mappin_${i}_Lat $AdmixPOP.1_Lat2.2_${i} | sed 's/\t/ /g' > $AdmixPOP.1_Lat3_${i}; done

###########
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
