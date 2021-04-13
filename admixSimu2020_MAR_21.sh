### Downloads
####
# genetic maps from 1KG hg38
wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip

################ 1KG HG38 ##################
# 1) Download the Reference haplotypes sample from here http://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3/
# Download dos arquivos vcfs referência hg38: ALL.chr*.vcf.gz
for i in {1..22}; do wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz
for i in {1..22}; do wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz.tbi

################ HGDP HG38 ##################
## Anders shared the HGDP data in a single file at ftp://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/statphase/
# Break in chrs
for i in {1..22}; do bcftools view hgdp_wgs.20190516.statphase.autosomes.vcf.gz --regions chr${i} -Oz -o /genetica_1/Marq/hgdp_phased_chr${i}_b38.vcf.gz; done

## Subset of samples from HGDP: EUROPE (Adygei, Basque, French, North Italian, Orcadian, Russian, Sardinian, Tuscan); AMERICA (Colombian, Karitiana, Maya, Pima, Surui); AFRICA (Bantu, Biaka, Mandenka, Mbuti, Mozabite, San, Yoruba)

# All metadata
egrep "EUROPE|AMERICA|AFRICA" ./hgdp_wgs.20190516.metadata.txt > Haps_interesse_HGDP

# ID only
awk ‘{print $2}’ Haps_interesse_HGDP >  samples_interesse

# Filter the variants with missing info
for i in {1..22}; do bcftools view ./hgdp_phased_chr${i}_b38.vcf.gz -S ./samples_interesse -e 'GT[*] = "mis"' -Oz -o hgdp_phased_chr${i}_QCed_b38.vcf.gz; done

## Break into two datasets, one for Admix-simu and another for RFMix. The selection was manually done. First half of each subpop for admixsimu and the remaining for RFmix.
awk '{print $1}' Sample_HGDP_admixSimu.txt > Sample_HGDP_admixSimu2.txt
awk '{print $1}' Sample_HGDP_RFmix.txt > Sample_HGDP_RFmix2.txt

for i in {1..22}; do bcftools view hgdp_phased_chr${i}_QCed_b38.vcf.gz -S Sample_HGDP_admixSimu2.txt -Oz -o hgdp_phased_chr${i}_QCed_b38_AdmixSimu.vcf.gz; bcftools view hgdp_phased_chr${i}_QCed_b38.vcf.gz -S Sample_HGDP_RFmix2.txt -Oz -o hgdp_phased_chr${i}_QCed_b38_RFMix.vcf.gz; done

#####
# merging both datasets:
#####

#
for i in {1..22}; do bcftools index hgdp_phased_chr${i}_QCed_b38.vcf.gz

# Intersection of variants between HGDP and 1KG: took a lot of time. Find a faster way next time.
for i in {1..22}; do echo "bcftools isec -n +2 -p curto_chr${i} ./1KG/ALL.chr${i}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz ../../Marq/hgdp_phased_chr${i}_QCed_b38.vcf.gz -Oz; cd curto_chr${i} ; bcftools merge 0000.vcf.gz 0001.vcf.gz -o 1kg_hgdp_short_chr${i}.vcf.gz -Oz; rm 0000.vcf.gz ; rm 0001.vcf.gz"; done| parallel -j 3


# adicionei em 5 de outubro para fazer o QC. Besteira pq se vai pro plink perde o phasing
#
# best way to QC the vcf and keep it phased
bcftools view -q 0.01:minor ../Jessica/admixsimu/curto_chr1/1kg_hgdp_short_chr1.vcf -Ov -o n1kg_hgdp_short_chr1.vcf

# GSA: extracting vars. I used a .bim from our cohort the were genotype with Global Screening array.
awk -v OFS="\t" '{print "chr"$1,$4}' Childs_BI20.bim > GSAvars.txt
bcftools view -q 0.01:minor -T GSAvars.txt ../Jessica/admixsimu/curto_chr1/1kg_hgdp_short_chr1.vcf -Ov -o nn1kg_hgdp_short_chr1.vcf

########################################################
###### Admixsimu Step, Preparing reference data and generating admixed genomes ######
########################################################
# Break into Pops for Admixsimu
## Separate each pop from HGDP
grep "AFRICA" Sample_HGDP_admixSimu.txt | awk '{print $1}' > AFR_HGDP_samp_AdmixSimu2 # if using all AFR from HGDP as ref
grep "AMERICA" Sample_HGDP_admixSimu.txt | awk '{print $1}' > NAT_HGDP_samp_AdmixSimu2
grep "EUROPE" Sample_HGDP_admixSimu.txt |awk '{print $1}' > EUR_HGDP_samp_AdmixSimu2 # if using all EUR from HGDP as ref

## Separate each pop from 1kg
egrep "IBS" 1000GP_Phase3.sample |wc
egrep "IBS" 1000GP_Phase3.sample |head -n 30 | awk '{print $1}' > novoEUR_IBR # I used only 30 samples to keep the same size of NAT, you may change
egrep "YRI" 1000GP_Phase3.sample |wc
egrep "YRI" 1000GP_Phase3.sample |head -n 31 | awk '{print $1}' > novoAFR_YRI # I used only 30 samples to keep the same size of NAT, you may change. One sample is not present in the 1kg_hgdp vcf file, removed manually.


# Extract the samples for each pop
## Note that I misslabled the EUR and AFR as "hgdp_*_chr${i}_AdmixSimu.vcf.gz" but their samples were from 1000 genomes
for i in {1..1};do bcftools view n1kg_hgdp_short_chr1.vcf -S NAT_HGDP_samp_AdmixSimu2 -Oz -o hgdp_NAT_chr${i}_AdmixSimu.vcf.gz; bcftools view n1kg_hgdp_short_chr1.vcf -S novoEUR_IBR2 -Oz -o hgdp_EUR_chr${i}_AdmixSimu.vcf.gz;bcftools view n1kg_hgdp_short_chr1.vcf -S novoAFR_YRI -Oz -o hgdp_AFR_chr${i}_AdmixSimu.vcf.gz;done


## Convert into haps format (Just as a remainder for the future that cost me a lot of time here, the option is --hapsample and not -hapsample, if only with one dash bcftools understand it as -h which is --haplegendsample option, at the end hap files hadn't the variant info)
for i in {1..1}; do bcftools convert hgdp_AFR_chr${i}_AdmixSimu.vcf.gz --hapsample hgdp_AFR_chr${i}_AdmixSimu;bcftools convert hgdp_EUR_chr${i}_AdmixSimu.vcf.gz --hapsample hgdp_EUR_chr${i}_AdmixSimu; bcftools convert hgdp_NAT_chr${i}_AdmixSimu.vcf.gz --hapsample hgdp_NAT_chr${i}_AdmixSimu; done

## convert into phgeno format
for i in {1..1}; do gunzip hgdp_*_chr${i}_AdmixSimu.hap.gz;done

for i in {1..1};do cut -d' ' -f 6- hgdp_NAT_chr${i}_AdmixSimu.hap | sed 's/\s//g' > NAT${i}.phgeno; cut -d' ' -f 6- hgdp_AFR_chr${i}_AdmixSimu.hap | sed 's/\s//g' > AFR${i}.phgeno; cut -d' ' -f 6- hgdp_EUR_chr${i}_AdmixSimu.hap | sed 's/\s//g' > EUR${i}.phgeno;done

## create the .snp files
for i in {1..1};do ~/Documents/Marquinhos/bino/plink2 --vcf hgdp_NAT_chr${i}_AdmixSimu.vcf.gz --const-fid 0 --max-alleles 2 --make-bed --out hgdp_pBim2_${i};done

for i in {1..1};do perl ../SGDP/admix-simu-master/insert-map.pl hgdp_pBim2_${i}.bim ../Jessica/admixsimu/genetic_map/genetic_map_chr${i}_hg38.txt > chr${i}.pos; awk -F' ' '{ print $1":"$4"_"$5"_"$6, $1, $3, $4, $5, $6 }' chr${i}.pos > hgdp.chr${i}.snp; done

# MINU: Puxei isso lá de baixo... Checar se precisa subir mais ou descer
for i in {1..1}; do awk '{print $4,$3*100 }' hgdp.chr${i}.snp > hgdp_chr${i}.snp_loc ;done
head hgdp_chr${i}.snp_loc

##################
# RFmix 18/06/2020 - 02/07/2020
#################

#### splitting the refs: 1) 1kg only: EUR = IBS, TSI, CEU, GBR;  NAT = PEL 2) 1kg e hgdp: EUR = IBS, TSI, CEU, GBR , HGDP-EUR ; NAT = HGDP-NAT 3) 1kg only: EUR = IBS, TSI, CEU, GBR;  NAT = PEL + EAS
pasta Jessica/admixsimu/ID_REF

# EUR_1KG, Iberic sample only #MUDAR AQUI TESRTE egrep "eur" ../1000GP_Phase3.sample | awk '{print $1}' > EUR_1kg
egrep "IBS" ../1000GP_Phase3.sample | awk '{print $1}' > EUR_1kg

tail -n 77 EUR_1kg > EUR_1kg2 # to remove the other 30 that were included in the admix-simu
for i in {1..1}; do bcftools view n1kg_hgdp_short_chr1.vcf -S ./EUR_1kg2 > EUR_1kg_chr1.vcf; bgzip EUR_1kg_chr1.vcf; bcftools convert --hapsample REF_EUR_1KG_chr1 EUR_1kg_chr1.vcf.gz; done

for i in {1..1}; do bcftools view n1kg_hgdp_short_chr1.vcf -S ./EUR_1kg2 -T GSAvars.txt > gsaEUR_1kg_chr1.vcf; bgzip gsaEUR_1kg_chr1.vcf; bcftools convert --hapsample gsaREF_EUR_1KG_chr1 gsaEUR_1kg_chr1.vcf.gz; done

# AFR - YRI alterado em 1810
egrep "YRI" ../SGDP/1000GP_Phase3.sample | awk '{print $1}' > AFR_1kg
tail -n 77 AFR_1kg > AFR_1kg2 # to remove the other 30+1 that were included in the admix-simu
for i in {1..1}; do bcftools view n1kg_hgdp_short_chr1.vcf -S ./AFR_1kg2 > AFR_1kg_chr1.vcf; bgzip AFR_1kg_chr1.vcf; bcftools convert --hapsample REF_AFR_1KG_chr1 AFR_1kg_chr1.vcf.gz; done

for i in {1..1}; do bcftools view n1kg_hgdp_short_chr1.vcf -S ./AFR_1kg2 -T GSAvars.txt > gsaAFR_1kg_chr1.vcf; bgzip gsaAFR_1kg_chr1.vcf; bcftools convert --hapsample gsaREF_AFR_1KG_chr1 gsaAFR_1kg_chr1.vcf.gz; done


# PEL_1KG
egrep "PEL" ../1000GP_Phase3.sample | awk '{print $1}' > PEL_1kg

bcftools view n1kg_hgdp_short_chr${i}.vcf -S ../Jessica/admixsimu/ID_REF/PEL_1kg > PEL_1kg_chr${i}.vcf; bgzip PEL_1kg_chr${i}.vcf; bcftools convert --hapsample REF_PEL_1KG_chr${i} PEL_1kg_chr${i}.vcf.gz

bcftools view n1kg_hgdp_short_chr${i}.vcf -S ../Jessica/admixsimu/ID_REF/PEL_1kg -T GSAvars.txt > gsaPEL_1kg_chr${i}.vcf; bgzip gsaPEL_1kg_chr${i}.vcf; bcftools convert --hapsample gsaREF_PEL_1KG_chr${i} gsaPEL_1kg_chr${i}.vcf.gz


# PEL_EAS_1KG
egrep "PEL|EAS" ../SGDP/1000GP_Phase3.sample | awk '{print $1}' > PEL_EAS_1KG

# só adicionei o n na frente do 1kg, 0510, 0
#0910
for i in {1..1}; do bcftools view n1kg_hgdp_short_chr${i}.vcf -S ./PEL_EAS_1KG > PEL_EAS_1KG_chr${i}.vcf; bgzip PEL_EAS_1KG_chr${i}.vcf; bcftools convert --hapsample REF_PEL_EAS_1KG_chr${i} PEL_EAS_1KG_chr${i}.vcf.gz; done

for i in {1..1}; do bcftools view n1kg_hgdp_short_chr${i}.vcf -S ./PEL_EAS_1KG -T GSAvars.txt > gsaPEL_EAS_1KG_chr${i}.vcf; bgzip gsaPEL_EAS_1KG_chr${i}.vcf; bcftools convert --hapsample gsaREF_PEL_EAS_1KG_chr${i} gsaPEL_EAS_1KG_chr${i}.vcf.gz; done

#0810
# dei uma editada para tirar o Colombian ***70 que entrou no admixsimu
egrep "Karitiana|Pima|Colombian" ../Jessica/admixsimu/HGDP/hgdp_wgs.20190516.metadata.txt | awk '{print $1}' > NAT_HGDP
for i in {1..1}; do bcftools view n1kg_hgdp_short_chr${i}.vcf -S NAT_HGDP > NAT_HGDP_chr${i}.vcf; bgzip NAT_HGDP_chr${i}.vcf; bcftools convert --hapsample REF_NAT_HGDP_chr${i} NAT_HGDP_chr${i}.vcf.gz;done

for i in {1..1}; do bcftools view n1kg_hgdp_short_chr${i}.vcf -S NAT_HGDP -T GSAvars.txt > gsaNAT_HGDP_chr${i}.vcf; bgzip gsaNAT_HGDP_chr${i}.vcf; bcftools convert --hapsample gsaREF_NAT_HGDP_chr${i} gsaNAT_HGDP_chr${i}.vcf.gz;done


# Note que tem um erro de typo no _nc europeu e afr (1810) chamando ele de HGDP tbm
#0810
awk '{gsub(/^chr/,""); print}' REF_NAT_HGDP_chr${i}.hap > REF_NAT_HGDP_chr${i}_nc.hap
awk '{gsub(/^chr/,""); print}' REF_EUR_1KG_chr1.hap > REF_EUR_1kg_HGDP_chr${i}_nc.hap
awk '{gsub(/^chr/,""); print}' REF_PEL_1KG_chr${i}.hap > REF_PEL_1KG_chr${i}_nc.hap
awk '{gsub(/^chr/,""); print}' REF_AFR_1KG_chr1.hap > REF_AFR_1kg_HGDP_chr${i}_nc.hap

awk '{gsub(/^chr/,""); print}' gsaREF_NAT_HGDP_chr${i}.hap > gsaREF_NAT_HGDP_chr${i}_nc.hap
awk '{gsub(/^chr/,""); print}' gsaREF_EUR_1KG_chr1.hap > gsaREF_EUR_1kg_HGDP_chr${i}_nc.hap
awk '{gsub(/^chr/,""); print}' gsaREF_PEL_1KG_chr${i}.hap > gsaREF_PEL_1KG_chr${i}_nc.hap
awk '{gsub(/^chr/,""); print}' gsaREF_AFR_1KG_chr1.hap > gsaREF_AFR_1kg_HGDP_chr${i}_nc.hap
# 0910
for i in {1..1}; do awk '{gsub(/^chr/,""); print}' REF_PEL_EAS_1KG_chr${i}.hap > REF_PEL_EAS_1KG_chr${i}_nc.hap ; done

for i in {1..1}; do awk '{gsub(/^chr/,""); print}' gsaREF_PEL_EAS_1KG_chr${i}.hap > gsaREF_PEL_EAS_1KG_chr${i}_nc.hap ; done

# printar as amostras ref em um mesmo arquivo.
# 2 way NAT eEUR

awk 'FNR>2 {print $1}' REF_PEL_EAS_1KG_chr1.sample > REF_PEL_EAS_1KG; cat REF_PEL_EAS_1KG REF_EUR_1kg > ./REF_PEL_EAS_1KG_REF_EUR_1kg

awk 'FNR>2 {print $1}' REF_NAT_HGDP_chr1.sample > REF_NAT_HGDP; awk 'FNR>2 {print $1}' REF_EUR_1KG_chr1.sample > REF_EUR_1kg; cat REF_NAT_HGDP REF_EUR_1kg > ../REF_NAT_EUR_1kg

awk 'FNR>2 {print $1}' REF_PEL_1KG_chr1.sample > REF_PEL_1KG; awk 'FNR>2 {print $1}' REF_EUR_1KG_chr1.sample > REF_EUR_1kg; cat REF_PEL_1KG REF_EUR_1kg > ../REF_PEL_1KG_REF_EUR_1kg

# 3 way
awk 'FNR>2 {print $1}' REF_AFR_1KG_chr1.sample > REF_AFR_1KG; cat REF_PEL_1KG REF_EUR_1kg REF_AFR_1KG > ./REF_PEL_1KG_REF_EUR_1kg_REF_AFR_1KG
cat REF_NAT_HGDP REF_EUR_1kg REF_AFR_1KG > ./REF_NAT_EUR_1kg_AFR_1KG
cat REF_PEL_EAS_1KG REF_EUR_1kg REF_AFR_1KG > ./REF_PEL_EAS_1KG_REF_EUR_1kg_AFR_1KG

# 2 way AFR e EUR

cat REF_EUR_1kg REF_AFR_1KG > ./REF_EUR_1kg_REF_AFR_1KG

#gsa

awk 'FNR>2 {print $1}' gsaREF_PEL_EAS_1KG_chr1.sample > gsaREF_PEL_EAS_1KG; cat gsaREF_PEL_EAS_1KG gsaREF_EUR_1kg > ./gsaREF_PEL_EAS_1KG_REF_EUR_1kg

awk 'FNR>2 {print $1}' gsaREF_NAT_HGDP_chr1.sample > gsaREF_NAT_HGDP; awk 'FNR>2 {print $1}' gsaREF_EUR_1KG_chr1.sample > gsaREF_EUR_1kg; cat gsaREF_NAT_HGDP gsaREF_EUR_1kg > ./gsaREF_NAT_EUR_1kg
# MINU: achei um erro aqui faltando um "gsa na frente do REF_EUR_1kg. Checar consequencias.
awk 'FNR>2 {print $1}' gsaREF_PEL_1KG_chr1.sample > gsaREF_PEL_1KG; awk 'FNR>2 {print $1}' gsaREF_EUR_1KG_chr1.sample > gsaREF_EUR_1kg; cat gsaREF_PEL_1KG REF_EUR_1kg > ./gsaREF_PEL_1KG_REF_EUR_1kg

awk 'FNR>2 {print $1}' gsaREF_AFR_1KG_chr1.sample > gsaREF_AFR_1KG; cat gsaREF_EUR_1kg gsaREF_AFR_1KG > ./gsaREF_EUR_1kg_REF_AFR_1KG

# alvos, arquivo criado na mão apenas com nomes ficticios da nossa amostra alvo
cp ../curto_chr1/noref.noref ../PEL.noref
cp /genetica/Jessica/admixsimu/PEL.noref $AdmixPOP.no
sed 's/PEL/Baiano/g' $AdmixPOP.no > $AdmixPOP.noref

# ajustar o .sample como o programa pede
#0910 ?
awk '{print "0",$1,"0","0","0","1","1"}' REF_PEL_EAS_1KG_chr1.sample > ./REF_PEL_EAS_1KG.sample

#0810
awk '{print "0",$1,"0","0","0","1","1"}' REF_NAT_HGDP_chr1.sample > ./REF_NAT_HGDP.sample
awk '{print "0",$1,"0","0","0","1","1"}' REF_EUR_1KG_chr1.sample > ./REF_EUR_1kg.sample
awk '{print "0",$1,"0","0","0","1","1"}' REF_PEL_1KG_chr1.sample > ./REF_PEL_1KG.sample
awk '{print "0",$1,"0","0","0","1","1"}' REF_AFR_1KG_chr1.sample > ./REF_AFR_1KG.sample

#gsa
awk '{print "0",$1,"0","0","0","1","1"}' gsaREF_PEL_EAS_1KG_chr1.sample > ./gsaREF_PEL_EAS_1KG.sample

#0810
awk '{print "0",$1,"0","0","0","1","1"}' gsaREF_NAT_HGDP_chr1.sample > ./gsaREF_NAT_HGDP.sample
awk '{print "0",$1,"0","0","0","1","1"}' gsaREF_EUR_1KG_chr1.sample > ./gsaREF_EUR_1kg.sample
awk '{print "0",$1,"0","0","0","1","1"}' gsaREF_PEL_1KG_chr1.sample > ./gsaREF_PEL_1KG.sample

# inserir a informacao do cabecalho na mao
vi ../REF_EUR_1kg_HGDP.sample
vi ../REF_PEL_EAS_1KG.sample

vi ../REF_EUR_HGDP.sample
vi ../REF_NAT_HGDP.sample

vi ../REF_PEL_1KG.sample
vi ./REF_AFR_1KG.sample

#exemplo do .sample final
ID_1 ID_2 missing father mother sex plink_pheno
0 0 0 D D D B
0 PEL1 0 0 0 1 1
0 PEL2 0 0 0 1 1
0 PEL3 0 0 0 1 1
0 PEL4 0 0 0 1 1

### Jessica, começa aqui
#####################
# Generating the Admixed genomes
#####################

cd /genetica_1/Marq

## Edit here the ancestry orders and which population will you test, remember to manually write the .dat file.
export  anc1=
export  anc2=
export  anc3=
export  AdmixPOP=

# 0/50/50 = Baiano
# 0/70/30 = Pardo
# 0/30/70 = Afam
# 33/33/34 = Mix

# editar o arquivo Baiano.dat

# 2 way admixture
for i in {1..1}; do ../SGDP/admix-simu-master/simu-mix.pl ./$AdmixPOP.dat hgdp.chr${i}.snp $AdmixPOP${i} -$anc1 $anc1${i}.phgeno -$anc2 $anc2${i}.phgeno ; done

# 3 way admixture
for i in {1..1}; do ../SGDP/admix-simu-master/simu-mix.pl ./$AdmixPOP.dat hgdp.chr${i}.snp $AdmixPOP${i} -$anc1 $anc1${i}.phgeno -$anc2 $anc2${i}.phgeno -$anc3 $anc3${i}.phgeno ; done


# convert .bp from the simulation to .hanc
for i in {1..1}; do ../Jessica/admixsimu/admix-simu-master/bp2anc.pl $AdmixPOP${i}.bp > $AdmixPOP${i}.hanc; done

# hanc1
for i in {1..1}; do sed 's/./& /g' $AdmixPOP${i}.hanc > $AdmixPOP${i}.hanc1 ;done

# transpose hanc1
for i in {1..1};do python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < $AdmixPOP${i}.hanc1 > $AdmixPOP${i}.hanc2; done

head $AdmixPOP${i}.hanc2
rm *.hanc1

# From Elizabeth's scripts: Now I need to separate the .phgeno file to have a space between each character
for i in {1..1}; do sed 's/./& /g' ./$AdmixPOP${i}.phgeno > temp_Genos_chr${i}; done
for i in {1..1}; do awk -F' ' '{ print $2, $1, $4, $5, $6 }' hgdp.chr${i}.snp > temp_SNP_chr${i} ; done
for i in {1..1}; do paste temp_SNP_chr${i} temp_Genos_chr${i} > temp_admix_chr${i}.haps ; done
for i in {1..1}; do sed 's/\t/ /g' temp_admix_chr${i}.haps > $AdmixPOP.chr${i}.haps ; done
head $AdmixPOP.chr1.haps

# Manually created the .sample file with fake names
# Jessica cuidado aqui
# MINU: MANO,
# nao consegui substituir o sed, fazer manualmente
sed 's/PEL/Baiano/g' ../Jessica/admixsimu/PEL.sample.txt > $AdmixPOP.sample.txt
for i in {1..1}; do cp $AdmixPOP.sample.txt $AdmixPOP.chr${i}.sample; done


#####################
#Prepare the haps files to RFmix (from Alicia Martin):
#############
# Jessica recomece aqui

cp /genetica/Jessica/admixsimu/PEL.noref $AdmixPOP.no
# Jessica cuidado aqui, precisa trocar o baiano
sed 's/PEL/Baiano/g' $AdmixPOP.no > $AdmixPOP.noref

###
# Jessica rode o Afam e o Pardo aqui (acho que nao precis mudar nada
# 2 way AFR e EUR
for i in {1..1}; do python /genetica_1/Jessica/admixsimu/ancestry_pipeline-master/shapeit2rfmix.py --shapeit_hap_ref REF_EUR_1kg_HGDP_chr${i}_nc.hap,REF_AFR_1kg_HGDP_chr${i}_nc.hap --shapeit_hap_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.haps --shapeit_sample_ref ./REF_AFR_1KG.sample,./REF_EUR_1kg.sample --shapeit_sample_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.sample --chr ${i} --genetic_map ../Jessica/admixsimu/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt --ref_keep /./REF_EUR_1kg_REF_AFR_1KG --admixed_keep ./$AdmixPOP.noref --out $AdmixPOP.admixsimu_RFmixed_Ref1kgYRI_Ref1kgIBS; done


####
# Jessica aqui é 3 way, um pouco mais de cuidado

# 3 way para Brasa e Mix
for i in {1..1}; do python /genetica_1/Jessica/admixsimu/ancestry_pipeline-master/shapeit2rfmix.py --shapeit_hap_ref REF_PEL_1KG_chr${i}_nc.hap,REF_EUR_1kg_HGDP_chr${i}_nc.hap,REF_AFR_1kg_HGDP_chr${i}_nc.hap --shapeit_hap_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.haps --shapeit_sample_ref ./REF_PEL_1KG.sample,./REF_EUR_1kg.sample,REF_AFR_1KG.sample --shapeit_sample_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.sample --chr ${i} --genetic_map ../Jessica/admixsimu/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt --ref_keep /genetica_1/Marq/REF_PEL_1KG_REF_EUR_1kg_REF_AFR_1KG --admixed_keep ./$AdmixPOP.noref --out $AdmixPOPadmixsimu_RFmixed_Ref1kgPEL_Ref1kgIBS_Ref1kgYRI; done

#1)
#1810
for i in {1..1}; do python /genetica_1/Jessica/admixsimu/ancestry_pipeline-master/shapeit2rfmix.py --shapeit_hap_ref REF_NAT_HGDP_chr${i}_nc.hap,REF_EUR_1kg_HGDP_chr${i}_nc.hap,REF_AFR_1kg_HGDP_chr${i}_nc.hap --shapeit_hap_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.haps --shapeit_sample_ref ./REF_NAT_HGDP.sample,./REF_EUR_1kg.sample,REF_AFR_1KG.sample --shapeit_sample_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.sample --chr ${i} --genetic_map ../Jessica/admixsimu/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt --ref_keep ./REF_NAT_EUR_1kg_AFR_1KG --admixed_keep ./$AdmixPOP.noref --out $AdmixPOPadmixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI; done

# 3) 1810
for i in {1..1}; do python /genetica_1/Jessica/admixsimu/ancestry_pipeline-master/shapeit2rfmix.py --shapeit_hap_ref REF_PEL_EAS_1KG_chr${i}_nc.hap,REF_EUR_1kg_HGDP_chr${i}_nc.hap,REF_AFR_1kg_HGDP_chr${i}_nc.hap --shapeit_hap_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.haps --shapeit_sample_ref ./REF_PEL_EAS_1KG.sample,./REF_EUR_1kg.sample,REF_AFR_1KG.sample --shapeit_sample_admixed /genetica_1/Marq/$AdmixPOP.chr${i}.sample --chr ${i} --genetic_map ../Jessica/admixsimu/genetic_map/genetic_map_chr${i}_hg38_firstLine.txt --ref_keep ./REF_PEL_EAS_1KG_REF_EUR_1kg_AFR_1KG --admixed_keep ./$AdmixPOP.noref --out $AdmixPOPadmixsimu_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_Ref1kgYRI; done

###
# Finally Run RFmix. Tem que rodar da pasta do RFmix. O que a Jessiaca baixou ainda nao esta instalado. Estou rodando no
cd /genetica_1/SGDP/RFMix_v1.5.4/

#

# 2 way Afam e Pardo e Baiano
for i in {1..1}; do python RunRFMix.py -e 2 -w 0.2 -n 5 -G 9 --num-threads 30 --forward-backward --use-reference-panels-in-EM TrioPhased /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_Ref1kgYRI_Ref1kgIBS_chr1.alleles /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_Ref1kgYRI_Ref1kgIBS.classes /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_Ref1kgYRI_Ref1kgIBS_chr1.snp_locations -o /genetica_1/Marq/$AdmixPOP.RFmixed_Ref1kgYRI_Ref1kgIBS_1; done

# outros testes de Referencia para AFRxEUR. Pensei em AFR-AME e Yoruba x sudeste AFR

#####
# 3 way Mix e Brasa

for i in {1..1}; do nohup python RunRFMix.py -e 2 -w 0.2 -n 5 -G 9 --num-threads 30 --forward-backward --use-reference-panels-in-EM TrioPhased /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_Ref1kgPEL_Ref1kgIBS_Ref1kgYRI_chr1.alleles /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_Ref1kgPEL_Ref1kgIBS_Ref1kgYRI.classes /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_Ref1kgPEL_Ref1kgIBS_Ref1kgYRI_chr1.snp_locations -o /genetica_1/Marq/$AdmixPOP_RFmixed_Ref1kgPEL_Ref1kgIBS_Ref1kgYRI${i}; done

#1810 # inseri o use in ref panel in EM em 1910
for i in {1..1}; do nohup python RunRFMix.py -e 2 -w 0.2 -n 5 -G 9 --num-threads 30 --forward-backward --use-reference-panels-in-EM TrioPhased /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_chr1.alleles /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI.classes /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_chr1.snp_locations -o /genetica_1/Marq/$AdmixPOP_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_${i}; done &

#18,b10
for i in {1..1}; do nohup python RunRFMix.py -e 2 -w 0.2 -n 5 -G 9 --num-threads 30 --forward-backward --use-reference-panels-in-EM TrioPhased /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_Ref1kgYRI_chr1.alleles /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_Ref1kgYRI.classes /genetica_1/Marq/$AdmixPOP.admixsimu_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_Ref1kgYRI_chr1.snp_locations -o /genetica_1/Marq/$AdmixPOP_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_Ref1kgYRI_${i}; done &


# Jessica Parar aqui



### I am here on 3rd July 2020
### restart on 8th July

## Following the post RFmix steps from Admix-simu_130519
# mudei o script da Alicia para incluir só chr 1 e 2
## nao estou usando isso!!
python2.7 /genetica_1/SGDP/ancestry_pipeline-master/collapse_ancestry.py --rfmix PEL_RFmixed_RefHGDP_1.2.Viterbi.txt --snp_locations PELadmixsimu_RFmixed_RefHGDP_chr1.snp_locations --fbk PEL_RFmixed_RefHGDP_1.2.ForwardBackward.txt --fbk_threshold 0.9 --ind PEL1 --ind_info ../PEL.noref --pop_labels NAT,EUR,AFR --out PEL1
3
sed 's/EUR/1/g' PEL1_A.bed > temp1_PEL1_A.bed
sed 's/NAT/0/g' temp1_PEL1_A.bed > temp2_PEL1_A.bed
sed 's/UNK/2/g' temp2_PEL1_A.bed > NA-PEL1.bed

sed 's/EUR/1/g' PEL1_B.bed > temp1_PEL1_B.bed
sed 's/NAT/0/g' temp1_PEL1_B.bed > temp2_PEL1_B.bed
sed 's/UNK/2/g' temp2_PEL1_B.bed > NB-PEL1.bed

# preparar Simu Real para verificar a acurácia ## 20 Julho

#Jessica começ ar aqui. o Exemplo abaixo é do Brasa, precisa alterar os radicais para o Mix e para o Baiano - OK 13 Mar 2021

#Brasa
for i in {1..1}; do awk '{print $1, $4}' ./chr${i}.pos > chr${i}.posmini; paste chr${i}.posmini ./Mix${i}.hanc2 | sed 's/\t/ /g' > SimuReal_Mix_chr${i}.AncCalls; done

for i in {1..1}; do awk '{print $1, $4}' ./chr${i}.pos > chr${i}.posmini; paste chr${i}.posmini ./Baiano${i}.hanc2 | sed 's/\t/ /g' > SimuReal_Baiano_chr${i}.AncCalls; done

# Prepare data to accuracy script # Mix

#####
# Brasa
# chamei de 1
for i in {1..1}; do sed 's/1/0/g' Brasa_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_${i}.2.Viterbi.txt > Brasa_1_Lat2_${i}; done
for i in {1..1}; do sed 's/2/1/g' Brasa_1_Lat2_${i} > Brasa_1_Lat2.1_${i}; done
for i in {1..1}; do sed 's/3/2/g' Brasa_1_Lat2.1_${i} > Brasa_1_Lat2.2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' Brasaadmixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_chr${i}.map | sed 's/:.*//g' > Brasa_1_mappin_${i}_Lat; done
for i in {1..1}; do paste Brasa_1_mappin_${i}_Lat Brasa_1_Lat2.2_${i} | sed 's/\t/ /g' > Brasa_1_Lat3_${i}; done


# Mix
# chamei de 1
for i in {1..1}; do sed 's/1/0/g' Mix.RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_${i}.2.Viterbi.txt > Mix_1_Lat2_${i}; done
for i in {1..1}; do sed 's/2/1/g' Mix_1_Lat2_${i} > Mix_1_Lat2.1_${i}; done
for i in {1..1}; do sed 's/3/2/g' Mix_1_Lat2.1_${i} > Mix_1_Lat2.2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' Mix.admixsimu_RFmixed_RefHGDPNAT_Ref1kgIBS_Ref1kgYRI_chr${i}.map | sed 's/:.*//g' > Mix_1_mappin_${i}_Lat; done
for i in {1..1}; do paste Mix_1_mappin_${i}_Lat Mix_1_Lat2.2_${i} | sed 's/\t/ /g' > Mix_1_Lat3_${i}; done

# chamei de 2
for i in {1..1}; do sed 's/1/0/g' Mix.RFmixed_Ref1kgPEL_Ref1kgIBS_Ref1kgYRI${i}.2.Viterbi.txt > Mix_2_Lat2_${i}; done
for i in {1..1}; do sed 's/2/1/g' Mix_2_Lat2_${i} > Mix_2_Lat2.1_${i}; done
for i in {1..1}; do sed 's/3/2/g' Mix_2_Lat2.1_${i} > Mix_2_Lat2.2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' Mix.admixsimu_RFmixed_Ref1kgPEL_Ref1kgIBS_Ref1kgYRI_chr${i}.map | sed 's/:.*//g' > Mix_2_mappin_${i}_Lat; done
for i in {1..1}; do paste Mix_2_mappin_${i}_Lat Mix_2_Lat2.2_${i} | sed 's/\t/ /g' > Mix_2_Lat3_${i}; done

# chamei de 3
for i in {1..1}; do sed 's/1/0/g' Mix.RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_Ref1kgYRI_${i}.2.Viterbi.txt > Mix_3_Lat2_${i}; done
for i in {1..1}; do sed 's/2/1/g' Mix_3_Lat2_${i} > Mix_3_Lat2.1_${i}; done
for i in {1..1}; do sed 's/3/2/g' Mix_3_Lat2.1_${i} > Mix_3_Lat2.2_${i}; done
for i in {1..1}; do awk '{print $1, $3}'  Mix.admixsimu_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_Ref1kgYRI_chr${i}.map | sed 's/:.*//g' > Mix_3_mappin_${i}_Lat; done
for i in {1..1}; do paste Mix_3_mappin_${i}_Lat Mix_3_Lat2.2_${i} | sed 's/\t/ /g' > Mix_3_Lat3_${i}; done

# Baiano
#
for i in {1..1}; do sed 's/1/0/g' Baiano.RFmixed_Ref1kgYRI_Ref1kgIBS_${i}.2.Viterbi.txt > Baiano_Lat2_${i}; done
for i in {1..1}; do sed 's/2/1/g' Baiano_Lat2_${i} > Baiano_Lat2.1_${i}; done
for i in {1..1}; do sed 's/3/2/g' Baiano_Lat2.1_${i} > Baiano_Lat2.2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' Baiano.admixsimu_RFmixed_Ref1kgYRI_Ref1kgIBS_chr${i}.map | sed 's/:.*//g' > Baiano_mappin_${i}_Lat; done
for i in {1..1}; do paste Baiano_mappin_${i}_Lat Baiano_Lat2.2_${i} | sed 's/\t/ /g' > Baiano_Lat3_${i}; done

#PEL_IBS
export AdmixPOP=MAM
for i in {1..1}; do sed 's/2/0/g' MAM_RFmixed_Ref1kgPEL_Ref1kgIBS_1.2.Viterbi.txt > $AdmixPOP.2_Lat2_${i}; done
for i in {1..1}; do awk '{print $1, $3}'  MAMadmixsimu_RFmixed_Ref1kgPEL_Ref1kgIBS_chr1.map | sed 's/:.*//g' > $AdmixPOP.2_mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.2_mappin_${i}_Lat $AdmixPOP.2_Lat2_${i} | sed 's/\t/ /g' > $AdmixPOP.2_Lat3_${i}; done

export AdmixPOP=PEL
for i in {1..1}; do sed 's/2/0/g' PEL_RFmixed_Ref1kgPEL_Ref1kgIBS_1.2.Viterbi.txt > $AdmixPOP.2_Lat2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' PELadmixsimu_RFmixed_Ref1kgPEL_Ref1kgIBS_chr1.map | sed 's/:.*//g' > $AdmixPOP.2_mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.2_mappin_${i}_Lat $AdmixPOP.2_Lat2_${i} | sed 's/\t/ /g' > $AdmixPOP.2_Lat3_${i}; done


export AdmixPOP=exNATadmix
for i in {1..1}; do sed 's/2/0/g' exNAT_RFmixed_Ref1kgPEL_Ref1kgIBS_1.2.Viterbi.txt > $AdmixPOP.2_Lat2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' exNATadmixsimu_RFmixed_Ref1kgPEL_Ref1kgIBS_chr1.map | sed 's/:.*//g' > $AdmixPOP.2_mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.2_mappin_${i}_Lat $AdmixPOP.2_Lat2_${i} | sed 's/\t/ /g' > $AdmixPOP.2_Lat3_${i}; done

export AdmixPOP=exEURadmix
for i in {1..1}; do sed 's/2/0/g' exEUR_RFmixed_Ref1kgPEL_Ref1kgIBS_1.2.Viterbi.txt > $AdmixPOP.2_Lat2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' exEURadmixsimu_RFmixed_Ref1kgPEL_Ref1kgIBS_chr1.map | sed 's/:.*//g' > $AdmixPOP.2_mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.2_mappin_${i}_Lat $AdmixPOP.2_Lat2_${i} | sed 's/\t/ /g' > $AdmixPOP.2_Lat3_${i}; done

export AdmixPOP=Brasa
for i in {1..1}; do sed 's/1/0/g'  Brasa_RFmixed_Ref1kgPEL_Ref1kgIBS_Ref1kgYRI1.2.Viterbi.txt > Brasa_2_Lat2_${i}; done
for i in {1..1}; do sed 's/2/1/g' Brasa_2_Lat2_${i} > Brasa_2_Lat2.1_${i}; done
for i in {1..1}; do sed 's/3/2/g' Brasa_2_Lat2.1_${i} > Brasa_2_Lat2.2_${i}; done
for i in {1..1}; do awk '{print $1, $3}'  Brasaadmixsimu_RFmixed_Ref1kgPEL_Ref1kgIBS_Ref1kgYRI_chr1.map | sed 's/:.*//g' > Brasa_2_mappin_${i}_Lat; done
for i in {1..1}; do paste Brasa_2_mappin_${i}_Lat Brasa_2_Lat2.2_${i} | sed 's/\t/ /g' > Brasa_2_Lat3_${i}; done

#PEL_EAS_IBS
export AdmixPOP=MAM
for i in {1..1}; do sed 's/2/0/g' MAM_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_1.2.Viterbi.txt > $AdmixPOP.3_Lat2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' MAMadmixsimu_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_chr1.map | sed 's/:.*//g' > $AdmixPOP.3_mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.3_mappin_${i}_Lat $AdmixPOP.3_Lat2_${i} | sed 's/\t/ /g' > $AdmixPOP.3_Lat3_${i}; done

export AdmixPOP=PEL
for i in {1..1}; do sed 's/2/0/g' PEL_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_1.2.Viterbi.txt > $AdmixPOP.3_Lat2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' PELadmixsimu_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_chr1.map | sed 's/:.*//g' > $AdmixPOP.3_mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.3_mappin_${i}_Lat $AdmixPOP.3_Lat2_${i} | sed 's/\t/ /g' > $AdmixPOP.3_Lat3_${i}; done

export AdmixPOP=exNATadmix
for i in {1..1}; do sed 's/2/0/g'  exNAT_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_1.2.Viterbi.txt > $AdmixPOP.3_Lat2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' exNATadmixsimu_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_chr1.map | sed 's/:.*//g' > $AdmixPOP.3_mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.3_mappin_${i}_Lat $AdmixPOP.3_Lat2_${i} | sed 's/\t/ /g' > $AdmixPOP.3_Lat3_${i}; done

export AdmixPOP=exEURadmix
for i in {1..1}; do sed 's/2/0/g'  exEUR_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_1.2.Viterbi.txt > $AdmixPOP.3_Lat2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' exEURadmixsimu_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_chr1.map | sed 's/:.*//g' > $AdmixPOP.3_mappin_${i}_Lat; done
for i in {1..1}; do paste $AdmixPOP.3_mappin_${i}_Lat $AdmixPOP.3_Lat2_${i} | sed 's/\t/ /g' > $AdmixPOP.3_Lat3_${i}; done

export AdmixPOP=Brasa
for i in {1..1}; do sed 's/1/0/g' Brasa_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_Ref1kgYRI_1.2.Viterbi.txt > Brasa_3_Lat2_${i}; done
for i in {1..1}; do sed 's/2/1/g' Brasa_3_Lat2_${i} > Brasa_3_Lat2.1_${i}; done
for i in {1..1}; do sed 's/3/2/g' Brasa_3_Lat2.1_${i} > Brasa_3_Lat2.2_${i}; done
for i in {1..1}; do awk '{print $1, $3}' Brasaadmixsimu_RFmixed_Ref1kgPEL_EAS_Ref1kgIBS_Ref1kgYRI_chr1.map | sed 's/:.*//g' > Brasa_3_mappin_${i}_Lat; done
for i in {1..1}; do paste Brasa_3_mappin_${i}_Lat Brasa_3_Lat2.2_${i} | sed 's/\t/ /g' > Brasa_3_Lat3_${i}; done