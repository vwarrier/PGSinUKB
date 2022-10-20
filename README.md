# PGSinUKB

This is a script to generate a PGS file in the UKB

Please follow the UKB file to download the data and do some basic QC. We will generate PGS using PRSice2, save the polygenic scores 
and run a seperate regression as and when needed. 


## Documents to read
Please read the following documents:

1. Genotype data: https://biobank.ctsu.ox.ac.uk/crystal/crystal/docs/ukbgene_instruct.html
2. WES data: https://www.ukbiobank.ac.uk/wp-content/uploads/2019/08/UKB-50k-Exome-Sequencing-Data-Release-July-2019-FAQs.pdf

I have currently downloaded only the joint called files for WES. Why is joint call better than sample-level calling? https://gatkforums.broadinstitute.org/gatk/discussion/3686/why-do-joint-calling-rather-than-single-sample-calling-retired

Please note, the variant calls are not good at all for the downloaded UKB data. 


## Step 1: Dowloading files from the UKB

Have a look at this to understand genomic files available: http://biobank.ctsu.ox.ac.uk/crystal/label.cgi?id=100314

Ensure the key is downloaded in the file and saved as .ukbkey

The scripts below help you download the genetics file

```bash
## genotype data
for i in {1..22}; do ./ukbgene imp -c$i; done  #this downloads the BGEN files
for i in {1..22}; do ./ukbgene imp -m -c$i; done

./ukbgene imp -cX
./ukbgene imp -cXY

./ukbgene imp -m -cX
./ukbgene imp -m -cXY

##whole exome data
for i in {1..22}; do ./ukbgene efe -c$i; done
```

To download the phenotypic file, please see the email you have, and download the file and save it in the R format. 


## Step 2: Converting BGEN to plink files

To run this command, let's first identify people who belong to British European subset 22006 (https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=22006), subset it to that and then run all the commands (both maf and hwe are sensitive to genetic grouping).

```bash
for i in {1..21}; do ./plink2 --bgen ukb_imp_chr${i}_v3.bgen --sample ukb20904_imp_chr${i}_v3_s487334.sample --keep europeansonlyfile --make-bed -out ukbchr${i} --maf 0.01 --geno 0.05 --threads 10 --hwe 0.000001 --mind 0.05; done

```


## Step 3: Create a merged bfile with all autosomes

```bash

./plink --bfile ukbchr1 --merge-list bmergfile.txt --extract snpextractfile.txt --make-bed --out UKB2_autosomes

for i in {1..22}; do ./plink --bfile ukbchr${i} --exclude UKB2_autosomes-merge.missnp --make-bed --out ukbchr_v2_${i}; done

./plink --bfile ukbchr_v2_1 --merge-list bmerg2.txt --make-bed --out UKB2_autosomes

./plink --bfile UKB2_autosomes --extract snpextractfile.txt --make-bed --out UKB2_prsicefile

```


## Step 4: Generate PGS

```bash
#casecontrol gwas
for i in autismdanerforprsice sczprsice; do Rscript PRSice2.R --dir . \
--prsice ./PRSice2_linux --base ./pgssumstats/${i}.txt \
--target ~/UKB_v2/Plink_files/UKB2_autosomes \
--thread 20 \
--stat OR \
--binary-target T \
--out ./PRSice2results/UKB2_autosomes_${i} \
--bar-levels 0.0001,0.001,0.01,0.1,0.25,0.5,0.75,1 \
--no-regress \
--fastscore; done

#quantitative trait gwas
for i in eduyears SQprsice cognitionprsice rmetprsice empathyprsice friendshipmtagprsice familymtagprsice EQ; \
do Rscript PRSice2.R \
--dir . \--prsice ./PRSice2_linux \
--base ./pgssumstats/${i}.txt \
--target ~/UKB_v2/Plink_files/UKB2_autosomes \
--thread 10 \
--stat BETA \
--binary-target T \
--out ./PRSice2results/UKB2_autosomes_${i} \
--bar-levels 0.0001,0.001,0.01,0.1,0.25,0.5,0.75,1 \
--no-regress \
--fastscore; done

```

## Step 5: Run regression analysis

I provide an example script below, all done in R. 

### Step 5a: Load the phenotypic data

```R
bd <- fread("~/UKB_v2/ukb20904.tab", header=TRUE, sep="\t")
#bd$f.53.0.0 <- as.Date(bd$f.53.0.0)
#bd$f.53.1.0 <- as.Date(bd$f.53.1.0)
#bd$f.53.2.0 <- as.Date(bd$f.53.2.0)
#bd$f.20400.0.0 <- as.Date(bd$f.20400.0.0) # date of cancer diagnosis
#bd$f.40005.0.0 <- as.Date(bd$f.40005.0.0)
#bd$f.40005.1.0 <- as.Date(bd$f.40005.1.0)
#bd$f.40005.2.0 <- as.Date(bd$f.40005.2.0)
#bd$f.40005.3.0 <- as.Date(bd$f.40005.3.0)
#bd$f.40005.4.0 <- as.Date(bd$f.40005.4.0)
#bd$f.40005.5.0 <- as.Date(bd$f.40005.5.0)
#bd$f.40005.6.0 <- as.Date(bd$f.40005.6.0)
#bd$f.40005.7.0 <- as.Date(bd$f.40005.7.0)
#bd$f.40005.8.0 <- as.Date(bd$f.40005.8.0)
#bd$f.40005.9.0 <- as.Date(bd$f.40005.9.0)
#bd$f.40005.10.0 <- as.Date(bd$f.40005.10.0)
#bd$f.40005.11.0 <- as.Date(bd$f.40005.11.0)
#bd$f.40005.12.0 <- as.Date(bd$f.40005.12.0)
#bd$f.40005.13.0 <- as.Date(bd$f.40005.13.0)
#bd$f.40005.14.0 <- as.Date(bd$f.40005.14.0)
#bd$f.40005.15.0 <- as.Date(bd$f.40005.15.0)
#bd$f.40005.16.0 <- as.Date(bd$f.40005.16.0)
#bd$f.40005.17.0 <- as.Date(bd$f.40005.17.0)
#bd$f.40005.18.0 <- as.Date(bd$f.40005.18.0)
#bd$f.40005.19.0 <- as.Date(bd$f.40005.19.0)
#bd$f.40005.20.0 <- as.Date(bd$f.40005.20.0)
#bd$f.40005.21.0 <- as.Date(bd$f.40005.21.0)
#bd$f.40005.22.0 <- as.Date(bd$f.40005.22.0)
#bd$f.40005.23.0 <- as.Date(bd$f.40005.23.0)
#bd$f.40005.24.0 <- as.Date(bd$f.40005.24.0)
#bd$f.40005.25.0 <- as.Date(bd$f.40005.25.0)
#bd$f.40005.26.0 <- as.Date(bd$f.40005.26.0)
#bd$f.40005.27.0 <- as.Date(bd$f.40005.27.0)
#bd$f.40005.28.0 <- as.Date(bd$f.40005.28.0)
#bd$f.40005.29.0 <- as.Date(bd$f.40005.29.0)
#bd$f.40005.30.0 <- as.Date(bd$f.40005.30.0)
#bd$f.40005.31.0 <- as.Date(bd$f.40005.31.0)
#bd$f.42006.0.0 <- as.Date(bd$f.42006.0.0) # first trike outcome
#bd$f.42008.0.0 <- as.Date(bd$f.42008.0.0) # ischemic stroke outcome
#bd$f.42010.0.0 <- as.Date(bd$f.42010.0.0) # introcerebral heamorrhage
#bd$f.42012.0.0 <- as.Date(bd$f.42012.0.0) # subarachnoid heamorrhage

#How to create a file for analysis
pheno  = bd[,c("f.eid", "yourpheno", "f.31.0.0", "f.22009.0.1",
               "f.22009.0.2", "f.22009.0.3", "f.22009.0.4", "f.22009.0.5",
               "f.22009.0.6", "f.22009.0.7", "f.22009.0.8", "f.22009.0.9",
               "f.22009.0.10", "f.22009.0.11", "f.22009.0.12", "f.22009.0.13",
               "f.22009.0.14", "f.22009.0.15", "f.22009.0.16", "f.22009.0.17",
               "f.22009.0.18", "f.22009.0.19", "f.22009.0.20", "f.22009.0.21",
               "f.22009.0.22", "f.22009.0.23", "f.22009.0.24", "f.22009.0.25",
               "f.22009.0.26", "f.22009.0.27", "f.22009.0.28", "f.22009.0.29",
               "f.22009.0.30", "f.22009.0.31", "f.22009.0.32", "f.22009.0.33",
               "f.22009.0.34", "f.22009.0.35", "f.22009.0.36", "f.22009.0.37",
               "f.22009.0.38", "f.22009.0.39", "f.22009.0.40", "f.34.0.0", 
               "f.22006.0.0","f.22000.0.0", "f.22007.0.0", "f.22027.0.0", 
               "f.22001.0.0", "f.22021.0.0")]

pheno2 = pheno[!is.na(pheno$yourpheno),] #remove items that are na in the pheno
pheno2 = pheno2[!is.na(pheno2$f.22006.0.0),] #remove non-europeans
pheno2$checksex = ifelse(pheno2$f.31.0.0 == pheno2$f.22001.0.0, "correct", "incorrect")
pheno2 = subset(pheno2, checksex == "correct") # remove sex mismatches
pheno2 = pheno2[is.na(pheno2$f.22027.0.0),] #remove excessive heterozygosity

#The dreaded removal of related individuals. This is a basic script, and removes more individuals than needed. 
rel = fread("~/UKB_v2/ukb20904_rel_s488302.dat")
rel2 = subset(rel, Kinship >  0.0884) #equivalent to 3rd degree relatives
alpha = count(rel2, vars = "ID1")
setnames(alpha, 2, "ID1_freq")
rel2 = merge(rel2, alpha)
one = subset(rel2, ID1_freq > 1)
pheno2 = pheno2[!(pheno2$f.eid %in% one$ID1),]
two = subset(rel2, ID1_freq < 2)
pheno2 = pheno2[!(pheno2$f.eid %in% two$ID2),]
```
