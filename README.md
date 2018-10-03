# PGSinUKB

This is a script to generate a PGS file in the UKB

Please follow the UKB file to download the data and do some basic QC. We will generate PGS using PRSice2, save the polygenic scores 
and run a seperate regression as and when needed. 


# Step 1: Create a merged bfile with all autosomes






## Generate PGS

```bash
#casecontrol gwas
for i in autismdanerforprsice sczprsice; do Rscript PRSice2.R --dir . \
--prsice ./PRSice2_linux --base ./pgssumstats/${i}.txt \
--target ~/SFARI/liftOverPlink/SSC_1Mv3/SSCimputationfile-updated \
--thread 10 \
--stat OR \
--binary-target T \
--out ./PRSice2results/SSC_1Mv3_${i} \
--bar-levels 0.0001,0.001,0.01,0.1,0.25,0.5,0.75,1 \
--no-regress \
--fastscore; done

#quantitative trait gwas
for i in eduyears SQprsice cognitionprsice rmetprsice empathyprsice friendshipmtagprsice familymtagprsice EQ; \
do Rscript PRSice2.R \
--dir . \--prsice ./PRSice2_linux \
--base ./pgssumstats/${i}.txt \
--target ~/SFARI/liftOverPlink/SSC_1Mv3/SSCimputationfile-updated \
--thread 10 \
--stat BETA \
--binary-target T \
--out ./PRSice2results/SSC_1Mv3_${i} \
--bar-levels 0.0001,0.001,0.01,0.1,0.25,0.5,0.75,1 \
--no-regress \
--fastscore; done

```
