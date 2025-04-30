#!/bin/bash -l

#SBATCH --job-name=gemma_second_pass
#SBATCH -D /group/jrigrp11/juliagh/GEMMA
#SBATCH -o /home/jgharenc/Raphanus_scripts/slurm_log/gemma_second_pass.%A_%a.out
#SBATCH -e /home/jgharenc/Raphanus_scripts/slurm_log/gemma_second_pass.%A_%a.err
#SBATCH -A jrigrp
#SBATCH -p high2
#SBATCH --time=10-00:00:00 #arbitratily large, no idea how long it will take
#SBATCH -c 2
#SBATCH --mem=150G

# load module
module load gemma/0.98.5

# set vairables
PLINK_PREFIX=./input/plink_filtered_58.2_138.1_rmd_0.2geno_0.5mind
KINSHIM_M=raph_gemma_kinship_0.2geno
MISS=1 # set to 1 so none are filtered since I am using plink to do the desired filtering
# genetate the kinship matrix
gemma -bfile "$PLINK_PREFIX" -gk 1 -miss $MISS -o "${KINSHIM_M}"
# -gk has options 1 and 2, 1 is centered (default) and 2 is standardized
#echo "done running kinship generation line"

# run gemma for anthocyanins
gemma -bfile "$PLINK_PREFIX" -k ./output/${KINSHIM_M}.cXX.txt -miss $MISS -lm 1 -n 2 -o anthocyanins_0.2geno_gemma_output
# -lm specifies which linear model; 1 is the reccomended Wald test
# -miss specifies a missingness threshold; using 0.1 now to allow 10% missingness, the default is 0.05 (in plink I already filtered for missingness of 0.1)
# -maf specifies minor allele freq threshold (default is 0.01; here not specifying bc plink filter I used is MORE stringent at 0.05 bc low coverage data)
# -r2 specifies r-squred threshold; any SNPs with r2 correlation with any of the covariates above the given value will not be included; default is 0.9999 but doesn't matter here as no covarites in this first pass

# run gemma for carotenoids 
gemma -bfile "$PLINK_PREFIX" -k ./output/${KINSHIM_M}.cXX.txt -miss $MISS -lm 1 -n 3 -o carotenoids_0.2geno_gemma_output

# run gemma for Glucobrassicin
gemma -bfile "$PLINK_PREFIX" -k ./output/${KINSHIM_M}.cXX.txt -miss $MISS -lm 1 -n 10 -o Glucobrassicin_RR_0.2geno_gemma_output

# run gemma for Glucoerucin
gemma -bfile "$PLINK_PREFIX" -k ./output/${KINSHIM_M}.cXX.txt -miss $MISS -lm 1 -n 8 -o Glucoerucin_RR_0.2geno_gemma_output

# run gemma for Glucoraphanin
gemma -bfile "$PLINK_PREFIX" -k ./output/${KINSHIM_M}.cXX.txt -miss $MISS -lm 1 -n 5 -o Glucoraphanin_RR_0.2geno_gemma_output

# run gemma for Glucoraphenin
gemma -bfile "$PLINK_PREFIX" -k ./output/${KINSHIM_M}.cXX.txt -miss $MISS -lm 1 -n 6 -o Glucoraphenin_RR_0.2geno_gemma_output

# run gemma for 4.Methylthio.3.butenyl.glucosinolate
gemma -bfile "$PLINK_PREFIX" -k ./output/${KINSHIM_M}.cXX.txt -miss $MISS -lm 1 -n 9 -o X4.Methylthio.3.butenyl.glucosinolate_RR_0.2geno_gemma_output

# run gemma for unknown1
gemma -bfile "$PLINK_PREFIX" -k ./output/${KINSHIM_M}.cXX.txt -miss $MISS -lm 1 -n 4 -o unknown1_RR_0.2geno_gemma_output

# run gemma for unknown2
gemma -bfile "$PLINK_PREFIX" -k ./output/${KINSHIM_M}.cXX.txt -miss $MISS -lm 1 -n 7 -o unknown2_RR_0.2geno_gemma_output

# run gemma for GS_OX_function
gemma -bfile "$PLINK_PREFIX" -k ./output/${KINSHIM_M}.cXX.txt -miss $MISS -lm 1 -n 11 -o GS_OX_function_RR_0.2geno_gemma_output

# run gemma for DOG_function
gemma -bfile "$PLINK_PREFIX" -k ./output/${KINSHIM_M}.cXX.txt -miss $MISS -lm 1 -n 12 -o DOG_function_RR_0.2geno_gemma_output

# run gemma for prop_indolic
gemma -bfile "$PLINK_PREFIX" -k ./output/${KINSHIM_M}.cXX.txt -miss $MISS -lm 1 -n 13 -o prop_indolic_RR_0.2geno_gemma_output
