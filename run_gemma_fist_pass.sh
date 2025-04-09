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

# genetate the kinship matrix
gemma -bfile ./input/plink_filtered_58.2_138.1_rmd_0.1geno_0.5mind -gk 1 -o raph_gemma_kinship
# -gk has options 1 and 2, 1 is centered (default) and 2 is standardized
#echo "done running kinship generation line"

# run gemma for anthocyanins
gemma -bfile ./input/plink_filtered_58.2_138.1_rmd_0.1geno_0.5mind -k ./output/raph_gemma_kinship.cXX.txt -miss 0.1 -lm 1 -o anthocyanins_raph1_gemma_output
# -lm specifies which linear model; 1 is the reccomended Wald test
# -miss specifies a missingness threshold; using 0.1 now to allow 10% missingness, the default is 0.05 (in plink I already filtered for missingness of 0.1)
# -maf specifies minor allele freq threshold (default is 0.01; here not specifying bc plink filter I used is MORE stringent at 0.05 bc low coverage data)
# -r2 specifies r-squred threshold; any SNPs with r2 correlation with any of the covariates above the given value will not be included; default is 0.9999 but doesn't matter here as no covarites in this first pass

# run gemma for carotenoids 
gemma -bfile ./input/plink_filtered_58.2_138.1_rmd_0.1geno_0.5mind -k ./output/raph_gemma_kinship.cXX.txt -miss 0.1 -lm 1 -n 2 -o carotenoids_raph1_gemma_output

# run gemma for Glucobrassicin_treatment_1 
gemma -bfile ./input/plink_filtered_58.2_138.1_rmd_0.1geno_0.5mind -k ./output/raph_gemma_kinship.cXX.txt -miss 0.1 -lm 1 -n 9 -o Glucobrassicin_treatment_1_raph1_gemma_output

# run gemma for Glucoerucin_treatment_1
gemma -bfile ./input/plink_filtered_58.2_138.1_rmd_0.1geno_0.5mind -k ./output/raph_gemma_kinship.cXX.txt -miss 0.1 -lm 1 -n 7 -o Glucoerucin_treatment_1_raph1_gemma_output

# run gemma for Glucoraphanin_treatment_1
gemma -bfile ./input/plink_filtered_58.2_138.1_rmd_0.1geno_0.5mind -k ./output/raph_gemma_kinship.cXX.txt -miss 0.1 -lm 1 -n 4 -o Glucoraphanin_treatment_1_raph1_gemma_output

# run gemma for Glucoraphenin_treatment_1
gemma -bfile ./input/plink_filtered_58.2_138.1_rmd_0.1geno_0.5mind -k ./output/raph_gemma_kinship.cXX.txt -miss 0.1 -lm 1 -n 5 -o Glucoraphenin_treatment_1_raph1_gemma_output

# run gemma for 4.Methylthio.3.butenyl.glucosinolate_trea
gemma -bfile ./input/plink_filtered_58.2_138.1_rmd_0.1geno_0.5mind -k ./output/raph_gemma_kinship.cXX.txt -miss 0.1 -lm 1 -n 8 -o X4.Methylthio.3.butenyl.glucosinolate_treatment_1_raph1_gemma_output

# run gemma for unknown1_treatment_1 
gemma -bfile ./input/plink_filtered_58.2_138.1_rmd_0.1geno_0.5mind -k ./output/raph_gemma_kinship.cXX.txt -miss 0.1 -lm 1 -n 3 -o unknown1_treatment_1_raph1_gemma_output

# run gemma fori unknown2_treatment_1
gemma -bfile ./input/plink_filtered_58.2_138.1_rmd_0.1geno_0.5mind -k ./output/raph_gemma_kinship.cXX.txt -miss 0.1 -lm 1 -n 6 -o unknown2_treatment_1_raph1_gemma_output
