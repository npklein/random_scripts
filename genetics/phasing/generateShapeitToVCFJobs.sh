phasedGeneChunks="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/phasedGeneChunks.21062017.csv"

for i in {1..22}
do

mkdir -p /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/jobs/convertShapeit/chr$i/

done

while read line
do
echo $line
CHR=`echo $line | awk '{print $1}' FS=","`
START=`echo $line | awk '{print $2}' FS=":" | awk '{print $1}' FS="-"`
END=`echo $line | awk '{print $2}' FS=":" | awk '{print $2}' FS="-"`
RESULTSDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/shapeitVCF/chr$CHR/"
INPUTDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/shapeit/chr$CHR/"

echo "#!/bin/bash
#SBATCH --job-name=chr$CHR.$START.$END.ConvertShapeitToVCF
#SBATCH --output=ConvertShapeitToVCF.chr$CHR.$START.$END.out
#SBATCH --error=ConvertShapeitToVCF.chr$CHR.$START.$END.err
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task 4
#SBATCH --mem 8gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=30L
#SBATCH --qos=leftover

set -e
set -u

ENVIRONMENT_DIR='.'


module load tabix/0.2.6-foss-2015b
module load shapeit/v2.r837-static
module list

echo \"## \"\$(date)\" Start \$0\"

mkdir -p ${RESULTSDIR}


zcat ${INPUTDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHR.$START.$END.shapeit.phased.hap.gz \\
        > ${INPUTDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHR.$START.$END.shapeit.phased.haps

if [ -f  ${INPUTDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHR.$START.$END.shapeit.phased.hap.gz.sample ];
then
    mv ${INPUTDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHR.$START.$END.shapeit.phased.hap.gz.sample \\
        ${INPUTDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHR.$START.$END.shapeit.phased.sample
fi

shapeit \\
 -convert \\
 --input-haps ${INPUTDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHR.$START.$END.shapeit.phased \\
 --output-vcf ${RESULTSDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHR.$START.$END.shapeit.phased.vcf.gz

 echo \"returncode: \$?\";
 cd ${RESULTSDIR}
 bname=\$(basename ${RESULTSDIR}/genotypes_BIOS_LLDeep_Diagnostics_merged.chr$CHR.$START.$END.shapeit.phased.vcf.gz)
 
 # has to be bgzipped
 gunzip \${bname}
 bgzip \${bname%.gz}
 tabix \${bname}
 echo \"making md5sums\"
 md5sum \${bname} > \${bname}.md5

 cd -
 echo \"succes moving files\";


echo \"## \"\$(date)\" ##  \$0 Done \"

">./jobs/convertShapeit/chr$CHR/ConvertShapeitToVCF.chr$CHR.$START.$END.sh

done<$phasedGeneChunks


