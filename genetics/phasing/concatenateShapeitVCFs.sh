module load BCFtools
module load tabix
# Concat per chromosome
cd /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/results/shapeitVCF/
for i in chr*; 
do
   echo $i
   vcfs=$(ls -v $i/*.vcf.gz)
   bcftools concat $vcfs -O z -o BIOS_LLDeep_Diagnositcs.shapeitPhased.chr$i.vcf.gz
   tabix BIOS_LLDeep_Diagnositcs.shapeitPhased.chr$i.vcf.gz
done
