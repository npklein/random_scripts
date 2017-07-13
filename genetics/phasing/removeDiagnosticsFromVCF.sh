module load BCFtools
INPUTDIR=/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/results/filterGQ20_callRate50/

for vcf in ${INPUTDIR}/*filter.GQ20_callRate50.PASSonly.BiallelicSNVsOnly.gg.vcf.gz
do
  echo $vcf;
  bcftools view \
  --samples-file ^/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/diagnosticSamples.txt \
  -o ${vcf%.gg.vc.fz}.noDiagnostics.gg.vcf.gz \
  -O z $vcf
done
