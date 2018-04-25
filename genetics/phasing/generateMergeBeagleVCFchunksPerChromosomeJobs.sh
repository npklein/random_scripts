#!/bin/bash



#/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing
ROOT="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/readbackedPhasing/"

#For every chromosome do
for CHR in {1..22}
do

JOBDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/jobs/mergeBeagleVCFchunksPerChromosome/"

echo "Processing chr$CHR .."
cd $ROOT/chr$CHR/


echo "#!/bin/bash


module load GATK/3.7-Java-1.8.0_74
module list



java -Xmx18g -XX:ParallelGCThreads=1 -cp \$EBROOTGATK/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \\
-R /apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \\" >>$JOBDIR/mergeBeagleVCFchunks.chr$CHR.sh

for CHUNK in $( ls . | sort -k1n )
do

echo "-V /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/beagleChunksNoRNAedit/chr$CHR/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.$CHUNK.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.probs.gg.vcf.gz \\" >>$JOBDIR/mergeBeagleVCFchunks.chr$CHR.sh

done

cd -

echo "-out /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/beagleChunksNoRNAedit/genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.probs.gg.vcf.gz \\
-assumeSorted

cd /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/beagleChunksNoRNAedit/
md5sum genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.probs.gg.vcf.gz > genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.probs.gg.vcf.gz.md5
md5sum genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.probs.gg.vcf.gz.tbi > genotypes_BIOSfreeze2.1_LLDeep_Diagnostics_merged.chr$CHR.filter.GQ20_callRate50.BiallelicSNVsOnly.noRnaEditSites.probs.gg.vcf.gz.tbi.md5
cd -


" >>$JOBDIR/mergeBeagleVCFchunks.chr$CHR.sh

done