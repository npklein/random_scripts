#!/bin/bash



#/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing
ROOT="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/readbackedPhasing/"


##For each chromosome check for every chunk if all samples are existing, and if they have been readbackedphased or not

#For every chromosome do
for CHR in {10..10}
do

JOBDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/jobs/phasing/mergeReadbackPhasingSampleVCFsPerChunk/chr$CHR/"

mkdir -p $JOBDIR

echo "Processing chr$CHR .."
cd $ROOT/chr$CHR

#For every chunk do
for CHUNK in $( ls . )
do

echo "Processing chunk $CHUNK .."

echo "#!/bin/bash
#SBATCH --job-name=combineVariants.chr$CHR.$CHUNK
#SBATCH --output=combineVariants.chr$CHR.$CHUNK.out
#SBATCH --error=combineVariants.chr$CHR.$CHUNK.err
#SBATCH --time=5:59:59
#SBATCH --cpus-per-task 2
#SBATCH --mem 20gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


ENVIRONMENT_DIR=\".\"
set -e
set -u

echo \"## \"\$(date)\" Start \$0\"




module load GATK/3.7-Java-1.8.0_74
module list


java -Xmx18g -XX:ParallelGCThreads=1 -Djava.io.tmpdir=\${TMPDIR} -jar \$EBROOTGATK/GenomeAnalysisTK.jar \\
-T CombineVariants \\
-R /apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta \\" >$JOBDIR/combineVariants.chr$CHR.$CHUNK.sh






while read line
do


SAM=$( echo $line | awk '{ print $1 }' )

if [ -f $ROOT/chr$CHR/$CHUNK/vcf_per_sample/BIOS_LLDeep_noRNAeditSites_phASER.$SAM.chr$CHR.$CHUNK.vcf.gz ]
then
    OUT="$ROOT/chr$CHR/$CHUNK/vcf_per_sample/BIOS_LLDeep_noRNAeditSites_phASER.$SAM.chr$CHR.$CHUNK.vcf.gz";

elif [ -f $ROOT/chr$CHR/$CHUNK/vcf_per_sample_no_counts/BIOS_LLDeep_noRNAeditSites_phASER.$SAM.chr$CHR.$CHUNK.vcf.gz ]
then
    OUT="$ROOT/chr$CHR/$CHUNK/vcf_per_sample_no_counts/BIOS_LLDeep_noRNAeditSites_phASER.$SAM.chr$CHR.$CHUNK.vcf.gz";

else

echo "Chr$CHR, Chunk $CHUNK, for sample $SAM, cannot be found!";
exit 1

fi

echo "--variant $OUT \\" >> $JOBDIR/combineVariants.chr$CHR.$CHUNK.sh

done</groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/samples.txt
#/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/freeze2_complete_GTE_Groningen_07092016.txt

echo "-o $ROOT/chr$CHR/$CHUNK/BIOS_LLDeep_noRNAeditSites_phASER.chr$CHR.$CHUNK.vcf.gz \\
-genotypeMergeOptions REQUIRE_UNIQUE

cd $ROOT/chr$CHR/$CHUNK/
md5sum $ROOT/chr$CHR/$CHUNK/BIOS_LLDeep_noRNAeditSites_phASER.chr$CHR.$CHUNK.vcf.gz > $ROOT/chr$CHR/$CHUNK/BIOS_LLDeep_noRNAeditSites_phASER.chr$CHR.$CHUNK.vcf.gz.md5
md5sum $ROOT/chr$CHR/$CHUNK/BIOS_LLDeep_noRNAeditSites_phASER.chr$CHR.$CHUNK.vcf.gz.tbi > $ROOT/chr$CHR/$CHUNK/BIOS_LLDeep_noRNAeditSites_phASER.chr$CHR.$CHUNK.vcf.gz.tbi.md5
cd -

echo \"## \"\$(date)\" ##  \$0 Done \"

" >> $JOBDIR/combineVariants.chr$CHR.$CHUNK.sh


done

cd -

done
