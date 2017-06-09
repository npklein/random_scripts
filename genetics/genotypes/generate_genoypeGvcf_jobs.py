import math
import glob  

##### CHANGE BELOW VARIABLES #####
outdir="path/to/write/output/"
jobs_dir = "path/to/write/jobs/"
merged_vcf_dir = "directory/with/merged/vcf/files"
ref_genome = "/apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.fasta"
batches = 19
dbsnp = '/apps/data/dbSNP/dbsnp_138.b37.vcf'
project_name = 'BIOS_freeze2'

template = """#!/bin/bash
#SBATCH --job-name=GenotypeGvcf_chrREPLACECHROMOSOME
#SBATCH --output=GenotypeGvcf_chrREPLACECHROMOSOME.out
#SBATCH --error=GenotypeGvcf_chrREPLACECHROMOSOME.err
#SBATCH --qos=leftover
#SBATCH --time=3-23:59:59
#SBATCH --cpus-per-task 3
#SBATCH --mem 40gb
#SBATCH --nodes 1
#SBATCH --export=NONE
#SBATCH --get-user-env=L


ENVIRONMENT_DIR="."
set -e
set -u

echo "## "$(date)" Start $0"

#Generate input files, according to number of batches
for i in {0..REPLACENUMBEROFBATCHES}
do
    echo "getFile file=REPLACEINPREFIX.batch${i}_REPLACECHROMOSOME.g.vcf.gz"
    inputs+=" --variant REPLACEINPREFIX.batch${i}_chrREPLACECHROMOSOME.g.vcf.gz"
done


#Load gatk module
module load GATK/3.4-0-Java-1.7.0_80
module list

mkdir -p REPLACEOUTDIR

if java -Xmx38g -XX:ParallelGCThreads=2 -Djava.io.tmpdir=${TMPDIR} -jar $EBROOTGATK/GenomeAnalysisTK.jar \
 -T GenotypeGVCFs \
 -R REPLACEREFGENOME \
 --dbsnp REPLACEDBSNP \
 -o REPLACEOUTDIR/REPLACEOUTPUT \
 $inputs \
 -stand_call_conf 10.0 \
 -stand_emit_conf 20.0 \
 -L /apps/data/ftp.broadinstitute.org/bundle/2.8/b37/human_g1k_v37.chrREPLACECHROMOSOME.interval_list
then
 echo "returncode: $?"; 
 
cd REPLACEOUTDIR
md5sum REPLACEOUTPUT > REPLACEOUTPUT.md5
cd -
 echo "succes moving files";
else
 echo "returncode: $?";
 echo "fail";
fi


"""

chromosomes = ['1','2','3','4','5','6','7',
        '8','9','10','11','12','13','14',
         '15','16','17','18','19','20',
         '21','22','23','24','25']

for chr in chromosomes:
    outfile = jobs_dir+'GenotypeGvcf_chr'+chr+'.sh'
    with open(outfile,'w') as out:
        new_template = template.replace('REPLACENUMBEROFBATCHES',str(batches))
        new_template = new_template.replace('REPLACECHROMOSOME',chr)
        new_template = new_template.replace('REPLACEOUTDIR',outdir)
        new_template = new_template.replace('REPLACEOUTPUT',project_name+'_chr'+chr+'.gg.vcf.gz')
        new_template = new_template.replace('REPLACEINPREFIX', merged_vcf_dir+project_name)
        new_template = new_template.replace('REPLACEREFGENOME',ref_genome)
        new_template=new_template.replace('REPLACEDBSNP', dbsnp)
        out.write(new_template)
    print(outfile)

