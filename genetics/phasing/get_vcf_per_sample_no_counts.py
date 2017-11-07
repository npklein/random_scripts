import glob
import subprocess
import re
import os
import sys

cwd = os.getcwd()
print("WARNING: simple script using hard-coded paths")

with open('/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/samples.txt') as input_file:
    samples = input_file.read().split('\n')

for f in glob.glob('haplotypic_counts/*'):
    command = subprocess.Popen(["wc","-l",f], stdout=subprocess.PIPE)
    lines = command.communicate()[0].decode('ASCII').split(' ')[0]
    if int(lines) == 1:
        shapeitVCF_base = '/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/shapeitVCFwithAC/'
        chunk = re.match('.*(chr\d+\.\d+\.\d+).*', f).group(1)
        chr = chunk.split('.')[0]
        shapeitVCF = shapeitVCF_base+chr+'/genotypes_BIOS_LLDeep_noRNAedit_merged.'+chunk+'.shapeit.phased.withAC.vcf.gz'
        outVCF = cwd+'/vcf_per_sample_no_counts/BIOS_LLDeep_noRNAeditSites_phASER.'+sample+'.'+chunk+'.vcf.gz'
        sample = re.match('.*phASER\.(.*?).chr.*',f).group(1)
        bcftoolsProcess = subprocess.Popen(['bcftools','view','--samples',sample,'-Oz','-o',outVCF, shapeitVCF])
        bcftoolsProcess.wait()
        subprocess.Popen(['tabix',outVCF])
        print(shapeitVCF+' -> '+outVCF)
        sys.stdout.flush()
