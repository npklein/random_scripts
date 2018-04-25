import glob
import subprocess
import re
import os
import sys

cwd = os.getcwd()
print("WARNING: simple script using hard-coded paths")

with open('/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/samples.txt') as input_file:
    samples = input_file.read().split('\n')
chr = cwd.split('/')[-2]
chunk = chr+'.'+cwd.split('/')[-1]

haplotypic_files = glob.glob('haplotypic_counts/*')
if len(haplotypic_files) > 0:
    for f in glob.glob('haplotypic_counts/*'):
        sample = re.match('.*phASER\.(.*?)\.chr.*',f).group(1)
        samples.remove(sample)
#    command = subprocess.Popen(["wc","-l",f], stdout=subprocess.PIPE)
#    lines = command.communicate()[0].decode('ASCII').split(' ')[0]

    if not os.path.exists(cwd+'/vcf_per_sample_no_counts/'):
        os.mkdir(cwd+'/vcf_per_sample_no_counts/')

samples = filter(None, samples) # fastest
for sample in samples:    
    shapeitVCF_base = '/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/shapeitVCFwithAC/'
    shapeitVCF = shapeitVCF_base+chr+'/genotypes_BIOS_LLDeep_noRNAedit_merged.'+chunk+'.shapeit.phased.withAC.vcf.gz'
    outVCF = cwd+'/vcf_per_sample_no_counts/BIOS_LLDeep_noRNAeditSites_phASER.'+sample+'.'+chunk+'.vcf.gz'
    print(shapeitVCF+' -> '+outVCF)
    sys.stdout.flush()

    bcftoolsProcess = subprocess.Popen(['bcftools','view','--samples',sample,'-Oz','-o',outVCF, shapeitVCF])
    bcftoolsProcess.wait()
    p = subprocess.Popen(['tabix',outVCF])
    p.wait()
