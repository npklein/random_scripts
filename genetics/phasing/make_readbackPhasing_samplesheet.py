import os
seen = []
with open('sample_individual_idLink.txt') as input_file, open('individual_bam_link.txt','w') as out:
    out.write('individualID,sampleName,bam\n')
    input_file.readline()
    for line in input_file:
        sample_id = line.split('\t')[0]
        individual_id = line.strip().split('\t')[1]
        if individual_id in seen:
            continue
        seen.append(individual_id)
        path = '/groups/umcg-bios/tmp03/projects/masked_BAMs/diagnostic/mergeAndReheaderBAMs/results/'+individual_id+'.mdup.sorted.readGroupsAdded.bam'
        if os.path.exists(path):
            bam_file = path
        else:
            print(individual_id,' bam not found')
            print(path)
            continue
        out.write(individual_id+','+individual_id+','+bam_file+'\n')

converter = {}
with open('/groups/umcg-bios/tmp03/projects/bbmriSampleInfo/sampleSheetDirectlyFromMdb26-01-2016.txt') as lldeepConverter:
    for line in lldeepConverter:
        line = line.split('\t')
        id = line[0]
        runID = line[1]    
        converter[id] = runID  

#LL-LLDeep_0043    BD1NR9ACXX-4-19
with open('freeze2_complete_GTE_Groningen_07092016.txt') as input_file,open('individual_bam_link.txt','a') as out:
    input_file.readline()
    for  line in    input_file:
        line = line.strip().split('\t')
        path = '/groups/umcg-bios/tmp03/projects/masked_BAMs/BAMsReadGroupsAdded/'+line[1]+'.mdup.sorted.readGroupsAdded.bam'
        if os.path.exists(path):
            bam_file = path
        else:
            lldeep_id = '_'.join(line[0].split('_')[1:])
            newID = converter['LL-'+lldeep_id]
            line[0] = newID
            path = '/groups/umcg-bios/tmp03/projects/masked_BAMs/BAMsReadGroupsAdded/'+newID+'.mdup.sorted.readGroupsAdded.bam'

            if os.path.exists(path):
                bam_file = path
            else:
                print(path+' bam not found')
                continue
        out.write(line[0]+','+line[1]+','+bam_file+'\n')


with open('lldeepNotInBiosSamples.txt') as input_file, open('individual_bam_link.txt','a') as out:
    for line in input_file:
        line = line.strip().split('\t')
        sampleID = line[0]
        genotypeID = line[1]
        path = '/groups/umcg-bios/tmp03/projects/masked_BAMs/lldeepNotInBIOS/'+sampleID+'-lib1.bam'
        if os.path.exists(path):
            bam_file = path
        else:
            print(path+' bam not found')
            continue
        out.write(sampleID+','+genotypeID+','+bam_file+'\n')
