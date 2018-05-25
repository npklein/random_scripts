# ['probeID', 'BETA', 'SE', 'P_VAL', 'N_samp', 'mean', 'SD', 'pheno\n']
cpg_info = {}
for input_fileName in ['/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/7pheno_8cg_ewas_model1_cpg_technic_cov_10112017.txt',
                       '/groups/umcg-lld/tmp03/umcg-ndeklein/EWAS/output/7pheno_8cg_ewas_model2_cpg_technic_cov_age_sex_smoking_10112017.txt']:
    print(input_fileName)
    outfile = input_fileName.replace('.txt', '.resultsFile.txt')
    with open(input_fileName) as input_file:
        print(input_file.readline().split('\t'))
        for line in input_file:
            line = line.strip().split('\t')
            cpg = line[0]
            if cpg not in cpg_info:
                cpg_info[cpg] = {}
            phenotype = line[-1]
            estimate = line[1]
            se = line[2]
            pval = line[3]
            n = line[4]
            cpg_info[cpg][phenotype] = [estimate, se, pval, n]
        
    pheno_order = ['BMI','Triglycerides','Waist-circumference','HDL-C','Fasting Glucose','Systolic BP','DiastolicBP']
    pheno_map = {'BMI':'antrop_BMI',
                   'Triglycerides':'Biochem_TG_log',
                   'Waist-circumference':'waist_circumference',
                   'HDL-C':'Biochem_HDL',
                   'Fasting Glucose':'Biochem_Glucose_log',
                   'Systolic BP':'antrop_SBP',
                   'DiastolicBP':'antrop_DBP'}
    cpg_order = ['cg04535902', 'cg09662411', 'cg09935388', 'cg10399789', 'cg12876356','cg18146737','cg14179389','cg18316974']
    with open(outfile,'w') as out:
        out.write('CpG_ID\tMetabolic Phenotype\tEstimate\tStandard Error\tConfidence Interval\tP-value\tN')
        for cpg in cpg_order:
            for pheno in pheno_order:
                data = cpg_info[cpg][pheno_map[pheno]]
                ci_upper = str(float(data[0])+(1.96*float(data[1])))
                ci_lower = str(float(data[0])-(1.96*float(data[1])))
                out.write('\n'+cpg+'\t'+pheno+'\t'+data[0]+'\t'+data[1]+'\t['+ci_lower+','+ci_upper+']\t'+data[2]+'\t'+data[3])
            out.write('\n')


    print('written to '+outfile)
