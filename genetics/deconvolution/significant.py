pval = 0.05
with open('BaseCellPerc_merged_filtered.csv') as input_file, open('significant.txt','w') as out, open('directions_if_significant.txt','w') as out2, open('pvalue_if_significant.txt') out3:
    input_file.readline()
    out.write('neut\tlymph\tmono\teos\twhole_blood\n')
    x = 1
    for line in input_file:
        line = line.strip().split('\t')
        neut_pval = float(line[1])
        lymph_pval = float(line[2])
        mono_pval = float(line[3])
        eos_pval = float(line[4])
        whole_blood_pval = float(line[-1])
        neut_direction = line[13]
        lymph_direction = line[14]
        mono_direction = line[15]
        eos_direction = line[16]
        whole_blood_direction = line[-2]
        if neut_pval < pval:
            out.write(str(x))
            out2.write(neut_direction)
        else:
            out.write('NA')
            out2.write('NA')
        out.write('\t')
        out2.write('\t')
        if lymph_pval < pval:
            out.write(str(x))
            out2.write(lymph_direction)
        else:
            out.write('NA')
            out2.write('NA')
        out.write('\t')
        out2.write('\t')
        if mono_pval < pval:
            out.write(str(x))
            out2.write(mono_direction)
        else:
            out.write('NA')
            out2.write('NA')
        out.write('\t')
        out2.write('\t')
        if eos_pval < pval:
            out.write(str(x))
            out2.write(eos_direction)
        else:
            out.write('NA')
            out2.write('NA')
        out.write('\t')
        out2.write('\t')
        if whole_blood_pval < pval:
            out.write(str(x))
            out2.write(whole_blood_direction)
        else:
            out.write('NA')
            out2.write('NA')
        x += 1
        out.write('\n')
        out2.write('\n')
