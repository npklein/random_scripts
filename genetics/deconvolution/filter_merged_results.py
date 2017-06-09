name = 'BaseDeconindexed_new_merged.csv'
with open(name) as input_file, open(name+'_filtered.csv','w') as out:
    header = input_file.readline()
    out.write(header)
    pvalue_index = []
    elements_to_write = []
    for index, element in enumerate(header.split('\t')):
        if element.endswith('pvalue'):
            pvalue_index.append(index)
            elements_to_write.append(index)
        if element == 'Spearman correlation p-value':
            pvalue_index.append(index)
            elements_to_write.append(index)
        if element == 'Spearman correlation expression~GT':
            elements_to_write.append(index)
    for line in input_file:
        elements = line.strip().split('\t')
        add_row = False
        for index in pvalue_index:
            if float(elements[index]) < 0.05:
               add_row = True
        if add_row:
            out.write(line)
