import glob

#for name in ['BaseCellPerc','Round1CellPerc','Round2CellPerc','Round3CellPerc']:
for name in ['BaseDeconindexed_new']:
    x = 0
    print(name)
    with open(name+'_merged.csv','w') as out:
        for csv in glob.glob('chunk_*/'+name+'/deconvolutionResults.csv'):
            print(csv)
            with open(csv) as input_file:
                header = input_file.readline()
                if x == 0:
                    out.write(header)
                else:
                    if header != prev_header:
                        raise RuntimeError('Headers not the same')
                for line in input_file:
                    out.write(line)
                prev_header = header        
            x += 1

