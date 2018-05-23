import gzip


exon_info = {}
prev_start = None
prev_stop = None
prev_chr = None
current_exon_ids = []
outfile = None

with gzip.open('/apps/data/ftp.ensembl.org/pub/release-75/bed/homo_sapiens/Homo_sapiens.GRCh37.75.bed.gz') as input_file:
    for line in input_file:
        line = line.decode('utf-8')
        line = line.strip().split('\t')
        if line[7] == 'exon':
            exon_id = line[9].split('exon_id "')[1].split('"')[0]
            start = line[1]
            chr = line[0]
            chr_range = [str(x) for x in range(1,23)]
            if chr not in chr_range:
                continue
            # reset the previous start/stop when new chr is found
            if not prev_chr or chr != prev_chr:
                if outfile:
                    outfile.close()
                print('Start chr'+chr)
                outfile = open('Homo_sapiens.GRCh37.75.chr'+chr+'.metaExons.bed','w')
                prev_start = None
                prev_stop = None
            # make sure that the bed file is sorted by start position
            if prev_start and int(start) < int(prev_start):
                print('start: '+start+' prev start: '+prev_start)
                print('chr: '+chr+ ' prev chr: '+prev_chr)
                raise RuntimeError('bed file should be sorted by start, but is not')
            stop = line[2]
            # Because bed file is ordered by start position, if the start position < previous stop position, the exons overlap. Merge
            if prev_stop and int(start) < int(prev_stop):
                start = prev_start
                # stop might not always be bigger than previous stop (if the exon completely within the previous exon), so check first
                if int(stop) < int(prev_stop):
                    stop = prev_stop
                if exon_id not in current_exon_ids:
                    current_exon_ids.append(exon_id)
            else:
                # if prev_stop == None, it is the first line and should not write to file yet
                if prev_stop:
                    outfile.write(chr+'\t'+prev_start+'\t'+prev_stop+'\t'+';'.join(current_exon_ids)+'\n')
                    current_exon_ids = [exon_id]
                else:
                    current_exon_ids.append(exon_id)
            prev_start = start
            prev_stop = stop
            prev_chr = chr
    outfile.close()

