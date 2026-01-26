import sys
import pandas as pd

'''
#!genome-build GRCh38.p12
#!genome-version GRCh38
#!genome-date 2013-12
#!genome-build-accession NCBI:GCA_000001405.27
#!genebuild-last-updated 2019-03
1       havana  gene    11869   14409   .       +       .       gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; gene_biotype "transcribed_unprocessed_pseudogene";
'''

'''
Geneid  Chr     Start   End     Strand  Length  ai09166 aq04898 au02019 bc07966 bd09479 bi07249 bm09152 br02535 ca07040 cd01173 dv00351 eh02957 en09652 eq01118      fl04721 gc02452 gp06952 gz03095 ik09540 io05883 ir07471 ji04434 jm00238 kn05854 lh01569 lr04577 ly04595 ma02224 mc09212 na08228 OX02514 ox06988 pa09915      pl04872 qf04849 qo05187 qy04712 rb06824 ri02200 rw08819 sb03953 sf03295 sl00765 sz01638 tf06738 tp04076 tp04704 ul01308 up04513 vj05882 vl04509 wc08450      wf00978 wm05009 ws09552 yc09610 zg00030 zj05762 zp06755 zt00066
ENSG00000223972 1;1;1;1;1;1;1;1;1       11869;12010;12179;12613;12613;12975;13221;13221;13453   12227;12057;12227;12721;12697;13052;13374;14409;13670   +;+;+;+;+;+;+;+;+    1735    0       0       1       1       0       0       0       0       0       0       0  

##### Additional Note: 
##### ENSG00000145362.21 has enough different transcript isoforms that fields get malformed when output files are opened in excel. 
##### logic is included to this behaviour by replcing chr, start, end, strand values with only those of first isoform from GTF, and 
##### some additional text alerting reader to issue 
'''

gtf_file = open(sys.argv[1],'r')

ensg_to_gene = {}

for line in gtf_file:
    sline = line.strip('\n').split('\t')
    if len(sline) < 5:
        continue

    if sline[2] != 'gene':
        continue

    info = sline[8].split(';')

    for i in info:
        if 'gene_id' in i:
            ensg = i.split('"')[1]
        if 'gene_name' in i:
            gene = i.split('"')[1]

    ensg_to_gene[ensg] = gene

# read counts file 
counts_file_path = sys.argv[2]
counts_file = pd.read_csv(counts_file_path, sep='\t')

# extract chr, start, end, and strand for ENSG00000145362.21
chr = counts_file['Chr'][counts_file['Geneid']=='ENSG00000145362.21']
start = counts_file['Start'][counts_file['Geneid']=='ENSG00000145362.21']
end = counts_file['End'][counts_file['Geneid']=='ENSG00000145362.21']
strand = counts_file['Strand'][counts_file['Geneid']=='ENSG00000145362.21']

# additional text to add to replaced fields 
add_text = ";additional_isoform_details_removed_for_processing"
chr_1 = chr.str.split(';').tolist()[0][0] + add_text
start_1 = start.str.split(';').tolist()[0][0] + add_text
end_1 = end.str.split(';').tolist()[0][0] + add_text
strand_1 = strand.str.split(';').tolist()[0][0] + add_text

# replace chr, start, and end values for ENSG00000145362.21
counts_file.loc[counts_file['Geneid']=='ENSG00000145362.21', 'Chr'] = chr_1
counts_file.loc[counts_file['Geneid']=='ENSG00000145362.21', 'Start'] = start_1
counts_file.loc[counts_file['Geneid']=='ENSG00000145362.21', 'End'] = end_1
counts_file.loc[counts_file['Geneid']=='ENSG00000145362.21', 'Strand'] = strand_1

header = list(counts_file.columns)
print ('\t'.join(['Ensembl ID', 'Gene Name'] + header[1:]))

for idx, row in counts_file.iterrows():
    ensg_name = row['Geneid']
    to_gene = ensg_to_gene[ensg_name]
    print ('\t'.join([ensg_name, to_gene] + row[1:].astype(str).tolist()))









