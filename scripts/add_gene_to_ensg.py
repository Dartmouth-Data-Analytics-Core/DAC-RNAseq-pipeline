import sys

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

counts_file = open(sys.argv[2],'r')

header = counts_file.readline().strip('\n').split('\t')
print ('\t'.join(['Ensembl ID', 'Gene Name'] + header[1:]))

for line in counts_file:
    sline = line.strip('\n').split('\t')
    ensg_name = sline[0]
    to_gene = ensg_to_gene[ensg_name]

    print ('\t'.join([ensg_name, to_gene] + sline[1:]))









