import sys, os
from collections import namedtuple
from xlsxwriter.workbook import Workbook

annot_tsv = sys.argv[1]

EH = namedtuple('EH', 'pos, motif, gene, size, gt')
eh_gt = {}
with open(annot_tsv) as f:
    annot = {e:l for e,l in enumerate(f.readline().strip().split("\t")[1:]) }
    for i in f:
        i = i.strip().split("\t")
        sample = i[0]
        if not sample in eh_gt:
            eh_gt[sample] = []

        for  e,l in enumerate(i[1:]):
            pos, motif, gene, size = annot[e].rsplit(":",3)
            gt = l
            eh_gt[sample].append(EH(pos,motif,gene,size,gt))

trf = {}
for i in eh_gt:
    for items in eh_gt[i]:
        if not items.pos in trf:
            trf[items.pos] = {}
        trf[items.pos].update({i:items})

samples = list(eh_gt.keys())

xlsx = os.path.splitext(annot_tsv)[0] + "per-sample.xlsx"
workbook = Workbook(xlsx)

worksheet = workbook.add_worksheet("Allsamples")
header = ["#location", "repeat motif", "gene", "disease threshold"] + ["GT."+i for i in samples]
worksheet.write_row(0, 0, header)
row = 1
for i in trf:
        # for i in eh_gt[sample]:
        info = trf[i][samples[0]]
        content = [info.pos, info.motif, info.gene, info.size]
        content += [ trf[i][s].gt for s in samples ]            
        # print(i,content)
        worksheet.write_row(row, 0, content)
        row += 1
workbook.close()

    


        
