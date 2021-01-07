import re
import pandas as pd
import glob

## Format of the config file
# | chromName | faFile | MapFile                 | readPosFile | binFileNorm    |
# |     1     |  1.fa  | hg19.CRC.50mer.chr1.txt |     1.seq   | 1.norm.bin     |

import argparse

parser = argparse.ArgumentParser()


parser.add_argument("lib_id", help = "Library ID, string")
parser.add_argument("seqpath", help = "Path to .seq files")
parser.add_argument("outfile", help = "path to outfile")

parser = parser.parse_args()

print(parser)

chroms = range(1, 23)
chroms = [str(i) for i in chroms] + ["X"]

chromName = chroms

faFile = ["/projects/lculibrk_prj/CNV/bicseq2/resources/" + chrom + ".fa" for chrom in chroms]

MapFile = ["/projects/lculibrk_prj/CNV/bicseq2/resources/hg19.50mer.CRC.chr" + chrom + ".txt" for chrom in chroms]

readPosFile = [parser.seqpath + chrom + ".seq" for chrom in chroms]

binFileNorm = [parser.seqpath + chrom + ".norm.bin" for chrom in chroms]

tab = pd.DataFrame({'chromName':chromName, 'faFile':faFile, 'MapFile':MapFile, 'readPosFile':readPosFile, 'binFileNorm':binFileNorm})

with open(parser.outfile, "w") as f:
    tab.to_csv(f, sep = "\t", index = False)