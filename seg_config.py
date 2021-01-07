import re
import pandas as pd
import glob

## Format of the config file
# | chromName | binFileNorm.Case | binFileNorm.Control
# |     1     |  chr1.norm.bin   | chr1.norm.bin

import argparse

parser = argparse.ArgumentParser()


parser.add_argument("normpath", help = "Path to normal files")
parser.add_argument("tumpath", help = "Path to tumor files")
parser.add_argument("outfile", help = "path to outfile")

parser = parser.parse_args()

print(parser)

chroms = range(1, 23)
chroms = [str(i) for i in chroms] + ["X"]

chromName = chroms

binFileNorm_control = [parser.normpath + chrom + ".norm.bin" for chrom in chroms]
binFileNorm_case = [parser.tumpath + chrom + ".norm.bin" for chrom in chroms]

tab = pd.DataFrame({'chromName':chromName, 'binFileNorm.Case':binFileNorm_case, 'binFileNorm.Control':binFileNorm_control})

with open(parser.outfile, "w") as f:
    tab.to_csv(f, sep = "\t", index = False)