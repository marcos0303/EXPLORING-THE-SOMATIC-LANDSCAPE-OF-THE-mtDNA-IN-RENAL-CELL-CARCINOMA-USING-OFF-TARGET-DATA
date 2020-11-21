import sys
import os
import re
import argparse
import numpy as np

parser=argparse.ArgumentParser(description='Arguments to hsmetrics_parse.py')
parser.add_argument('--file', required=True)
parser.add_argument('--outdir', required=True)
parser.add_argument('--outhsmetrics', required=True) # something in .csv format for example (several results will be appended)
parser.add_argument('--outhistogram', required=True) # something in .csv format for example (single file per sample)


args=parser.parse_args()
filename=args.file
outputdir = args.outdir
outhistogram = args.outhistogram
outhsmetrics = args.outhsmetrics
infile=open(filename).read()

# Split the input file into the different sections which are separated by '##' headers
text = infile.split('##')


# Work the metrics part --> metrics = text[3]
metrics = text[3].split('\n')
units = metrics[1]
values = metrics[2]

    # Output

if not os.path.exists(str(outputdir) + '/' + str(outhsmetrics)):
    with open(str(outputdir) + '/' + str(outhsmetrics), 'w') as f:
        pass
        units = 'Sample' + '\t' + metrics[1] + '\n'
        values = filename + '\t' + values + '\n'
        f.write(units)
        f.write(values)
        f.close()

else:
    with open(str(outputdir) + '/' + str(outhsmetrics), 'a') as f:
        values = filename + '\t' + values + '\n'
        f.write(values)
        f.close()

# Work the histogram part --> histo = text[4]
histo = text[4].split('\n')
histo = list(filter(None, histo))

histo_var = histo[1] # histo_var 'coverage_or_base_quality\thigh_quality_coverage_count\tunfiltered_baseq_count\thigh_quality_coverage_count\tunfiltered_baseq_count\thigh_quality_coverage_count\tunfiltered_baseq_count'
histo_val= histo[2:] # histo_val[0] =  '0\t400290\t0\t400290\t0\t400290\t0'  ; histo_val[500] = '500\t313800\t0\t313800\t0\t313800\t0'

x = list()
y = list()
for i in histo_val: # extract the histogram
    aux = i.split('\t')
    x.append(aux[0])
    y.append(aux[1])

y_rev = list(reversed(y))
cdf_rev = list(np.cumsum(list(map(int,y_rev))))
cdf = np.asarray(list(reversed(cdf_rev)))
pct = cdf/sum(list(map(int,y)))


# Output parsed and calculated information

with open(str(outputdir) + '/' + str(outhistogram) + '_histogram.csv', 'w') as f:
    x = '\t'.join(x)
    f.write('coverage' + '\t' + x + '\n')
    f.write('count' + '\t' + '\t'.join(y) + '\n')
    f.write('cdf' + '\t' + '\t'.join(list(map(str,cdf))) + '\n')
    f.write('pct' + '\t' + '\t'.join(list(map(str,pct))) + '\n')
f.close()
