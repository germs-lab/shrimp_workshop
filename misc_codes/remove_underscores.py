import sys

for line in open(sys.argv[1]):
    dat = line.rstrip().split('\t')
    filename = dat[0].rsplit('_',2)[0]
    print filename.split('_seqs')[1][1:] + '\t' + dat[1]
