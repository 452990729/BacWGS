#!/usr/bin/env python2


import sys
import re
import os
from HTSeq import FastaReader

def RefFa(file_in):
    dict_fa = {}
    for item in FastaReader(file_in):
        list_lb = re.split('_', item.name)
        lb = '_'.join(list_lb[:2])
        dict_fa[lb] = str(item.seq)
    return dict_fa

def ReadCulste(file_in):
    dict_cls = {}
    with open(file_in, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                list_split = re.split('\t', line.strip('\n'))
                list_seq = [re.split('\(', list_split[0])[0],]+[re.split('\(', i)[0] for i in re.split('\s+', list_split[2])]
                lb = re.split('\(', list_split[0])[0]
                dict_cls[lb] = list_seq
    return dict_cls

def GetSeq(dict_fa, dict_cls, pan_id):
    out = open('AllSequence.fa', 'w')
    print dict_cls[pan_id]
    for seq in dict_cls[pan_id]:
        out.write('>{}\n{}\n'.format(re.split('_', seq)[0], dict_fa[seq]))
    out.close()

def main():
    dict_fa = RefFa(sys.argv[1])
    dict_cls = ReadCulste(sys.argv[2])
    GetSeq(dict_fa, dict_cls, sys.argv[3])


if __name__ == '__main__':
    main()

