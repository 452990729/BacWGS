#!/usr/bin/env python2


import sys
import re
import os
import numpy as np
import pandas as pd
from HTSeq import FastaReader


def HandleData(dir_in):
    Input = os.path.join(dir_in, 'input.fa')
    Output = os.path.join(dir_in, 'output.fa')
    with open(Input, 'r') as f:
        for line in f:
            Consus = line.strip().lstrip('>')
            break
    dict_fa = {}
    for item in FastaReader(Output):
        dict_fa[item.name] = str(item.seq)
    total = len(dict_fa[Consus])
    dict_burdon = {}
    for key in dict_fa:
        m = 0
        for i in range(total):
            if dict_fa[key][i] != dict_fa[Consus][i]:
                m += 1
        lb = re.findall('(\w+_\w+)GM', key)[0]
        dict_burdon[lb] = str(round(float(m*1000/total), 4))
    return dict_burdon, Consus

def SelectAll(path_in, outfile):
    root, dirs, files = next(os.walk(path_in))
    dict_data = {}
    set_all = set()
    for dr in dirs:
        dict_burdon, Consus = HandleData(os.path.join(root, dr))
        dict_data[Consus] = dict_burdon
        set_all = set_all|set(dict_burdon.keys())
    list_all = sorted(list(set_all))
    pd_data = pd.DataFrame(columns=list_all)
    label = []
    for key in dict_data:
        list_tmp = []
        for lb in list_all:
            if lb in dict_data[key]:
                list_tmp.append(dict_data[key][lb])
            else:
                list_tmp.append('1000.0')
        label.append(key)
        sub = pd.DataFrame(np.array(list_tmp).reshape(1, len(list_tmp)), columns=list_all)
        pd_data = pd_data.append(sub, ignore_index=True)
    pd_data.T.to_csv(outfile, sep='\t', index=True, header=label)

def main():
    SelectAll(sys.argv[1], sys.argv[2])


if __name__ == '__main__':
    main()




