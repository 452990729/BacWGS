#!/usr/bin/env python2


import sys
import os
from glob import glob

BasePath = os.path.split(os.path.realpath(__file__))[0]
delete_dupLine = BasePath+'/delete_dupLine.pl'
delete_duprRNA = BasePath+'/delete_duprRNA.pl'

gff = glob(sys.argv[1]+'/*.gff')[0]
fa = glob(sys.argv[1]+'/*.fa')[0]

os.system("{} {} {} '^#'".format(delete_dupLine, gff, sys.argv[1]+'/gff2'))
os.system('{} {} {}'.format(delete_duprRNA, fa, sys.argv[1]+'/faa'))
os.system('mv {} {}'.format(sys.argv[1]+'/gff2', gff))
os.system('mv {} {}'.format(sys.argv[1]+'/faa', fa))
