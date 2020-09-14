#!/usr/bin/env python2

import sys
import re
import os
import argparse
import ConfigParser
from glob import glob

BasePath = os.path.split(os.path.realpath(__file__))[0]
config = ConfigParser.ConfigParser()
config.read(BasePath+'/../Bin/config.ini')

#### SOFT
PYTHON = config.get('SOFTWARE', 'python')
PERL = config.get('SOFTWARE', 'perl')
SNAKEMAKE = config.get('SOFTWARE', 'snakemake')
FASTP = config.get('SOFTWARE', 'fastp')
SPADE = config.get('SOFTWARE', 'spade')
BWA = config.get('SOFTWARE', 'bwa')
SAMTOOLS = config.get('SOFTWARE', 'samtools')

### SCRIPT
HandlerRNA = config.get('SCRIPT', 'HandlerRNA')

#### DATABASE


class ReadList(object):
    def __init__(self, line_in):
        list_split = re.split('\s+', line_in)
        self.Sample = list_split[0]
        self.Group = list_split[1]
        self.Name = '_'.join(list_split[:2])
        if ',' in list_split[2]:
            self.paired = True
            list_fq = re.split(',', list_split[2])
            self.fq1 = list_fq[0]
            self.fq2 = list_fq[1]
        else:
            self.paired = False
            self.fq = list_split[2]

class Snake(object):
    def __init__(self, process):
        self.process = process
        self.input = ''
        self.output = ''
        self.params = ''
        self.log = ''
        self.threads = ''
        self.shell = ''

    def UpdateInput(self, line_in):
        self.input = line_in

    def UpdateOutput(self, line_in):
        self.output = line_in

    def UpdateParams(self, line_in):
        self.params = line_in

    def UpdateLog(self, line_in):
        self.log = line_in

    def UpdateThreads(self, line_in):
        self.threads = line_in

    def UpdateShell(self, line_in):
        self.shell = line_in

    def WriteStr(self, fn):
        fn.write('rule '+self.process+':\n')
        fn.write('\tinput:\n\t\t'+self.input+'\n')
        if self.output:
            fn.write('\toutput:\n\t\t'+self.output+'\n')
        if self.params:
            fn.write('\tparams:\n\t\t'+self.params+'\n')
        if self.log:
            fn.write('\tlog:\n\t\t'+self.log+'\n')
        if self.threads:
            fn.write('\tthreads: '+self.threads+'\n')
        if self.shell:
            fn.write('\tshell:\n\t\t'+self.shell+'\n')
        fn.write('\n')


def Argparse():
    parser = argparse.ArgumentParser(description="Micro WGS Assemble pipeline")
    parser.add_argument('-c', help='the input fasta list', required=True)
    parser.add_argument('-o', help='the abs output path', required=True)
    parser.add_argument('-a1', help='the read1 adapter', default='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA')
    parser.add_argument('-a2', help='the read2 adapter', default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
    parser.add_argument('-t', help='the core used', default=10)
    parser.add_argument('-sge', help='use sge',  action='store_true')
    parser.add_argument('-lsf', help='use lsf',  action='store_true')
    parser.add_argument('-r', help='run now',  action='store_true')
    argv=vars(parser.parse_args())
    return argv

def HandleRawdata(argv):
    outpath = argv['o']
    RawData = os.path.join(outpath, '0.RawData')
    list_ob = []
    if not os.path.exists(RawData):
        os.mkdir(RawData)
    else:
        os.system('rm -rf '+RawData)
        os.mkdir(RawData)
    with open(argv['c'], 'r') as f:
        os.chdir(RawData)
        for line in f:
            if not line.startswith('#'):
                ob = ReadList(line.strip())
                list_ob.append(ob)
                if ob.paired:
                    Paired = True
                    if ob.fq1.endswith('.gz'):
                        lb = '.fastq.gz'
                    else:
                        lb = '.fastq'
                    os.system('ln -s {} {}'.format(ob.fq1, ob.Name+'_1'+lb))
                    os.system('ln -s {} {}'.format(ob.fq2, ob.Name+'_2'+lb))
                else:
                    Paired = False
                    if ob.fq.endswith('.gz'):
                        lb = '.fastq.gz'
                    else:
                        lb = '.fastq'
                    os.system('ln -s {} {}'.format(ob.fq, ob.Name+lb))
    os.chdir(outpath)
    return list_ob, Paired, lb

def WriteSnake(argv, list_ob, Paired, lb):
    outpath = argv['o']
    snakefile = open(os.path.join(argv['o'], 'snakefile.txt'), 'w')
    ### config file
    snakefile.write('Samples = "{}".split()\n'.format(' '.join([i.Name for i in\
                                                                list_ob])))
    snakefile.write('adapter1 = "{}"\nadapter2 = "{}"\n'.format(argv['a1'], argv['a2']))

    ###all
    All = Snake('All')
    All.UpdateInput('expand("'+outpath+'/4.Depth/{sample}/{sample}.depth", sample=Samples)')
    All.WriteStr(snakefile)

    ###QC
    QC = Snake('QC')
    if Paired:
        QC.UpdateInput('A = "'+outpath+'/0.RawData/{sample}_1'+lb+'", B = "'+outpath+'/0.RawData/{sample}_2'+lb+'"')
        QC.UpdateOutput('A = "'+outpath+'/1.QC/{sample}_1.clean.fastq.gz", B = "'+outpath+'/1.QC/{sample}_2.clean.fastq.gz"')
        QC.UpdateThreads('2')
        QC.UpdateLog('e = "'+outpath+'logs/{sample}.qc.e", o = "'+outpath+'logs/{sample}.qc.o"')
        QC.UpdateShell(r'"'+FASTP+r' -i {input.A} -o {output.A} -I {input.B} -O {output.B} --adapter_sequence {adapter1} --adapter_sequence_r2 {adapter2} -w {threads} -j '+outpath+'/1.QC/{wildcards.sample}_QC_report.json -h '+outpath+'/1.QC/{wildcards.sample}_QC_report.html 1>{log.o} 2>{log.e}"')
    else:
        QC.UpdateInput('"'+outpath+'/0.RawData/{sample}'+lb+'"')
        QC.UpdateOutput('"'+outpath+'/1.QC/{sample}.clean.fastq.gz"')
        QC.UpdateThreads('2')
        QC.UpdateLog('e = "'+outpath+'logs/{sample}.qc.e", o = "'+outpath+'logs/{sample}.qc.o"')
        QC.UpdateShell(r'"'+FASTP+r' -i {input} -o {output}  --adapter_sequence {adapter1} -w {threads} -j '+outpath+'/1.QC/{wildcards.sample}_QC_report.json -h '+outpath+'/1.QC/{wildcards.sample}_QC_report.html 1>{log.o} 2>{log.e}"')
    QC.WriteStr(snakefile)

    ### Assemble
    Assemble = Snake('Assemble')
    if Paired:
        Assemble.UpdateInput('A = "'+outpath+'/1.QC/{sample}_1.clean.fastq.gz", B = "'+outpath+'/1.QC/{sample}_2.clean.fastq.gz"')
        Assemble.UpdateOutput('"'+outpath+'/2.Assemble/{sample}/scaffolds.fasta"')
        Assemble.UpdateThreads('5')
        Assemble.UpdateLog('e = "'+outpath+'logs/{sample}.align.e", o = "'+outpath+'logs/{sample}.align.o"')
        Assemble.UpdateShell(r'"'+SPADE+r' -1 {input.A} -2 {input.B} -t {threads} -o '+outpath+'/2.Assemble/{wildcards.sample}/ 1>{log.o} 2>{log.e}"')
    else:
        Assemble.UpdateInput('"'+outpath+'/1.QC/{sample}.clean.fq.gz"')
        Assemble.UpdateOutput('"'+outpath+'/2.Assemble/{sample}/scaffolds.fasta"')
        Assemble.UpdateThreads('5')
        Assemble.UpdateLog('e = "'+outpath+'logs/{sample}.assemble.e", o = "'+outpath+'logs/{sample}.assemble.o"')
        Assemble.UpdateShell(r'"'+SPADE+r' -s {input} -t {threads} -o '+outpath+'/2.Assemble/{wildcards.sample}/ 1>{log.o} 2>{log.e}"')
    Assemble.WriteStr(snakefile)

    ### Align
    Align = Snake('Align')
    Align.UpdateInput('fasta = "'+outpath+'/2.Assemble/{sample}/scaffolds.fasta", fq1 = "'+outpath+'/1.QC/{sample}_1.clean.fastq.gz", fq2 = "'+outpath+'/1.QC/{sample}_2.clean.fastq.gz"')
    Align.UpdateOutput('"'+outpath+'/3.Align/{sample}/{sample}.bam"')
    Align.UpdateThreads('1')
    Align.UpdateLog('e = "'+outpath+'logs/{sample}.align.e", o = "'+outpath+'logs/{sample}.align.o"')
    Align.UpdateShell('"cd '+outpath+'/3.Align/{wildcards.sample}/;"\n\t\t"'+BWA+' index -p reference {input.fasta};"\n\t\t"'+BWA+' mem -t {threads} -M reference {input.fq1} {input.fq2} |'+SAMTOOLS+r' view -Sb -@ 4 -T {input.fasta} -o {output} 1>{log.o} 2>{log.e}"')
    Align.WriteStr(snakefile)

    ### Depth
    Depth = Snake('Depth')
    Depth.UpdateInput('"'+outpath+'/3.Align/{sample}/{sample}.bam"')
    Depth.UpdateOutput('"'+outpath+'/4.Depth/{sample}/{sample}.depth"')
    Depth.UpdateThreads('1')
    Depth.UpdateLog('e = "'+outpath+'logs/{sample}.depth.e", o = "'+outpath+'logs/{sample}.depth.o"')
    Depth.UpdateShell('"cd '+outpath+'/4.Depth/{wildcards.sample}/;"\n\t\t"'+SAMTOOLS+' depth -a {input} 1>{output} 2>{log.e}"')
    Depth.WriteStr(snakefile)


def RunShell(argv):
    out = open(os.path.join(argv['o'], 'runsnake.sh'), 'w')
    if argv['lsf']:
        out.write(SNAKEMAKE+' '+' '.join(['--cores', argv['t'], '--cluster', "'bsub -q normal -n {threads} -o %J.o -e %J.e'",\
                                      '--printshellcmds', '--snakefile', 'snakefile.txt']))
    else:
        out.write(SNAKEMAKE+' '+' '.join(['--cores', argv['t'], '--printshellcmds', '--snakefile', 'snakefile.txt']))
    out.close()
    if argv['r']:
        os.system('nohup sh runsnake.sh& > snake.log')

def main():
    argv = Argparse()
    list_ob, Paired, lb = HandleRawdata(argv)
    WriteSnake(argv, list_ob, Paired, lb)
    RunShell(argv)


if __name__ == '__main__':
    main()
