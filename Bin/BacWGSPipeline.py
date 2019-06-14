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
CMSCAN = config.get('SOFTWARE', 'cmscan')
GMSN = config.get('SOFTWARE', 'gmsn')
GMHMMP = config.get('SOFTWARE', 'gmhmmp')
REPEATMASKER = config.get('SOFTWARE', 'RepeatMasker')
TRF = config.get('SOFTWARE', 'trf')
TRNASCAN = config.get('SOFTWARE', 'tRNAscan-SE')
RNAMMER = config.get('SOFTWARE', 'rnammer')
SEQKIT = config.get('SOFTWARE', 'seqkit')

### SCRIPT
HandlerRNA = config.get('SCRIPT', 'HandlerRNA')
genemark_convert_2 = config.get('SCRIPT', 'genemark_convert_2')
PGAP_stat = config.get('SCRIPT', 'PGAP_stat')
repeat_to_gff = config.get('SCRIPT', 'repeat_to_gff')
tRNAscan_to_gff3 = config.get('SCRIPT', 'tRNAscan_to_gff3')
sample_stat = config.get('SCRIPT', 'sample_stat')

#### DATABASE
RepeatmaskLib = config.get('DATABASE', 'RepeatmaskLib')
Rfam = config.get('DATABASE', 'Rfam')


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
    parser = argparse.ArgumentParser(description="Micro WGS pipeline")
    parser.add_argument('-c', help='the input fasta list', required=True)
    parser.add_argument('-o', help='the abs output path', required=True)
    parser.add_argument('-a1', help='the read1 adapter', default='AGATCGGAAGAGCACACGTCTGAACTCCAGTCA')
    parser.add_argument('-a2', help='the read2 adapter', default='AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT')
    parser.add_argument('-p', help='the type of spe',choices=['bac', 'euk', 'vir'], default='bac')
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
    specie= argv['p']
    outpath = argv['o']
    if specie == 'bac':
        spe = 'prok'
        spr = 'bac'
        spt = 'B'
    elif specie == 'euk':
        spe = 'euk'
        spr = 'euk'
        spt = 'E'
    elif specie == 'vir':
        spe = 'virus'
        spr = 'bac'
        spt = 'B'
    snakefile = open(os.path.join(argv['o'], 'snakefile.txt'), 'w')
    ### config file
    snakefile.write('Samples = "{}".split()\n'.format(' '.join([i.Name for i in\
                                                                list_ob])))
    snakefile.write('adapter1 = "{}"\nadapter2 = "{}"\n'.format(argv['a1'], argv['a2']))

    ###all
    All = Snake('All')
    All.UpdateInput('expand("'+outpath+'/6.Result/{sample}/{sample}.sRNA.gff", sample=Samples)')
    All.WriteStr(snakefile)

    ###QC
    QC = Snake('QC')
    if Paired:
        QC.UpdateInput('A = "'+outpath+'/0.RawData/{sample}_1'+lb+'", B = "'+outpath+'/0.RawData/{sample}_2'+lb+'"')
        QC.UpdateOutput('A = "'+outpath+'/1.QC/{sample}_1.clean.fastq.gz", B = "'+outpath+'/1.QC/{sample}_2.clean.fastq.gz"')
        QC.UpdateThreads('2')
        QC.UpdateLog('e = "logs/{sample}.qc.e", o = "logs/{sample}.qc.o"')
        QC.UpdateShell(r'"'+FASTP+r' -i {input.A} -o {output.A} -I {input.B} -O {output.B} --adapter_sequence {adapter1} --adapter_sequence_r2 {adapter2} -w {threads} -j '+outpath+'/1.QC/{wildcards.sample}_QC_report.json -h '+outpath+'/1.QC/{wildcards.sample}_QC_report.html 1>{log.o} 2>{log.e}"')
    else:
        QC.UpdateInput('"'+outpath+'/0.RawData/{sample}'+lb+'"')
        QC.UpdateOutput('"'+outpath+'/1.QC/{sample}.clean.fastq.gz"')
        QC.UpdateThreads('2')
        QC.UpdateLog('e = "logs/{sample}.qc.e", o = "logs/{sample}.qc.o"')
        QC.UpdateShell(r'"'+FASTP+r' -i {input} -o {output}  --adapter_sequence {adapter1} -w {threads} -j '+outpath+'/1.QC/{wildcards.sample}_QC_report.json -h '+outpath+'/1.QC/{wildcards.sample}_QC_report.html 1>{log.o} 2>{log.e}"')
    QC.WriteStr(snakefile)

    ### Assemble
    Assemble = Snake('Assemble')
    if Paired:
        Assemble.UpdateInput('A = "'+outpath+'/1.QC/{sample}_1.clean.fastq.gz", B = "'+outpath+'/1.QC/{sample}_2.clean.fastq.gz"')
        Assemble.UpdateOutput('"'+outpath+'/2.Assemble/{sample}/scaffolds.fasta"')
        Assemble.UpdateThreads('5')
        Assemble.UpdateLog('e = "logs/{sample}.align.e", o = "logs/{sample}.align.o"')
        Assemble.UpdateShell(r'"'+SPADE+r' -1 {input.A} -2 {input.B} -t {threads} -o '+outpath+'/2.Assemble/{wildcards.sample}/ 1>{log.o} 2>{log.e}"')
    else:
        Assemble.UpdateInput('"'+outpath+'/1.QC/{sample}.clean.fq.gz"')
        Assemble.UpdateOutput('"'+outpath+'/2.Assemble/{sample}/scaffolds.fasta"')
        Assemble.UpdateThreads('5')
        Assemble.UpdateLog('e = "logs/{sample}.assemble.e", o = "logs/{sample}.assemble.o"')
        Assemble.UpdateShell(r'"'+SPADE+r' -s {input} -t {threads} -o '+outpath+'/2.Assemble/{wildcards.sample}/ 1>{log.o} 2>{log.e}"')
    Assemble.WriteStr(snakefile)

    ### GenePredict1
    GenePredict1 = Snake('GenePredict1')
    GenePredict1.UpdateInput('"'+outpath+'/2.Assemble/{sample}/scaffolds.fasta"')
    GenePredict1.UpdateOutput('"'+outpath+'/3.GenePredict/{sample}/{sample}_hmm.mod"')
    GenePredict1.UpdateThreads('1')
    GenePredict1.UpdateLog('e = "logs/{sample}.gp1.e", o = "logs/{sample}.gp1.o"')
    GenePredict1.UpdateShell('"cd '+outpath+'/3.GenePredict/{wildcards.sample}/;"\n\t\t"'+GMSN+' -name {wildcards.sample} -clean -gcode 11 -shape partial --combine --'+spe+r' {input} 1>{log.o} 2>{log.e}"')
    GenePredict1.WriteStr(snakefile)

    ### GenePredict2
    GenePredict2 = Snake('GenePredict2')
    GenePredict2.UpdateInput('A = "'+outpath+'/3.GenePredict/{sample}/{sample}_hmm.mod", B = "'+outpath+'/2.Assemble/{sample}/scaffolds.fasta"')
    GenePredict2.UpdateOutput('"'+outpath+'/3.GenePredict/{sample}/{sample}.gmhmmp"')
    GenePredict2.UpdateThreads('1')
    GenePredict2.UpdateLog('e = "logs/{sample}.gp2.e", o = "logs/{sample}.gp2.o"')
    GenePredict2.UpdateShell('"cd '+outpath+'/3.GenePredict/{wildcards.sample}/;"\n\t\t"'+GMHMMP+' -m {input.A} -o {output} -a -d -p 1 -f L {input.B} 1>{log.o} 2>{log.e}"')
    GenePredict2.WriteStr(snakefile)

    ### GenePredict3
    GenePredict3 = Snake('GenePredict3')
    GenePredict3.UpdateInput('A = "'+outpath+'/3.GenePredict/{sample}/{sample}.gmhmmp", B = "'+outpath+'/2.Assemble/{sample}/scaffolds.fasta"')
    GenePredict3.UpdateOutput('A = "'+outpath+'/3.GenePredict/{sample}/{sample}.gmhmmp.cds", B = "'+outpath+'/3.GenePredict/{sample}/{sample}.gmhmmp.pep", C = "'+outpath+'/3.GenePredict/{sample}/{sample}.gmhmmp.gff"')
    GenePredict3.UpdateThreads('1')
    GenePredict3.UpdateLog('e = "logs/{sample}.gp3.e", o = "logs/{sample}.gp3.o"')
    GenePredict3.UpdateShell('"cd '+outpath+'/3.GenePredict/{wildcards.sample}/;"\n\t\t"'+genemark_convert_2+' --final {wildcards.sample}GM --gcode 11 --log --verbose {input.A} {input.B} 1>{log.o} 2>{log.e}"')
    GenePredict3.WriteStr(snakefile)

    ### Repeat1
    Repeat1 = Snake('Repeat1')
    Repeat1.UpdateInput('"'+outpath+'/2.Assemble/{sample}/scaffolds.fasta"')
    Repeat1.UpdateOutput('A = "'+outpath+'/4.Repeat/{sample}/scaffolds.fasta.out", B = "'+outpath+'/4.Repeat/{sample}/scaffolds.fasta.out.gff"')
    Repeat1.UpdateThreads('2')
    Repeat1.UpdateLog('e = "logs/{sample}.rp1.e", o = "logs/{sample}.rp1.o"')
    Repeat1.UpdateShell('"cd '+outpath+'/4.Repeat/{wildcards.sample}/;"\n\t\t"'+REPEATMASKER+' -nolow -no_is -norna -engine wublast -parallel {threads} -lib '+RepeatmaskLib+' {input} -dir '+outpath+'/4.Repeat/{wildcards.sample}/;"\n\t\t"'+repeat_to_gff+' --prefix {wildcards.sample} {output.A} 1>{log.o} 2>{log.e}"')
    Repeat1.WriteStr(snakefile)

    ### Repeat2
    Repeat2 = Snake('Repeat2')
    Repeat2.UpdateInput('"'+outpath+'/2.Assemble/{sample}/scaffolds.fasta"')
    Repeat2.UpdateOutput('A = "'+outpath+'/4.Repeat/{sample}/scaffolds.fasta.2.7.7.80.10.50.2000.dat", B = "'+outpath+'/4.Repeat/{sample}/scaffolds.fasta.2.7.7.80.10.50.2000.dat.gff"')
    Repeat2.UpdateLog('e = "logs/{sample}.rp2.e", o = "logs/{sample}.rp2.o"')
    Repeat2.UpdateShell('"cd '+outpath+'/4.Repeat/{wildcards.sample};"\n\t\t"'+TRF+' {input} 2 7 7 80 10 50 2000 -d -h;"\n\t\t"'+repeat_to_gff+' --prefix {wildcards.sample} {output.A} 1>{log.o} 2>{log.e}"')
    Repeat2.WriteStr(snakefile)

    ### rRNA
    rRNA = Snake('rRNA')
    rRNA.UpdateInput('"'+outpath+'/2.Assemble/{sample}/scaffolds.fasta"')
    rRNA.UpdateOutput('A = "'+outpath+'/5.RNA/rRNA/{sample}/{sample}.rRNAd.gff", B = "'+outpath+'/5.RNA/rRNA/{sample}/{sample}.rRNAd.fa"')
    rRNA.UpdateLog('e = "logs/{sample}.rRNA.e", o = "logs/{sample}.rRNA.o"')
    rRNA.UpdateShell('"cd '+outpath+'/5.RNA/rRNA/{wildcards.sample}/;"\n\t\t"'+PERL+' '+RNAMMER+' -S '+spr+' -m tsu,lsu,ssu -gff {output.A} -f {output.B} {input};"\n\t\t"'+HandlerRNA+' '+outpath+'/5.RNA/rRNA/{wildcards.sample} 1>{log.o} 2>{log.e}"')
    rRNA.WriteStr(snakefile)

    ### tRNA
    tRNA = Snake('tRNA')
    tRNA.UpdateInput('"'+outpath+'/2.Assemble/{sample}/scaffolds.fasta"')
    tRNA.UpdateOutput('A = "'+outpath+'/5.RNA/tRNA/{sample}/{sample}.tRNA", B = "'+outpath+'/5.RNA/tRNA/{sample}/{sample}.tRNA.structure", C = "'+outpath+'/5.RNA/tRNA/{sample}/{sample}.tRNA.gff"')
    tRNA.UpdateLog('e = "logs/{sample}.tRNA.e", o = "logs/{sample}.tRNA.o"')
    tRNA.UpdateShell('"cd '+outpath+'/5.RNA/tRNA/{wildcards.sample}/;"\n\t\t"'+TRNASCAN+' -'+spt+' -o {output.A} -f {output.B} {input};"\n\t\t"'+tRNAscan_to_gff3+' --prefix {wildcards.sample} {output.A} {output.B} > {output.C} 1>{log.o} 2>{log.e}"')
    tRNA.WriteStr(snakefile)

    ### sRNA
    sRNA = Snake('sRNA')
    sRNA.UpdateInput('"'+outpath+'/2.Assemble/{sample}/scaffolds.fasta"')
    sRNA.UpdateOutput('A = "'+outpath+'/5.RNA/sRNA/{sample}/{sample}.tblout", B = "'+outpath+'/5.RNA/sRNA/{sample}/{sample}.cmscan"')
    sRNA.UpdateLog('e = "logs/{sample}.sRNA.e", o = "logs/{sample}.sRNA.o"')
    sRNA.UpdateShell('"cd '+outpath+'/5.RNA/sRNA/{wildcards.sample}/;"\n\t\t"'+CMSCAN+' -Z `'+SEQKIT+' stats -T {input} | awk \'{{if(NR==2) print int($5/2000000)+1}}\'` --cut_ga --rfam --nohmmonly --tblout {output.A} --fmt 2 --clanin '+Rfam+'/Rfam.clanin '+Rfam+'/Rfam.cm {input} > {output.B} 1>{log.o} 2>{log.e}"')
    sRNA.WriteStr(snakefile)

    ### TotalGff
    TotalGff = Snake('TotalGff')
    TotalGff.UpdateInput('cds = "'+outpath+'/3.GenePredict/{sample}/{sample}.gmhmmp.cds",\
                         pep = "'+outpath+'/3.GenePredict/{sample}/{sample}.gmhmmp.pep",\
                         gene = "'+outpath+'/3.GenePredict/{sample}/{sample}.gmhmmp.gff",\
                         repbase = "'+outpath+'/4.Repeat/{sample}/scaffolds.fasta.out.gff",\
                         trf = "'+outpath+'/4.Repeat/{sample}/scaffolds.fasta.2.7.7.80.10.50.2000.dat.gff",\
                         rRNA_denovo = "'+outpath+'/5.RNA/rRNA/{sample}/{sample}.rRNAd.gff",\
                         tRNA = "'+outpath+'/5.RNA/tRNA/{sample}/{sample}.tRNA.gff",\
                         sRNA = "'+outpath+'/5.RNA/sRNA/{sample}/{sample}.tblout"')
    TotalGff.UpdateOutput('cds = "'+outpath+'/6.Result/{sample}/{sample}.cds",\
                          pep = "'+outpath+'/6.Result/{sample}/{sample}.pep",\
                          gene = "'+outpath+'/6.Result/{sample}/{sample}.gff",\
                          repbase = "'+outpath+'/6.Result/{sample}/{sample}.repbase.gff",\
                          trf = "'+outpath+'/6.Result/{sample}/{sample}.trf.gff",\
                          rRNA_denovo = "'+outpath+'/6.Result/{sample}/{sample}.rRNA.gff",\
                          tRNA = "'+outpath+'/6.Result/{sample}/{sample}.tRNA.gff",\
                          sRNA = "'+outpath+'/6.Result/{sample}/{sample}.sRNA.gff"')
    TotalGff.UpdateLog('e = "logs/{sample}.result.e", o = "logs/{sample}.result.o"')
    TotalGff.UpdateShell('"ln -s {input.cds} {output.cds};"\n\t\t"ln -s {input.pep} {output.pep};"\n\t\t"ln -s {input.gene} {output.gene};"\n\t\t"ln -s {input.repbase} {output.repbase};"\n\t\t"ln -s {input.trf} {output.trf};"\n\t\t"ln -s {input.rRNA_denovo} {output.rRNA_denovo};"\n\t\t"ln -s {input.tRNA} {output.tRNA};"\n\t\t"ln -s {input.sRNA} {output.sRNA}"')
    TotalGff.WriteStr(snakefile)



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
