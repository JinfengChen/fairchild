#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def usage():
    test="name"
    message='''
python revision_FC_fix_gaps.py --input test.fasta 

    '''
    print message

def convert_gaps(fastafile):
    ofile = open('%s.fix_gaps.fa' %(fastafile), 'w')
    fastaid = defaultdict(lambda : str())
    gap_100 = 100*"N"
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[str(record.id)] = str(record.seq)
        sequence   = str(record.seq)
        sequence = re.sub(r'N{100000,}', gap_100, sequence)
        sequence = re.sub(r'N{1,100}', gap_100, sequence)
        newrecord = SeqRecord(Seq(sequence), id = record.id, description='')
        SeqIO.write(newrecord, ofile, 'fasta') 
    ofile.close()

def readtable(infile):
    data = defaultdict(lambda : int())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[0]] = 1
    return data


#HEG4v1.0_188	90484	23596..23653	adaptor:NGB00360.1
def readid(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                pos  = re.split(r'\.\.', unit[2])
                if not data.has_key(unit[0]):
                    data[unit[0]] = [unit[1], pos[0], pos[1]]
                    #print unit[0], unit[1], pos[0], pos[1]
    return data

#chr01   479773  481012  23      W       A123v1.0_12     1       1240    +
#read agp to store seq infromatioin
def readagp(infile):
    data = defaultdict(lambda : list())
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith(r'chr') and len(line) > 2: 
                unit = re.split(r'\t',line)
                data[unit[5]] = unit
    return data


#chr04   26913152        27255910        1705    W       A123v1.0_3349   1       342759  +
#chr04   27255911        27256204        1706    N       294     scaffold        yes     paired-ends;align_genus
#chr04   27256205        27332325        1707    W       A123v1.0_3350   1       76121   +

def maskregion(prefix,id_mask, agp, contig, id_remove):
    ofasta= '%s.masked.contig.fa' %(prefix)
    oagp  = '%s.temp.agp' %(prefix)
    agp_ref = readagp(agp)
    agp_split = defaultdict(lambda: str())
    ofile = open(ofasta, 'w')
    for seqid in sorted(contig.keys()):
        #print '>%s' %(seqid)
        if id_mask.has_key(seqid):
            print seqid, id_mask[seqid][0], id_mask[seqid][1], id_mask[seqid][2]
            length = int(id_mask[seqid][2]) - int(id_mask[seqid][1])
            mask   = 'N'*length
            if int(id_mask[seqid][2]) == int(id_mask[seqid][0]):
                newseq = contig[seqid][:int(id_mask[seqid][1])]
                newrecord = SeqRecord(Seq(newseq), id = seqid, description='')
                SeqIO.write(newrecord, ofile, 'fasta')
                print 'Length changed***'
            elif int(id_mask[seqid][1]) == 0:
                print 'Length changed***'
                newseq = contig[seqid][int(id_mask[seqid][1]):]
                newrecord = SeqRecord(Seq(newseq), id = seqid, description='')
                SeqIO.write(newrecord, ofile, 'fasta')
            else:
                newseq1 = contig[seqid][:int(id_mask[seqid][1])]
                newseq2 = contig[seqid][int(id_mask[seqid][2]):]
                seqid1  = '%s.1' %(seqid)
                seqid2  = '%s.2' %(seqid)
                newrecord1 = SeqRecord(Seq(newseq1), id = seqid1, description='')
                newrecord2 = SeqRecord(Seq(newseq2), id = seqid2, description='')
                SeqIO.write(newrecord1, ofile, 'fasta')
                SeqIO.write(newrecord2, ofile, 'fasta')
                #print '>%s\n%s' %(seqid, newseq)
                if agp_ref.has_key(seqid):
                    part1_end = int(agp_ref[seqid][1])+int(id_mask[seqid][1])-1
                    part1_len = int(id_mask[seqid][1])
                    part1 = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(agp_ref[seqid][0],agp_ref[seqid][1],part1_end,agp_ref[seqid][3],agp_ref[seqid][4],seqid1,agp_ref[seqid][6],part1_len,agp_ref[seqid][8])
                    gap_start = int(agp_ref[seqid][1]) + int(id_mask[seqid][1])
                    gap_end   = int(agp_ref[seqid][1]) + int(id_mask[seqid][2]) - 1
                    gap_len   = gap_end - gap_start + 1
                    gap   = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(agp_ref[seqid][0], gap_start, gap_end, agp_ref[seqid][3], 'N', gap_len, 'scaffold', 'yes', 'paired-ends')
                    part2_start = int(agp_ref[seqid][1]) + int(id_mask[seqid][2])
                    part2_len   = int(agp_ref[seqid][2]) - part2_start + 1
                    part2 = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' %(agp_ref[seqid][0], part2_start, agp_ref[seqid][2], agp_ref[seqid][3], agp_ref[seqid][4], seqid2, agp_ref[seqid][6],part2_len, agp_ref[seqid][8])
                    agp_split[seqid] = '%s\n%s\n%s' %(part1, gap, part2)
                    print '\t'.join(agp_ref[seqid])
                    print agp_split[seqid]
                #print 'orignial length: %s, newlength length: %s' %(len(contig[seqid]),len(newseq))
        elif id_remove.has_key(seqid):
            continue
        else:
            newrecord = SeqRecord(Seq(contig[seqid]), id = seqid, description='')
            SeqIO.write(newrecord, ofile, 'fasta')
    ofile.close()

    ofile = open(oagp, 'w')
    with open (agp, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith(r'chr') and len(line) > 2: 
                unit = re.split(r'\t',line)
                if unit[6] == 'scaffold':
                    print >> ofile, line
                elif agp_split.has_key(unit[5]):
                    print >> ofile, agp_split[unit[5]]
                else:
                    print >> ofile, line
            else:
                print >> ofile, line
    ofile.close()

## Generated from SOAPdenovo assembly file /rhome/cjinfeng/BigData/01.Rice_genomes/HEG4/00.Genome/Chromosome/HEG4_ALLPATHLG_v1.chr.fasta using script fasta2apg.p
#chr01   1       33508   1       W       HEG4v1.0_1      1       33508   +
#chr01   33509   33608   2       N       100     scaffold        yes     paired-ends;align_genus
#chr01   33609   147049  3       W       HEG4v1.0_2      1       113441  +
#chr01   147050  147244  4       N       195     scaffold        yes     paired-ends;align_genus
def update_agp(prefix, infile, id_remove):
    oagp = '%s.masked.agp' %(prefix)
    ofile = open(oagp, 'w')
    lastchr = ''
    lastrank = 0
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if line.startswith(r'#'):
                print >> ofile, line
            elif(len(line) > 2): 
                unit = re.split(r'\t',line)
                if unit[0] == lastchr and int(unit[3]) != lastrank + 1:
                    unit[3] = str(lastrank + 1)
                lastchr = unit[0]
                lastrank = int(unit[3])
                if id_remove.has_key(unit[5]):
                    continue
                elif unit[6] == 'scaffold':
                    if int(unit[5]) == 100:
                        unit[4] = 'U'
                        newline = '\t'.join(unit)
                        print >> ofile, newline
                    else:
                        unit[8] = 'paired-ends'
                        newline = '\t'.join(unit)
                        print >> ofile, newline
                else:
                    newline = '\t'.join(unit)
                    print >> ofile, newline
    ofile.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()

    convert_gaps(args.input)

if __name__ == '__main__':
    main()

