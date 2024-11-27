
'''
Copyright (c) 2020 Children's Hospital of Philadelphia
Author: Li Fang (https://github.com/fangli08)

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

import os
import sys

import gzip
from datetime import datetime

TimeFormat = '%m/%d/%Y %H:%M:%S'

### IO ### 

def read_list_file(in_list_file, abspath = False):
    out_list = list()
    in_list_fp = open(in_list_file, 'r')
    lines = list(in_list_fp)
    in_list_fp.close()

    for line in lines:
        line = line.strip()
        if abspath == True:
            line = os.path.abspath(line)
            out_list.append(line)
        else:
            out_list.append(line)
    return out_list

def gzopen(in_file, mode = 'rt'):

    if '.gz' == in_file[-3:] or '.bgz' == in_file[-4:]:
        in_fp = gzip.open(in_file, 'rt')
    else:
        in_fp = open(in_file, 'rt')

    return in_fp

def create_dir(dir):

    if os.path.exists(dir) == False:
        try:
            os.makedirs(dir)
        except:
            eprint('ERROR: Failed to create directory: %s' % dir)
            sys.exit()
    return

def check_input_file_exists(in_file):

    if os.path.exists(in_file) == False:
        eprint('ERROR! file does not exist: %s' % in_file)
        sys.exit()

    return

def eprint(message):
    sys.stderr.write('[' + datetime.now().strftime(TimeFormat) + '] ' + message + '\n')

    return

def get_file_prefix(input_file):

    return os.path.splitext(os.path.split(input_file)[1])[0]

### FASTQ/FASTA ###

def read_fasta_file(fasta_file):

    fasta_fp = gzopen(fasta_file)

    fasta_name_list = list()
    fasta_seq_list = list()
    curr_name = ''
    curr_seq = ''

    while 1:
        line = fasta_fp.readline()
        if not line: break
        line = line.strip()
        if not line: continue 
        if line[0] ==  '>':
            if len(curr_seq) > 0 and len(curr_name) > 0:
                fasta_name_list.append(curr_name)
                fasta_seq_list.append(curr_seq)
            curr_name = line[1:].split()[0]
            curr_seq = ''
            continue

        curr_seq += line
     
    fasta_fp.close()

    if len(curr_seq) > 0 and len(curr_name) > 0:
        fasta_name_list.append(curr_name)
        fasta_seq_list.append(curr_seq)
    
    return fasta_name_list, fasta_seq_list

def split_fastq(in_fastq_file_list, num_out_file, out_prefix):

    out_readname_set = set() # used to skip duplicated reads
    out_fastq_file_list = list()
    out_fastq_fp_list = list()
    for i in range(0, num_out_file):
        out_fastq_file = out_prefix + '.part%d.fastq' % i
        out_fastq_file_list.append(out_fastq_file)
        out_fastq_fp = open(out_fastq_file, 'w')
        out_fastq_fp_list.append(out_fastq_fp)

    for in_fastq_file in in_fastq_file_list:
        in_fastq_fp = gzopen(in_fastq_file)

        i = 0
        status = 0
        while 1:
            if status == 0:
                line1 = in_fastq_fp.readline()
            line2 = in_fastq_fp.readline()
            line3 = in_fastq_fp.readline()
            line4 = in_fastq_fp.readline()
            if not line1: break
            if not line2: break
            if not line3: break
            if not line4: break

            if line1[0] != '@' or len(line2) != len(line4) or line3.strip() != '+':
                eprint('ERROR! Bad fastq file: %s' % in_fastq_file)
                break
                

            i += 1
            file_id = i % num_out_file
            readname = line1.strip().split()[0][1:]
            if readname not in out_readname_set:
                out_fastq_fp_list[file_id].write(line1)
                out_fastq_fp_list[file_id].write(line2)
                out_fastq_fp_list[file_id].write(line3)
                out_fastq_fp_list[file_id].write(line4)

            out_readname_set.add(readname)

        in_fastq_fp.close()

    for out_fastq_fp in out_fastq_fp_list:
        out_fastq_fp.close()

    return out_fastq_file_list

def extract_fastq_tail_seq(in_fastq_file, read_tail_length, left_tail_fastq_file, right_tail_fastq_file):

    in_fastq_fp = gzopen(in_fastq_file)
    fq_left_tail_fp  = open(left_tail_fastq_file, 'w')
    fq_right_tail_fp = open(right_tail_fastq_file, 'w')

    while 1:
        line1 = in_fastq_fp.readline()
        line2 = in_fastq_fp.readline()
        line3 = in_fastq_fp.readline()
        line4 = in_fastq_fp.readline()

        if not line1: break
        if not line2: break
        if not line3: break
        if not line4: break

        if line1[0] != '@' or len(line2) != len(line4) or line3.strip() != '+':
            eprint('ERROR! bad fastq record: ')
            eprint(line1 + line2 + line3 + line4)
            sys.exit()
        
        read_seq  = line2.strip()
        read_qual = line4.strip()

        left_tail_seq  = read_seq[0:read_tail_length]
        left_tail_qual = read_qual[0:read_tail_length]

        right_tail_seq  = rev_comp(read_seq[-read_tail_length:])
        right_tail_qual = ''.join(reversed(read_qual[-read_tail_length:]))

        fq_left_tail_fp.write(line1)
        fq_left_tail_fp.write(left_tail_seq + '\n')
        fq_left_tail_fp.write(line3)
        fq_left_tail_fp.write(left_tail_qual + '\n')

        fq_right_tail_fp.write(line1)
        fq_right_tail_fp.write(right_tail_seq + '\n')
        fq_right_tail_fp.write(line3)
        fq_right_tail_fp.write(right_tail_qual + '\n')

    in_fastq_fp.close()
    fq_left_tail_fp.close()
    fq_right_tail_fp.close()

    return




def rev_comp(seq):

    complement  = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n': 'n'}
    rev_seq = ''.join(reversed(seq))

    rev_comp_seq = ''
    for b in rev_seq:
        rev_comp_seq += complement[b]

    return rev_comp_seq



def format_parameters(input_para):

    input_para = input_para.strip('`')
    input_para = input_para.strip('\'')
    input_para = input_para.strip('\"')
    input_para = input_para.lower()

    return input_para

def compute_overlap_len(start1, end1, start2, end2):
    max_start = max(start1, start2)
    min_end   = min(end1, end2)

    overlap_len = min_end - max_start
    return max(0, overlap_len)



### SAM/BAM/Alignment ###

def minimap2_align(in_fastq_file, ref_fasta_file, minimap2, para, out_paf):

    cmd = '%s %s %s %s > %s 2> /dev/null' % (minimap2, para, ref_fasta_file, in_fastq_file, out_paf)
    ret = os.system(cmd)
    if ret != 0: 
        eprint('ERROR: Failed to run command: %s' % cmd)
        eprint('Return value is: %d' % ret)
        sys.exit()

    return

def analysis_cigar_string (cigar):
    
    opr_set   = set(['=', 'X', 'I', 'D', 'N', 'S', 'H', 'P', 'M'])
    digit_set = set(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9'])

    cigar_opr_list = list()
    cigar_opr_len_list = list()
    
    length_str = ''
    for i in range(0, len(cigar)):
        c = cigar[i]
        if c in digit_set:
            length_str += c
        elif c in opr_set:
            cigar_opr_list.append(c)
            cigar_opr_len_list.append(int(length_str))
            length_str = ''
        else:
            eprint('ERROR: unknown CIGAR operation: %s' % c)
            sys.exit()
    
    return cigar_opr_list, cigar_opr_len_list
    
### other ###

def run_system_cmd(cmd):

    eprint('CMD: %s' % cmd)
    ret = os.system(cmd)
    if ret != 0: 
        eprint('ERROR: Failed to run command: %s' % cmd)
        eprint('Return value is: %d' % ret)
        sys.exit()
    return 
