#!/usr/bin/env python

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

import argparse
from multiprocessing import Process
import random

import tk

tab  = '\t'
endl = '\n'
arg = sys.argv[1:]


def parse_user_arguments():

    parser = argparse.ArgumentParser(description='A barcode demultiplexer for Oxford Nanopore long-read sequencing data')

    ### required arguments ###
    parser.add_argument('--in_fq', required = False, metavar = 'FILE', type = str, default = '', help = 'input sequencing reads in one FASTQ(.gz) file')
    parser.add_argument('--in_fq_list', required = False, metavar = 'FILE', type = str, default = '', help = 'a list file specifying all input FASTQ(.gz) files, one file per line')
    parser.add_argument('--amp_seq_fasta', required = True, metavar = 'FILE', type = str, default = '', help = 'reference amplicon sequence in FASTA format')
    parser.add_argument('--out_dir', required = True, metavar = 'PATH', type = str, help ='output directory')
    parser.add_argument('--exp_name', required = True, metavar = 'STRING', type = str, help ='experimental name, used as prefix of output files')
    parser.add_argument('--fwd_barcode_fasta', required = False, metavar = 'FILE', type = str,  default = '', help ='barcode sequences of the forward primer (in FASTA format)')
    parser.add_argument('--rev_barcode_fasta', required = False, metavar = 'FILE', type = str,  default = '', help ='barcode sequences of the reverse primer (in FASTA format)')
    parser.add_argument('--require_two_barcodes', action='store_true', help='''require matched barcodes on both ends (default: False). Notice: this option is valid only if both '--fwd_barcode_fasta' and '--rev_barcode_fasta' are supplied.''')
    ### optional arguments ###
    parser.add_argument('--num_threads', required = False, metavar = 'INT', type = int, default = 1, help ='number of threads (default: 1)')
    parser.add_argument('--minimap2', required = False, metavar = 'FILE', type = str, default = 'minimap2', help ='path to minimap2 (default: using environment default)')    
    parser.add_argument('--version', action='version', version='%(prog)s 0.4.0')

    input_args = parser.parse_args()

    return input_args


def main():

    input_args = parse_user_arguments()

    if input_args.num_threads < 1:
        tk.eprint('ERROR: `--num_threads` should be a positive number.')
        sys.exit()
    if input_args.in_fq == '' and input_args.in_fq_list == '':
        tk.eprint('ERROR! No input file! Both `--in_fq` and in_fq_list were not supplied. ')
        sys.exit()
    if input_args.in_fq != '' and input_args.in_fq_list != '':
        tk.eprint('ERROR! `--in_fq` and `--in_fq_list` should not be supplied at the same time.')
        sys.exit()

    if input_args.fwd_barcode_fasta == '' and input_args.rev_barcode_fasta == '':
        tk.eprint('ERROR! Both `--fwd_barcode_fasta` and `--rev_barcode_fasta` are not supplied.')
        sys.exit()

    if input_args.minimap2 != 'minimap2':
        input_args.minimap2 = os.path.abspath(input_args.minimap2)

    input_args.out_dir  = os.path.abspath(input_args.out_dir)
    if '/' in input_args.exp_name or '\\' in input_args.exp_name:
        tk.eprint('''ERROR! `--exp_name` should not have special characters such as '/' or '\\'.''')
        sys.exit()


    AmpliconBinner(input_args)

    return

def AmpliconBinner(input_args):

    minimap2                 = input_args.minimap2
    amplicon_seq_fasta_file  = input_args.amp_seq_fasta
    n_threads                = input_args.num_threads

    out_prefix               = os.path.join(input_args.out_dir, input_args.exp_name)
    tmp_out_prefix           = out_prefix + '.tmp'

    tk.create_dir(input_args.out_dir)
    in_fastq_file = preprocessing_input_files(input_args.in_fq, input_args.in_fq_list, tmp_out_prefix)

    fwd_barcode_info       = BarcodeInfo()
    rev_barcode_info       = BarcodeInfo()
    fwd_only_out_dir    = os.path.join(input_args.out_dir, 'fwd_barcode_demultiplexing_results')
    rev_only_out_dir    = os.path.join(input_args.out_dir, 'rev_barcode_demultiplexing_results')

    fwd_barcode_out_prefix = os.path.join(fwd_only_out_dir, input_args.exp_name)
    rev_barcode_out_prefix = os.path.join(rev_only_out_dir, input_args.exp_name)

    if input_args.fwd_barcode_fasta != '' and input_args.rev_barcode_fasta == '':
        mode = 'fwd_only'
        tk.create_dir(fwd_only_out_dir)
    elif input_args.fwd_barcode_fasta == '' and input_args.rev_barcode_fasta != '':
        mode = 'rev_only'
        tk.create_dir(rev_only_out_dir)
    elif input_args.fwd_barcode_fasta != '' and input_args.rev_barcode_fasta != '':
        if input_args.require_two_barcodes:
            mode = 'both2'
        else:
            mode = 'both1'

    if input_args.fwd_barcode_fasta != '':
        fwd_barcode_info.init_from_file(input_args.fwd_barcode_fasta, amplicon_seq_fasta_file, 'fwd')

        demultiplex1side(fwd_barcode_info, in_fastq_file, minimap2, n_threads, out_prefix)
        if mode == 'fwd_only':
            out_fastq_file_list = output_binned_reads_for1side(in_fastq_file, fwd_barcode_info, fwd_barcode_out_prefix)
            remove_empty_out_fastq_file(out_fastq_file_list)

    if input_args.rev_barcode_fasta != '':

        rev_barcode_info.init_from_file(input_args.rev_barcode_fasta, amplicon_seq_fasta_file, 'rev')

        demultiplex1side(rev_barcode_info, in_fastq_file, minimap2, n_threads, out_prefix)
        if mode == 'rev_only':
            out_fastq_file_list = output_binned_reads_for1side(in_fastq_file, rev_barcode_info, rev_barcode_out_prefix)
            remove_empty_out_fastq_file(out_fastq_file_list)

    if mode == 'both1':
        fwd_rev_unmatch = check_fwd_rev_barcodes(fwd_barcode_info, rev_barcode_info)
        if fwd_rev_unmatch:
            tk.eprint('''ERROR! fwd barcodes and rev barcodes are different. Please supply the SAME barcode.fasta file if you have barcodes on both fwd and rev primers but only require barcode matching on one end.''')
            tk.eprint('''NOTICE: If you do have DIFFERENT barcodes on both ends and want to bin the reads for each side separately, you can run ampliconBinner twice and supply either '--fwd_barcode_fasta' or '--fwd_barcode_fasta' once a time.''')
            tk.eprint('''NOTICE: If you have barcodes on both ends and want to require barcode matching on both ends, please supply '--require_two_barcodes'. ''')
            sys.exit(1)

        both1_out_prefix = os.path.join(input_args.out_dir, input_args.exp_name)
        out_fastq_file_list = output_binned_reads_both1(in_fastq_file, fwd_barcode_info, rev_barcode_info, both1_out_prefix)
        remove_empty_out_fastq_file(out_fastq_file_list)

    elif mode == 'both2':
        both2_out_prefix = os.path.join(input_args.out_dir, input_args.exp_name)
        out_fastq_file_list = output_binned_reads_both2(in_fastq_file, fwd_barcode_info, rev_barcode_info, both2_out_prefix)
        remove_empty_out_fastq_file(out_fastq_file_list)
   
    cmd = 'rm %s*' % tmp_out_prefix
    tk.run_system_cmd(cmd)

    return

def remove_empty_out_fastq_file(out_fastq_file_list):

    for out_fastq_file in out_fastq_file_list:
        fsize = 0
        try:
            fsize = os.path.getsize(out_fastq_file)
        except:
            pass

        if fsize == 0:
            try: 
                os.remove(out_fastq_file)
            except:
                eprint('WARNING: failed to remove empty file: %s' % out_fastq_file)
    return

def output_binned_reads_both2(in_fastq_file, fwd_barcode_info, rev_barcode_info, out_prefix):


    readname_to_barcode_file = out_prefix + '.readname_to_barcode.txt' 
    summary_file             = out_prefix + '.summary.txt'

    readname_to_sample_idx_dict = dict()
    sample_read_count_dict = dict()

    out_fastq_file_list = list()
    sample_idx = -1
    barcode_key_to_sample_idx_dict = dict()
    sample_idx_to_barcode_dict = dict()
    for i in range(0, len(fwd_barcode_info.barcode_name_list)):
        for j in range(0, len(rev_barcode_info.barcode_name_list)):
            sample_name = 'fwd.%s.rev.%s' % (fwd_barcode_info.barcode_name_list[i], rev_barcode_info.barcode_name_list[j])
            sample_idx += 1
            barcode_key = '%d\t%d' % (i, j)
            barcode_key_to_sample_idx_dict[barcode_key] = sample_idx
            out_fastq_file = out_prefix + '.%s.fastq' % sample_name
            out_fastq_file_list.append(out_fastq_file)

    for readname in fwd_barcode_info.read_barcode_idx_dict:
        fwd_barcode_idx = fwd_barcode_info.read_barcode_idx_dict[readname]
        if readname not in rev_barcode_info.read_barcode_idx_dict: continue
        rev_barcode_idx = rev_barcode_info.read_barcode_idx_dict[readname]

        barcode_key = '%d\t%d' % (fwd_barcode_idx, rev_barcode_idx)
        sample_idx = barcode_key_to_sample_idx_dict[barcode_key]
        sample_idx_to_barcode_dict[sample_idx] = (fwd_barcode_idx, rev_barcode_idx)
        readname_to_sample_idx_dict[readname] = sample_idx

        if sample_idx not in sample_read_count_dict:
            sample_read_count_dict[sample_idx] = 1
        else:
            sample_read_count_dict[sample_idx] += 1

    sorted_sample_read_count_list = sorted(sample_read_count_dict.items(), key=lambda x: x[1], reverse=True)
    summary_fp = open(summary_file, 'w')
    summary_fp.write('#fwd_barcode_name\trev_barcode_name\tnum_reads\n')
    for x in sorted_sample_read_count_list:
        sample_idx = x[0]
        num_reads = x[1]
        fwd_barcode_idx, rev_barcode_idx = sample_idx_to_barcode_dict[sample_idx]
        fwd_barcode_name = fwd_barcode_info.barcode_name_list[fwd_barcode_idx]
        rev_barcode_name = rev_barcode_info.barcode_name_list[rev_barcode_idx]
        summary_fp.write('%s\t%s\t%d\n' % (fwd_barcode_name, rev_barcode_name, num_reads))
    summary_fp.close()

    readname_to_barcode_name_list = list()
    for readname in readname_to_sample_idx_dict:
        sample_idx = readname_to_sample_idx_dict[readname]
        fwd_barcode_idx, rev_barcode_idx = sample_idx_to_barcode_dict[sample_idx]
        fwd_barcode_name = fwd_barcode_info.barcode_name_list[fwd_barcode_idx]
        rev_barcode_name = rev_barcode_info.barcode_name_list[rev_barcode_idx]

        readname_to_barcode_name_list.append((sample_idx, readname, fwd_barcode_name, rev_barcode_name))
        
    readname_to_barcode_name_list.sort(key = lambda x:x[0])

    readname_to_barcode_fp = open(readname_to_barcode_file, 'w')
    readname_to_barcode_fp.write('#readname\tfwd_barcode_name\trev_barcode_name\n')

    for x in readname_to_barcode_name_list:
        readname_to_barcode_fp.write('%s\t%s\t%s\n' % (x[1], x[2], x[3]))
    readname_to_barcode_fp.close()

    output_binned_fastq(in_fastq_file, readname_to_sample_idx_dict, out_fastq_file_list)

    return out_fastq_file_list

def output_binned_reads_both1(in_fastq_file, fwd_barcode_info, rev_barcode_info, out_prefix):

    readname_to_barcode_file = out_prefix + '.readname_to_barcode.txt' 
    summary_file             = out_prefix + '.summary.txt'

    readname_to_sample_idx_dict = dict()
    barcode_read_count_dict = dict()

    discordant_readname_set = set()

    ## fwd_barcode ##
    for readname in fwd_barcode_info.read_barcode_idx_dict:
        fwd_barcode_idx = fwd_barcode_info.read_barcode_idx_dict[readname]
        if readname in rev_barcode_info.read_barcode_idx_dict:
            rev_barcode_idx = rev_barcode_info.read_barcode_idx_dict[readname]

            if fwd_barcode_idx != rev_barcode_idx: 
                discordant_readname_set.add(readname)
                continue

        if readname in readname_to_sample_idx_dict and readname_to_sample_idx_dict[readname] != fwd_barcode_idx:
            discordant_readname_set.add(readname)
            continue

        readname_to_sample_idx_dict[readname] = fwd_barcode_idx

    ## fwd_barcode ##
    for readname in rev_barcode_info.read_barcode_idx_dict:
        rev_barcode_idx = rev_barcode_info.read_barcode_idx_dict[readname]
        if readname in fwd_barcode_info.read_barcode_idx_dict: continue
        if readname in readname_to_sample_idx_dict and readname_to_sample_idx_dict[readname] != rev_barcode_idx:
            discordant_readname_set.add(readname)
            continue
        readname_to_sample_idx_dict[readname] = rev_barcode_idx

    tk.eprint('WARNING: %d reads have discordant barcodes on the two ends' % len(discordant_readname_set))

    readname_to_barcode_name_list = list()
    for readname in readname_to_sample_idx_dict:
        barcode_idx = readname_to_sample_idx_dict[readname]
        barcode_name = fwd_barcode_info.barcode_name_list[barcode_idx]
        if barcode_name not in barcode_read_count_dict:
            barcode_read_count_dict[barcode_name] = 1
        else:
            barcode_read_count_dict[barcode_name] += 1
        readname_to_barcode_name_list.append((barcode_idx, readname, barcode_name))
        
    readname_to_barcode_name_list.sort(key = lambda x:x[0])
    sorted_barcode_read_count_list = sorted(barcode_read_count_dict.items(), key=lambda x: x[1], reverse=True)

    readname_to_barcode_fp  = open(readname_to_barcode_file, 'w')
    readname_to_barcode_fp.write('#readname\tbarcode_name\n')
    for x in readname_to_barcode_name_list:
        readname_to_barcode_fp.write('%s\t%s\n' % (x[1], x[2]))
    readname_to_barcode_fp.close()

    summary_fp = open(summary_file, 'w')
    summary_fp.write('#barcode_name\tnum_reads\n')
    for x in sorted_barcode_read_count_list:
        summary_fp.write('%s\t%d\n' % (x[0], x[1]))
    summary_fp.close()


    out_fastq_file_list = list()
    
    for i in range(0, len(fwd_barcode_info.barcode_name_list)):
        out_fastq_file = out_prefix + '.%s.fastq' % fwd_barcode_info.barcode_name_list[i]
        out_fastq_file_list.append(out_fastq_file)

    output_binned_fastq(in_fastq_file, readname_to_sample_idx_dict, out_fastq_file_list)

    return out_fastq_file_list

def output_binned_fastq(in_fastq_file, readname_to_sample_idx_dict, out_fastq_file_list):

    out_fastq_fp_list = list()
    for i in range(0, len(out_fastq_file_list)):
        out_fastq_file = out_fastq_file_list[i]
        out_fastq_fp   = open(out_fastq_file, 'w')
        out_fastq_fp_list.append(out_fastq_fp)
    
    in_fastq_fp = tk.gzopen(in_fastq_file)

    while 1:
        line1 = in_fastq_fp.readline()
        line2 = in_fastq_fp.readline()
        line3 = in_fastq_fp.readline()
        line4 = in_fastq_fp.readline()

        if not line1: break
        if not line2: break
        if not line3: break
        if not line4: break

        readname = line1.strip().split()[0][1:]
        if readname not in readname_to_sample_idx_dict: continue

        sample_idx = readname_to_sample_idx_dict[readname]
        out_fastq_fp_list[sample_idx].write(line1 + line2 + line3 + line4)

    in_fastq_fp.close()
    for sample_idx in range(0, len(out_fastq_fp_list)):
        out_fastq_fp_list[sample_idx].close()

    return


def check_fwd_rev_barcodes(fwd_barcode_info, rev_barcode_info):

    fwd_barcode_list = list()
    rev_barcode_list = list()

    unmatch = False
    for i in range(0, len(fwd_barcode_info.barcode_name_list)):
        barcode_name = fwd_barcode_info.barcode_name_list[i]
        barcode_seq = fwd_barcode_info.barcode_seq_list[i]
        fwd_barcode_list.append(barcode_name + '\t' + barcode_seq)

    for i in range(0, len(rev_barcode_info.barcode_name_list)):
        barcode_name = rev_barcode_info.barcode_name_list[i]
        barcode_seq = rev_barcode_info.barcode_seq_list[i]
        rev_barcode_list.append(barcode_name + '\t' + barcode_seq)

    if len(fwd_barcode_list) != len(rev_barcode_list):
        unmatch = True

    for i in range(0, len(fwd_barcode_list)):
        if fwd_barcode_list[i] != rev_barcode_list[i]: unmatch = True

    return unmatch


def demultiplex1side(barcode_info, in_fastq_file, minimap2, n_threads, out_prefix):
    
    tmp_out_prefix  = out_prefix + '.tmp.%s' % barcode_info.side
    read_tail_add_length  = 256
    barcode_length = len(barcode_info.barcode_seq_list[0])

    barcode_plus_seq_file = tmp_out_prefix + '.barcode_plus%dbp.fasta' % barcode_info.anchor_seq_len
 
    generate_barcode_plus_tail_file(barcode_info, barcode_plus_seq_file)
    
    min_read_length = int( len(barcode_info.amplicon_seq) * 0.667 )
    max_read_length = int( (len(barcode_info.amplicon_seq) + barcode_length*2) * 1.5 )

    tk.eprint('NOTICE: length of amplicon is %d' % len(barcode_info.amplicon_seq))
    tk.eprint('NOTICE: reads shorter than %d bp would be skipped' % min_read_length)
    tk.eprint('NOTICE: reads longer  than %d bp would be skipped' % max_read_length)

    read_tail_length = barcode_length + barcode_info.anchor_seq_len + read_tail_add_length

    left_tail_fastq_file    = '%s.left%dbp_tail.fastq'  % (tmp_out_prefix, read_tail_length)
    right_tail_fastq_file   = '%s.right%dbp_tail.fastq' % (tmp_out_prefix, read_tail_length)

    extract_fastq_tail_seq (in_fastq_file, read_tail_length, min_read_length, max_read_length, left_tail_fastq_file, right_tail_fastq_file)

    left_tail_sam_file    = '%s.left%dbp_tail.sam'  % (tmp_out_prefix, read_tail_length)
    right_tail_sam_file   = '%s.right%dbp_tail.sam' % (tmp_out_prefix, read_tail_length)
   
    cmd = '%s -N 400 --cs -t %d -a -x map-ont %s %s > %s 2> /dev/null' % (minimap2, n_threads, barcode_plus_seq_file, left_tail_fastq_file, left_tail_sam_file)
    tk.run_system_cmd(cmd)

    cmd = '%s -N 400 --cs -t %d -a -x map-ont %s %s > %s 2> /dev/null' % (minimap2, n_threads, barcode_plus_seq_file, right_tail_fastq_file, right_tail_sam_file)
    tk.run_system_cmd(cmd)

    barcode_info.read_barcode_idx_dict  = dict() # read_barcode_idx_dict[readname] = barcode_idx
    extract_confident_reads_from_sam (left_tail_sam_file,  barcode_length, barcode_info.barcode_plus_seq_to_barcode_idx_dict, barcode_info.read_barcode_idx_dict)
    extract_confident_reads_from_sam (right_tail_sam_file, barcode_length, barcode_info.barcode_plus_seq_to_barcode_idx_dict, barcode_info.read_barcode_idx_dict)

    cmd = 'rm %s*' % tmp_out_prefix
    tk.run_system_cmd(cmd)
    return

def output_binned_reads_for1side(in_fastq_file, barcode_info, out_prefix):

    readname_to_barcode_file = out_prefix + '.readname_to_barcode.txt' 
    summary_file             = out_prefix + '.summary.txt'
    readname_to_barcode_name_list = list()
    barcode_read_count_dict = dict()
    for readname in barcode_info.read_barcode_idx_dict:
        barcode_idx  = barcode_info.read_barcode_idx_dict[readname]
        barcode_name = barcode_info.barcode_name_list[barcode_idx]
        if barcode_name not in barcode_read_count_dict:
            barcode_read_count_dict[barcode_name] = 1
        else:
            barcode_read_count_dict[barcode_name] += 1
        readname_to_barcode_name_list.append((barcode_idx, readname, barcode_name))
        
    readname_to_barcode_name_list.sort(key = lambda x:x[0])
    sorted_barcode_read_count_list = sorted(barcode_read_count_dict.items(), key=lambda x: x[1], reverse=True)

    readname_to_barcode_fp  = open(readname_to_barcode_file, 'w')
    readname_to_barcode_fp.write('#readname\tbarcode_name\n')
    for x in readname_to_barcode_name_list:
        readname_to_barcode_fp.write('%s\t%s\n' % (x[1], x[2]))
    readname_to_barcode_fp.close()

    summary_fp = open(summary_file, 'w')
    summary_fp.write('#barcode_name\tnum_reads\n')
    for x in sorted_barcode_read_count_list:
        summary_fp.write('%s\t%d\n' % (x[0], x[1]))
    summary_fp.close()

    out_fastq_file_list = list()
    out_fastq_fp_list = list()
    for i in range(0, len(barcode_info.barcode_name_list)):
        out_fastq_file = out_prefix + '.%s.fastq' % barcode_info.barcode_name_list[i]
        out_fastq_file_list.append(out_fastq_file)
        out_fastq_fp   = open(out_fastq_file, 'w')
        out_fastq_fp_list.append(out_fastq_fp)
    
  
    in_fastq_fp = tk.gzopen(in_fastq_file)

    while 1:
        line1 = in_fastq_fp.readline()
        line2 = in_fastq_fp.readline()
        line3 = in_fastq_fp.readline()
        line4 = in_fastq_fp.readline()

        if not line1: break
        if not line2: break
        if not line3: break
        if not line4: break

        readname = line1.strip().split()[0][1:]
        if readname not in barcode_info.read_barcode_idx_dict: continue

        out_file_id = barcode_info.read_barcode_idx_dict[readname]
        out_fastq_fp_list[out_file_id].write(line1 + line2 + line3 + line4)

    in_fastq_fp.close()
    for sample_id in range(0, len(out_fastq_fp_list)):
        out_fastq_fp_list[sample_id].close()

    return out_fastq_file_list

def preprocessing_input_files(in_fq, in_fq_list, tmp_out_prefix):

    tk.eprint('NOTICE: preprocessing the input fastq file')
    raw_input_fq_list = list()
    if in_fq != '': 
        in_fq = os.path.abspath(in_fq)
        raw_input_fq_list.append(in_fq)
    if in_fq_list != '': 
        raw_input_fq_list = tk.read_list_file(in_fq_list, abspath = True)

    fastq_file_list = tk.split_fastq(raw_input_fq_list, 1, tmp_out_prefix) # 1. split 2. remove duplicates
    in_fastq_file = fastq_file_list[0]
    return in_fastq_file


def generate_barcode_plus_tail_file(barcode_info, barcode_plus_seq_file):

    barcode_info.barcode_plus_seq_to_barcode_idx_dict = dict()
    barcode_plus_seq_fp  = open(barcode_plus_seq_file, 'w')

    for i in range(0, len(barcode_info.barcode_seq_list)):
        barcode_plus_seq      = barcode_info.barcode_seq_list[i] + barcode_info.downstream_seq
        barcode_plus_seq_name = '%s.%s.with_%dbp_downstream_seq' % (barcode_info.side, barcode_info.barcode_name_list[i], barcode_info.anchor_seq_len)
        barcode_plus_seq_fp.write('>%s\n' % (barcode_plus_seq_name))
        barcode_plus_seq_fp.write('%s\n' % barcode_plus_seq)

        barcode_info.barcode_plus_seq_to_barcode_idx_dict[barcode_plus_seq_name] = i

    barcode_plus_seq_fp.close()

    return

def extract_confident_reads_from_sam (in_sam_file, barcode_length, barcode_plus_seq_to_barcode_idx_dict, read_barcode_idx_dict):

    min_mapq = 20

    in_sam_fp = open(in_sam_file, 'r')
    num_error_alignments = 0
    num_aligned_reads = 0
    num_unmapped_reads = 0

    while 1:
        line = in_sam_fp.readline()
        if not line: break
        if line[0] == '@': continue

        line = line.strip().split('\t')
        if len(line) < 6:
            num_error_alignments += 1
            continue
        readname, flag, contig, left_pos, mapq = line[0:5]
        
        flag = int(flag)
        if flag & 4:
            num_unmapped_reads += 1
            continue

        if flag & 256 or flag & 1024 or flag & 2048: continue
        
        num_aligned_reads += 1

        mapq = int(mapq)
        if mapq < min_mapq: continue

        left_pos = int(left_pos)
        if left_pos >= barcode_length: continue
        
        if contig in barcode_plus_seq_to_barcode_idx_dict:
            barcode_idx = barcode_plus_seq_to_barcode_idx_dict[contig]
        else:
            tk.eprint('ERROR!! unknown template name in sam: %s' % contig)
            num_error_alignments += 1
            continue

        read_barcode_idx_dict[readname] = barcode_idx

    in_sam_fp.close()

    tk.eprint('STATISTICS: sam_file = %s, num_aligned_reads = %d, num_unmapped_reads = %d, num_of_confident_reads = %d' % (in_sam_file, num_aligned_reads, num_unmapped_reads, len(read_barcode_idx_dict) ))

    return 

def extract_fastq_tail_seq(in_fastq_file, read_tail_length, min_read_length, max_read_length, left_tail_fastq_file, right_tail_fastq_file):

    in_fastq_fp = tk.gzopen(in_fastq_file)

    fq_left_tail_fp  = open(left_tail_fastq_file, 'w')
    fq_right_tail_fp = open(right_tail_fastq_file, 'w')

    num_skipped_reads = 0

    num_processd_reads = 0
    while 1:
        line1 = in_fastq_fp.readline()
        line2 = in_fastq_fp.readline()
        line3 = in_fastq_fp.readline()
        line4 = in_fastq_fp.readline()

        if not line1: break
        if not line2: break
        if not line3: break
        if not line4: break

        read_seq  = line2.strip()
        if len(read_seq) < min_read_length or len(read_seq) > max_read_length:
            num_skipped_reads += 1
            continue

        read_qual = line4.strip()

        left_tail_seq  = read_seq[0:read_tail_length]
        left_tail_qual = read_qual[0:read_tail_length]

        right_tail_seq  = tk.rev_comp(read_seq[-read_tail_length:])
        right_tail_qual = ''.join(reversed(read_qual[-read_tail_length:]))

        fq_left_tail_fp.write(line1)
        fq_left_tail_fp.write(left_tail_seq + '\n')
        fq_left_tail_fp.write(line3)
        fq_left_tail_fp.write(left_tail_qual + '\n')

        fq_right_tail_fp.write(line1)
        fq_right_tail_fp.write(right_tail_seq + '\n')
        fq_right_tail_fp.write(line3)
        fq_right_tail_fp.write(right_tail_qual + '\n')

        num_processd_reads += 1

        if num_processd_reads % 100000 == 0:
            tk.eprint('processed %d reads' % num_processd_reads)

    in_fastq_fp.close()
    fq_left_tail_fp.close()
    fq_right_tail_fp.close()

    tk.eprint('NOTICE: finished extracting tail sequences from fastq. number of skipped reads = %d' % num_skipped_reads)
    return


class BarcodeInfo:
    def __init__(self):
        self.downstream_seq = ''
        self.barcode_fasta_file = ''
        self.amplicon_seq_fasta_file = ''
        self.amplicon_name = ''
        self.amplicon_seq = ''
        self.side = ''
        self.barcode_name_list = list()
        self.barcode_seq_list = list()
        self.anchor_seq_len = 256
        self.barcode_plus_seq_to_barcode_idx_dict = dict()
        self.read_barcode_idx_dict = dict()
        return
    
    def init_from_file(self, barcode_fasta_file, amplicon_seq_fasta_file, barcode_side):

        self.barcode_fasta_file = barcode_fasta_file
        self.amplicon_seq_fasta_file = amplicon_seq_fasta_file
        self.side = barcode_side

        self.barcode_name_list, self.barcode_seq_list = tk.read_fasta_file(barcode_fasta_file)
        fasta_name_list, fasta_seq_list = tk.read_fasta_file(amplicon_seq_fasta_file)
        if len(fasta_name_list) > 1:
            tk.eprint('ERROR: There are more than 1 sequence in the amp_seq_fasta file: %s' % amplicon_seq_fasta_file)
            sys.exit()
        if len(fasta_name_list) == 1 and len(fasta_seq_list) == 1:
            self.amplicon_name = fasta_name_list[0]
            self.amplicon_seq  = fasta_seq_list[0]
        else:
            tk.eprint('ERROR: No sequence was found in the amp_seq_fasta file: %s' % amplicon_seq_fasta_file)
            sys.exit()

        if barcode_side == 'fwd':
            self.downstream_seq = self.amplicon_seq[0:self.anchor_seq_len]
        elif barcode_side == 'rev':
            self.downstream_seq = tk.rev_comp(self.amplicon_seq[-self.anchor_seq_len:])
        return


if __name__ == '__main__':
    main()
