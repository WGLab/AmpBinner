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

import argparse
from multiprocessing import Process
import random

import tk


def parse_user_arguments():

    parser = argparse.ArgumentParser(description='A barcode demultiplexer for Oxford Nanopore long-read sequencing data with 10X Genomics Chromium barcodes')

    ### required arguments ###
    parser.add_argument('--in_fq', required = False, metavar = 'FILE', type = str, default = '', help = 'input sequencing reads in one FASTQ(.gz) file')
    parser.add_argument('--in_fq_list', required = False, metavar = 'FILE', type = str, default = '', help = 'a list file specifying all input FASTQ(.gz) files, one file per line')
    parser.add_argument('--barcode_list', required = True, metavar = 'FILE', type = str, default = '', help = 'a list file of all barcode sequences, one barcode sequence per line, no barcode name')
    parser.add_argument('--barcode_upstream_seq', required = True, metavar = 'STRING', type = str, default = '', help = 'known upstream sequence of the barcode')
    parser.add_argument('--out_prefix', required = True, metavar = 'PATH', type = str, help ='prefix of output files')

    ### optional arguments ###
    parser.add_argument('--num_threads', required = False, metavar = 'INT', type = int, default = 1, help ='number of threads (default: 1)')
    parser.add_argument('--minimap2', required = False, metavar = 'FILE', type = str, default = 'minimap2', help ='path to minimap2 (default: using environment default)')    
    parser.add_argument('--version', action='version', version='%(prog)s 0.4.0')

    input_args = parser.parse_args()

    return input_args

def main():
    
    input_args = parse_user_arguments()

    if input_args.num_threads < 1:
        tk.eprint('ERROR: --num_threads should be a positive number.')
        sys.exit()
    if input_args.in_fq == '' and input_args.in_fq_list == '':
        tk.eprint('ERROR! No input file! Both --in_fq and in_fq_list were not supplied. ')
        sys.exit()
    if input_args.in_fq != '' and input_args.in_fq_list != '':
        tk.eprint('ERROR! --in_fq and --in_fq_list should not be supplied at the same time.')
        sys.exit()

    
    if input_args.minimap2 != 'minimap2':
        tk.check_input_file_exists(input_args.minimap2)
        input_args.minimap2 = os.path.abspath(input_args.minimap2)
        
    tk.check_input_file_exists(input_args.barcode_list)
    input_args.barcode_list == os.path.abspath(input_args.barcode_list)
    input_args.out_prefix = os.path.abspath(input_args.out_prefix)

    AmpliconBinner_10X(input_args)

def AmpliconBinner_10X(input_args):
    
    tmp_out_prefix = input_args.out_prefix + '.tmp'

    num_threads   = input_args.num_threads
    minimap2      = input_args.minimap2

    barcode_info = BarcodeInfo()
    barcode_info.init(input_args.barcode_list, input_args.barcode_upstream_seq)

    raw_input_fq_list = list()
    if input_args.in_fq != '': 
        input_args.in_fq = os.path.abspath(input_args.in_fq)
        raw_input_fq_list.append(input_args.in_fq)
    if input_args.in_fq_list != '': 
        raw_input_fq_list = tk.read_list_file(input_args.in_fq_list, abspath = True)

    tk.eprint('NOTICE: preprocessing the input fastq file')
    fastq_file_list = tk.split_fastq(raw_input_fq_list, num_threads, tmp_out_prefix) # 1. split 2. remove duplicates
    
    process_list = list()

    for i in range(0, num_threads):
        p = Process(target=demultiplex1barcode, args=(i, fastq_file_list, barcode_info, minimap2, tmp_out_prefix) )
        process_list.append(p)
        

    for p in process_list:
        p.start()

    for p in process_list:
        p.join()


    merge_thread_summary_file(num_threads, input_args.out_prefix)

    
    cmd = 'rm %s*' % (tmp_out_prefix)
    ret = os.system(cmd)
    if ret != 0: 
        tk.eprint('ERROR: Failed to run command: %s' % cmd)
        tk.eprint('Return value is: %d' % ret)
        sys.exit()
    

    return

def merge_thread_summary_file(num_threads, out_prefix):
    out_file_list = list()
    out_stat_file_list = list()
    out_all_read_barcode_file_list = list()
    tmp_out_prefix = out_prefix + '.tmp'
    for i in range(0, num_threads):
        out_file      = tmp_out_prefix + '.thread%d' % i + '.demultiplexing.reads.txt'
        out_stat_file = tmp_out_prefix + '.thread%d' % i + '.demultiplexing.statistics.txt'
        out_all_read_barcode_file = tmp_out_prefix + '.thread%d' % i + '.all_reads.txt'

        out_file_list.append(out_file)
        out_stat_file_list.append(out_stat_file)
        out_all_read_barcode_file_list.append(out_all_read_barcode_file)

    final_out_file = out_prefix + '.demultiplexing.PASS.reads.txt'
    final_out_stat_file = out_prefix + '.demultiplexing.statistics.txt'
    final_all_read_barcode_file = out_prefix + '.all_reads.txt'

    header = '#readname\tbest_matched_barcode\tnum_edit_bases\tmismatch|insertion|deletion\tstrand\tsecond_best_matched_barcode\tnum_edit_bases\tmismatch|insertion|deletion\tstrand\n'
    final_out_fp = open(final_out_file, 'w')
    final_out_fp.write(header)
    final_out_fp.close()

    final_all_read_barcode_fp = open(final_all_read_barcode_file, 'w')
    final_all_read_barcode_fp.write(header)
    final_all_read_barcode_fp.close()
    
    cmd = 'cat '
    for f in out_file_list:
        cmd += ' %s ' % f
    cmd += ' >> %s' % final_out_file
    ret = os.system(cmd)
    if ret != 0: 
        tk.eprint('ERROR: Failed to run command: %s' % cmd)
        tk.eprint('Return value is: %d' % ret)
        sys.exit()

    cmd = 'cat '
    for f in out_all_read_barcode_file_list:
        cmd += ' %s ' % f
    cmd += ' > %s' % final_all_read_barcode_file
    ret = os.system(cmd)
    if ret != 0: 
        tk.eprint('ERROR: Failed to run command: %s' % cmd)
        tk.eprint('Return value is: %d' % ret)
        sys.exit()

    barcode_count_dict = dict()
    for stat_file in out_stat_file_list:
        stat_fp = open(stat_file, 'r')
        lines = list(stat_fp)
        stat_fp.close()
        for line in lines:
            col_list = line.strip().split('\t')
            barcode = col_list[0]
            count = int(col_list[1])
            if barcode not in barcode_count_dict:
                barcode_count_dict[barcode] = count
            else:
                barcode_count_dict[barcode] += count
    
    barcode_count_sorted_list = sorted(barcode_count_dict.items(), key=lambda x: x[1], reverse = True)
    final_out_stat_fp = open(final_out_stat_file, 'w')
    final_out_stat_fp.write('#cellular_barcode_seq\tnum_reads\n')
    for x in barcode_count_sorted_list:
        final_out_stat_fp.write('%s\t%d\n' % (x[0], x[1]))
    final_out_stat_fp.close()

    return

def demultiplex1barcode(thread_id, in_fastq_file_list, barcode_info, minimap2, tmp_out_prefix):

    in_fastq_file = in_fastq_file_list[thread_id]

    read_tail_length = len(barcode_info.upstream_seq + barcode_info.barcode_list[0]) + min(barcode_info.anchor_seq_len, len(barcode_info.downstream_seq))
    read_tail_length = int(read_tail_length * 1.5)

    tmp_out_prefix += '.thread%d' % thread_id
    left_tail_fastq_file    = tmp_out_prefix + '.left%dbp_tail.fastq'  % (read_tail_length)
    right_tail_fastq_file   = tmp_out_prefix + '.right%dbp_tail.fastq' % (read_tail_length)

    tk.eprint('NOTICE: (process %d) extracting tails from fastq reads' % thread_id)
    tk.extract_fastq_tail_seq (in_fastq_file, read_tail_length, left_tail_fastq_file, right_tail_fastq_file)
    
    tk.eprint('NOTICE: (process %d) locating anchors' % thread_id)
    left_tail_upstream_anchor_paf_file  = align_reads_to_anchors(thread_id, minimap2, 1, barcode_info, left_tail_fastq_file,  tmp_out_prefix)
    right_tail_upstream_anchor_paf_file = align_reads_to_anchors(thread_id, minimap2, 1, barcode_info, right_tail_fastq_file, tmp_out_prefix)

    upstream_anchor_avg_alignments = count_average_num_alignments(left_tail_upstream_anchor_paf_file)
   
    anchor_loc = 'none'
    if upstream_anchor_avg_alignments < 1.5 and len(barcode_info.upstream_seq) > 4:
        anchor_loc = 'upstream'

    if anchor_loc == 'upstream':
        demultiplex1barcode_method2(thread_id, left_tail_fastq_file, right_tail_fastq_file, minimap2, barcode_info, left_tail_upstream_anchor_paf_file, right_tail_upstream_anchor_paf_file, anchor_loc, tmp_out_prefix)
    else:
        if upstream_anchor_avg_alignments > 1.5 and len(barcode_info.upstream_seq) > 0:
            tk.eprint('WARNING: The UPSTREAM_SEQ (%s) have multiple alignments in reads! Try to supply a longer sequence!' % barcode_info.upstream_seq)
            tk.eprint('WARNING: AmpRepeat will try to demultiplex the reads without unique anchor sequence')

        demultiplex1barcode_method1(thread_id, left_tail_fastq_file, right_tail_fastq_file, minimap2, barcode_info, tmp_out_prefix)
                
    return

def demultiplex1barcode_method1(thread_id, left_tail_fastq_file, right_tail_fastq_file, minimap2, barcode_info, tmp_out_prefix):

    barcode_template_file = tmp_out_prefix + '.barcode_with_anchor.fasta' 
    
    left_tail_barcode_compare_paf  = tmp_out_prefix + '.left_tail_barcode_compare.paf'
    right_tail_barcode_compare_paf = tmp_out_prefix + '.right_tail_barcode_compare.paf'

    barcode_len = len(barcode_info.barcode_list[0])
    generate_barcode_template_file (barcode_info, barcode_template_file, 'none', 0)
    target_seq_len = len(barcode_info.upstream_seq) + barcode_len + len(barcode_info.downstream_seq)

    short_para  = ' -k 3 -w 2 -n 1 -m 10 -s 40 '
    mid_para    = ' -k 5 -w 3 -n 1 -m 10 -s 40 '
    normal_para = ' '

    general_para = ' -x map-ont -t 1 --for-only --eqx -c --cs -N 200 -K 1M '

    short_para  += general_para
    mid_para    += general_para
    normal_para += general_para
    
    if target_seq_len < 25:
        para = short_para
    elif target_seq_len < 50:
        para = mid_para
    else:
        para = normal_para
   
    tk.eprint('NOTICE: aligning reads to barcodes with anchors')
    tk.minimap2_align(left_tail_fastq_file, barcode_template_file, minimap2, para, left_tail_barcode_compare_paf)
    tk.minimap2_align(right_tail_fastq_file, barcode_template_file, minimap2, para, right_tail_barcode_compare_paf)

    tk.eprint('NOTICE: (process %d) assigning reads to barcodes' % thread_id)
    barcode_start_pos = len(barcode_info.upstream_seq)

    tk.eprint('DEBUG: barcode_start_pos = %d' % barcode_start_pos)
    read_barcode_info_dict = assign_reads_to_barcodes(thread_id, barcode_start_pos, barcode_start_pos + barcode_len, left_tail_barcode_compare_paf, right_tail_barcode_compare_paf)
    output_summary(barcode_info, read_barcode_info_dict, tmp_out_prefix)

    return

def demultiplex1barcode_method2(thread_id, left_tail_fastq_file, right_tail_fastq_file, minimap2, barcode_info, left_tail_anchor_paf_file, right_tail_anchor_paf_file, anchor_loc, tmp_out_prefix):

    tk.eprint('NOTICE: (process %d) demultiplexing using with anchor sequences' % thread_id)
    barcode_len = len(barcode_info.barcode_list[0])
    flank_len = 4
    
    left_tail_barcode_position_dict  = analysis_of_anchor_paf(left_tail_anchor_paf_file, barcode_len, anchor_loc, flank_len)
    right_tail_barcode_position_dict = analysis_of_anchor_paf(right_tail_anchor_paf_file, barcode_len, anchor_loc, flank_len)
    
    left_tail_barcode_candidate_fastq_file  = tmp_out_prefix + '.left_tail_barcode_candidate.fastq'  
    right_tail_barcode_candidate_fastq_file = tmp_out_prefix + '.right_tail_barcode_candidate.fastq'

    extract_region_from_fastq(left_tail_fastq_file, left_tail_barcode_position_dict, left_tail_barcode_candidate_fastq_file)
    extract_region_from_fastq(right_tail_fastq_file, right_tail_barcode_position_dict, right_tail_barcode_candidate_fastq_file)

    barcode_template_file = tmp_out_prefix + '.barcode_with_anchor.fasta'

    generate_barcode_template_file (barcode_info, barcode_template_file, anchor_loc, flank_len)
    
    left_tail_barcode_compare_paf  = tmp_out_prefix + '.left_tail_barcode_compare.paf' 
    right_tail_barcode_compare_paf = tmp_out_prefix + '.right_tail_barcode_compare.paf'

    tk.eprint('NOTICE: (process %d) aligning barcodes' % thread_id)

    barcode_compare_para = ' -t 1 --for-only --eqx -c --cs -N 200 -k 5 -w 3 -n 1 -m 10 -s 40 -A 4 -x map-ont '

    tk.minimap2_align(left_tail_barcode_candidate_fastq_file, barcode_template_file, minimap2, barcode_compare_para, left_tail_barcode_compare_paf)
    tk.minimap2_align(right_tail_barcode_candidate_fastq_file, barcode_template_file, minimap2, barcode_compare_para, right_tail_barcode_compare_paf)
    
    tk.eprint('NOTICE: (process %d) assigning reads to barcodes' % thread_id)
    if anchor_loc == 'upstream':
        barcode_start_pos = flank_len
    elif anchor_loc == 'downstream':
        barcode_start_pos = 0
   
    read_barcode_info_dict = assign_reads_to_barcodes(thread_id, barcode_start_pos, barcode_start_pos + barcode_len, left_tail_barcode_compare_paf, right_tail_barcode_compare_paf)
    output_summary(barcode_info, read_barcode_info_dict, tmp_out_prefix)

    return

def output_summary(barcode_info, read_barcode_info_dict, out_prefix):

    out_file = out_prefix + '.demultiplexing.reads.txt'
    all_read_barcode_file = out_prefix + '.all_reads.txt'
    all_read_barcode_fp = open(all_read_barcode_file, 'w')
    out_stat_file = out_prefix + '.demultiplexing.statistics.txt'
    out_fp = open(out_file, 'w')
    barcode_count_dict = dict()

    for readname in read_barcode_info_dict:
        align_info_list = read_barcode_info_dict[readname]
        if len(align_info_list) == 0: 
            continue
        elif len(align_info_list) == 1:
            best_align   = align_info_list[0]
            second_align = AlignmentInfo()
            second_align.barcode = 'N.A.'
            second_align.cigar = '*'
            all_read_barcode_fp.write('%s\t%s\t%s\n' % (readname, best_align.output(), second_align.output()))
            if best_align.num_edit_bases < 3:
                out_fp.write('%s\t%s\t%s\n' % (readname, best_align.output(), second_align.output()))
                if best_align.barcode not in barcode_count_dict:
                    barcode_count_dict[best_align.barcode] = 1
                else:
                    barcode_count_dict[best_align.barcode] += 1
        
        elif len(align_info_list) > 1:
            best_align = align_info_list[0]
            second_align = align_info_list[1]
            all_read_barcode_fp.write('%s\t%s\t%s\n' % (readname, best_align.output(), second_align.output()))
            if best_align.num_edit_bases < 3 and second_align.num_edit_bases - best_align.num_edit_bases > 2:
                out_fp.write('%s\t%s\t%s\n' % (readname, best_align.output(), second_align.output()))
            
                if best_align.barcode not in barcode_count_dict:
                    barcode_count_dict[best_align.barcode] = 1
                else:
                    barcode_count_dict[best_align.barcode] += 1


    out_fp.close()
    all_read_barcode_fp.close()
    out_stat_fp = open(out_stat_file, 'w')
    barcode_count_list = list()
    
    for barcode in barcode_count_dict:
        count = barcode_count_dict[barcode]
        barcode_count_list.append((barcode, count))
    
    barcode_count_list.sort(key = lambda x:x[1], reverse = True)
    for x in barcode_count_list:
        out_stat_fp.write('%s\t%d\n' % (x[0], x[1]))

    out_stat_fp.close()

    return

def generate_barcode_template_file (barcode_info, barcode_template_file, anchor_loc, flank_len):

    barcode_template_fp = open(barcode_template_file, 'w')
    for barcode in barcode_info.barcode_list:
        if anchor_loc == 'upstream':
            seq = barcode_info.upstream_seq[-flank_len:] + barcode
        elif anchor_loc == 'downstream':
            seq = barcode + barcode_info.downstream_seq[0:flank_len]
        else:
            seq = barcode_info.upstream_seq + barcode + barcode_info.downstream_seq
        
        barcode_template_fp.write('>%s\n' % barcode)
        barcode_template_fp.write(seq + '\n')
    
    barcode_template_fp.close()

    return
    
def extract_region_from_fastq(in_fastq_file, position_dict, out_file):

    in_fastq_fp = tk.gzopen(in_fastq_file)
    out_fp  = open(out_file, 'w')
 
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
        if readname not in position_dict: continue
 
        start, end = position_dict[readname]
        
        read_seq  = line2.strip()
        read_qual = line4.strip()

        out_seq  = read_seq[start:end]
        out_qual = read_qual[start:end]

        if len(out_seq) == 0 or len(out_qual) == 0: continue

        out_fp.write(line1)
        out_fp.write(out_seq + '\n')
        out_fp.write(line3)
        out_fp.write(out_qual + '\n')

    in_fastq_fp.close()
    out_fp.close()
  
    return

def analysis_of_anchor_paf(anchor_paf_file, barcode_len, anchor_loc, flank_len):

    barcode_position_dict = dict()
    anchor_paf_fp = open(anchor_paf_file, 'r')

    while 1:
        line = anchor_paf_fp.readline()
        if not line: break
        col_list = line.strip().split('\t')
        if len(col_list) < 12:
            tk.eprint('ERROR: There should be at least 12 columns in the PAF file: %s' % anchor_paf_file)
            sys.exit()

        readname     = col_list[0]
        read_len     = int(col_list[1])
        read_start   = int(col_list[2])
        read_end     = int(col_list[3])

        target_len   = int(col_list[6])
        target_start = int(col_list[7])
        target_end   = int(col_list[8])

        if anchor_loc == 'upstream' and target_end < target_len - 2: continue
        if anchor_loc == 'downstream' and target_start > 2: continue

        if anchor_loc == 'upstream':
            barcode_start = read_end - flank_len
            barcode_end = barcode_start + barcode_len + flank_len * 2 + target_len - target_end
        elif anchor_loc == 'downstream':
            barcode_end = read_start + flank_len
            barcode_start = barcode_end - barcode_len - flank_len * 2 - target_start

        if barcode_start < 0: barcode_start = 0
        if barcode_end > read_len: barcode_end = read_len

        if readname not in barcode_position_dict:
            barcode_position_dict[readname] = (barcode_start, barcode_end)
        
    anchor_paf_fp.close()

    return barcode_position_dict

def align_reads_to_anchors(thread_id, minimap2, num_threads, barcode_info, in_fastq_file, tmp_out_prefix):

    upstream_anchor_fasta_file = tmp_out_prefix + '.upstream_anchor.fasta'
    upstream_anchor_paf_file   = tmp_out_prefix + '.upstream_anchor.paf'

    short_para  = ' -k 3 -w 2 -n 1 -m 10 -s 40 -A 4 -x map-ont -t 1 --for-only --eqx -c --cs '
    mid_para    = ' -k 5 -w 3 -n 1 -m 10 -s 40 -A 4 -x map-ont -t 1 --for-only --eqx -c --cs '
    normal_para = ' -A 4 -x map-ont -t 1 --for-only --eqx -c --cs '

    if len(barcode_info.upstream_seq) > 0:
        target_seq_len = len(barcode_info.upstream_seq)
        upstream_anchor_fasta_fp = open(upstream_anchor_fasta_file, 'w')
        upstream_anchor_fasta_fp.write('>upstream_anchor\n')
        upstream_anchor_fasta_fp.write('%s\n' % (barcode_info.upstream_seq))
        upstream_anchor_fasta_fp.close()
        if target_seq_len < 30:
            para = short_para
        elif target_seq_len < 60:
            para = mid_para
        else:
            para = normal_para

        tk.minimap2_align(in_fastq_file, upstream_anchor_fasta_file, minimap2, para, upstream_anchor_paf_file)

    else:
        upstream_anchor_paf_fp = open(upstream_anchor_paf_file, 'w')
        upstream_anchor_paf_fp.close()

    
    return upstream_anchor_paf_file

def count_average_num_alignments(paf_file):

    total_num_align = 0
    readname_set = set()

    paf_fp = open(paf_file, 'r')
    while 1:
        line = paf_fp.readline()
        if not line: break
        total_num_align += 1
        col_list = line.strip().split('\t')
        readname_set.add(col_list[0])

    paf_fp.close()

    if len(readname_set) > 0:
        return float(total_num_align) / float(len(readname_set))
    else:
        return 1000000.0

def assign_reads_to_barcodes(thread_id, barcode_start_pos, barcode_end_pos, left_tail_paf_file, right_tail_paf_file):

    read_barcode_info_dict = dict()
    
    max_num_align_retain = 3
    left_tail_paf_fp = open(left_tail_paf_file, 'r')
    while 1:
        line = left_tail_paf_fp.readline()
        if not line: break
        col_list = line.strip().split('\t')
        readname = col_list[0]
        align_info = AlignmentInfo()
        try:
            align_info.target_len = int(col_list[6])
        except:
            tk.eprint('ERROR! file is: %s' % left_tail_paf_file)
            sys.exit()
        
        align_info.target_start = int(col_list[7])
        align_info.target_end = int(col_list[8])
        align_info.barcode = col_list[5]
        align_info.strand = 1
        align_info.mapq = int(col_list[11])
        for col in col_list[12:]:
            if col[0:5] == 'AS:i:':
                align_info.score = int(col[5:])
            if col[0:5] == 'cg:Z:':
                align_info.cigar = col[5:]
                if align_info.score > 0 and len(align_info.cigar) > 0: break

        align_info.calculate_barcode_mismatch(barcode_start_pos, barcode_end_pos)
        if readname not in read_barcode_info_dict:
            read_barcode_info_dict[readname] = list()
     
        read_barcode_info_dict[readname].append(align_info)
        read_barcode_info_dict[readname].sort(key = lambda align_info:align_info.total_num_edit_bases)
        read_barcode_info_dict[readname] = read_barcode_info_dict[readname][0:max_num_align_retain]
        
    left_tail_paf_fp.close()

    right_tail_paf_fp = open(right_tail_paf_file, 'r')
    while 1:
        line = right_tail_paf_fp.readline()
        if not line: break
        col_list = line.strip().split('\t')
        readname = col_list[0]
        align_info = AlignmentInfo()
        try:
            align_info.target_len = int(col_list[6])
        except:
            tk.eprint('ERROR! file is: %s' % right_tail_paf_file)
            sys.exit()
        align_info.target_start = int(col_list[7])
        align_info.target_end = int(col_list[8])
        align_info.barcode = col_list[5]
        align_info.strand = -1
        align_info.mapq = int(col_list[11])
        for col in col_list[12:]:
            if col[0:5] == 'AS:i:':
                align_info.score = int(col[5:])
            if col[0:5] == 'cg:Z:':
                align_info.cigar = col[5:]
                if align_info.score > 0 and len(align_info.cigar) > 0: break

        align_info.calculate_barcode_mismatch(barcode_start_pos, barcode_end_pos)
        if readname not in read_barcode_info_dict:
            read_barcode_info_dict[readname] = list()
     
        read_barcode_info_dict[readname].append(align_info)
        read_barcode_info_dict[readname].sort(key = lambda align_info:align_info.total_num_edit_bases)
        read_barcode_info_dict[readname] = read_barcode_info_dict[readname][0:max_num_align_retain]

    right_tail_paf_fp.close()
  
    for readname in read_barcode_info_dict:
        align_info_list = read_barcode_info_dict[readname]
        align_info_list.sort(key = lambda align_info:align_info.total_num_edit_bases)
    

    
    return read_barcode_info_dict


class BarcodeInfo:
    def __init__(self):
        self.upstream_seq = ''
        self.downstream_seq = ''
        self.barcode_list_file = ''
        self.barcode_list = list()
        self.anchor_seq_len = 256
        return

    def apply_anchor_len(self):

        if len(self.upstream_seq) > self.anchor_seq_len:
            self.upstream_seq = self.upstream_seq[-self.anchor_seq_len:]
        if len(self.downstream_seq) > self.anchor_seq_len:
            self.downstream_seq = self.downstream_seq[0:self.anchor_seq_len]
        return
    
    def init(self, barcode_list_file, barcode_upstream_seq):

        self.barcode_list_file = barcode_list_file
        self.read_barcode_list_file()
        self.upstream_seq = barcode_upstream_seq
        self.apply_anchor_len()

        return

    def read_barcode_list_file(self):
        
        tk.eprint('NOTICE: reading barcodes from BARCODE_LIST file: %s' % (self.barcode_list_file))
        self.barcode_list = list()
        barcode_fp = open(self.barcode_list_file, 'r')
        lines = list(barcode_fp)
        barcode_fp.close()
        
        for line in lines:
            if line[0] == '>': continue
            barcode = line.strip().split()[0]
            self.barcode_list.append(barcode)
        
        if len(self.barcode_list) == 0:
            tk.eprint('ERROR: No barcodes were found in the BARCODE_LIST file: %s' % self.barcode_list_file)
            sys.exit()
        
        self.barcode_list = list(set(self.barcode_list))
        tk.eprint('NOTICE: %d barcodes were found in the BARCODE_LIST file: %s' % (len(self.barcode_list), self.barcode_list_file))
        return

class AlignmentInfo:
    def __init__(self):
        self.barcode = ''
        self.target_len = 0
        self.target_start = 0
        self.target_end = 0
        self.score = 0
        self.mapq = 0
        self.cigar = ''
        self.strand = 0
        self.num_mismatch = 0
        self.num_ins = 0
        self.num_del = 0
        self.num_edit_bases = 0
        self.total_num_edit_bases = 0
        
    def output(self):
        return '%s\t%d\t%d|%d|%d\t%d' % (self.barcode, self.num_edit_bases, self.num_mismatch, self.num_ins, self.num_del, self.strand)

    '''
    def pmf(self):
        barcode_len = len(self.barcode)
        return scipy.stats.binom.pmf(self.num_edit_bases, barcode_len, 0.05)
    '''

    def calculate_total_num_mismatches(self):
        cigar_opr_list, cigar_opr_len_list = tk.analysis_cigar_string(self.cigar)
        if len(cigar_opr_list) != len(cigar_opr_len_list):
            tk.eprint('ERROR: len(cigar_opr_list) != len(cigar_opr_len_list)')
            sys.exit()
        
        self.total_num_edit_bases = 0
        for i in range(0, len(cigar_opr_list)):
            cigar_opr = cigar_opr_list[i]
            cigar_opr_len = cigar_opr_len_list[i]
            if cigar_opr == '=':  # match
                continue
            elif cigar_opr == 'X': # mismatch
                self.total_num_edit_bases += cigar_opr_len
            elif cigar_opr == 'I': # insertion
                self.total_num_edit_bases += cigar_opr_len
            elif cigar_opr == 'D': # deletion
                self.total_num_edit_bases += cigar_opr_len
            elif cigar_opr == 'S':
                continue
            else:
                tk.eprint('ERROR: unsupported cigar operation: %s' % cigar_opr)
                sys.exit()

        self.total_num_edit_bases += self.target_start
        self.total_num_edit_bases += self.target_len - self.target_end
        
        return

    def calculate_barcode_mismatch(self, barcode_start_pos, barcode_end_pos):
        cigar_opr_list, cigar_opr_len_list = tk.analysis_cigar_string(self.cigar)
        if len(cigar_opr_list) != len(cigar_opr_len_list):
            tk.eprint('ERROR: len(cigar_opr_list) != len(cigar_opr_len_list)')
            sys.exit()
        
        current_ref_pos = self.target_start

        for i in range(0, len(cigar_opr_list)):
            cigar_opr = cigar_opr_list[i]
            cigar_opr_len = cigar_opr_len_list[i]
            if cigar_opr == '=':  # match
                current_ref_pos += cigar_opr_len
            elif cigar_opr == 'X': # mismatch
                overlap_len = tk.compute_overlap_len(current_ref_pos, current_ref_pos+cigar_opr_len, barcode_start_pos, barcode_end_pos)
                if overlap_len > 0:
                    self.num_mismatch += overlap_len
                current_ref_pos += cigar_opr_len
            elif cigar_opr == 'I': # insertion
                if current_ref_pos > barcode_start_pos and current_ref_pos < barcode_end_pos - 1:
                    self.num_ins += cigar_opr_len
            elif cigar_opr == 'D': # deletion
                overlap_len = tk.compute_overlap_len(current_ref_pos, current_ref_pos+cigar_opr_len, barcode_start_pos, barcode_end_pos)
                if overlap_len > 0:
                    self.num_del += overlap_len
                current_ref_pos += cigar_opr_len
            elif cigar_opr == 'S':
                continue
            else:
                tk.eprint('ERROR: unsupported cigar operation: %s' % cigar_opr)
                sys.exit()
        if self.target_end < barcode_end_pos:
            self.num_mismatch += barcode_end_pos - self.target_end
        
        if self.target_start > barcode_start_pos:
            self.num_mismatch += self.target_start - barcode_start_pos
        
        self.num_edit_bases = self.num_ins + self.num_del + self.num_mismatch

        self.calculate_total_num_mismatches()
    
        return


if __name__ == '__main__':
    main()
