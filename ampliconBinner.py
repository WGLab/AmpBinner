#!/usr/bin/env python

import os
import sys
import gzip

tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + '<barcode.fasta> <in.amplicon_seq.fasta> <in.fastq> <out_dir> <n_threads>'
argc  = 5

def main():
    if len(arg) < argc:
        print (usage)
        sys.exit()

    in_barcode_fasta_file     = os.path.abspath(arg.pop(0))
    amplicon_seq_fasta_file   = os.path.abspath(arg.pop(0))
    in_fastq_file             = os.path.abspath(arg.pop(0))
    out_dir                   = os.path.abspath(arg.pop(0))
    n_threads                 = int(arg.pop(0))
    barcode_tail_length       = 256

    demultiplex_barcoded_pcr(in_barcode_fasta_file, amplicon_seq_fasta_file, in_fastq_file, n_threads, barcode_tail_length, out_dir)


def demultiplex_barcoded_pcr(in_barcode_fasta_file, amplicon_seq_fasta_file, in_fastq_file, n_threads, barcode_tail_length, out_dir):

    read_tail_add_length = 500

    barcode_name_list, barcode_seq_list = read_fasta_file(in_barcode_fasta_file)
    print('number of barcodes = %d' % len(barcode_seq_list))

    if len(barcode_name_list) == 0:
        print('ERROR!! FILE %s is empty!' % in_barcode_fasta_file)
        sys.exit()

    barcode_length = len(barcode_seq_list[0])

    fasta_name_list, fasta_seq_list = read_fasta_file(amplicon_seq_fasta_file)
    if len(fasta_name_list) > 0  and len(fasta_seq_list) > 0:
        amplicon_name = fasta_name_list[0]
        amplicon_seq = fasta_seq_list[0]
    else:
        print('ERROR!! FILE %s is empty!' % amplicon_seq_fasta_file)
        sys.exit()

    left_barcode_plus_seq_file  = os.path.join(out_dir, 'left_barcode_plus%dbp.fasta' % barcode_tail_length)
    right_barcode_plus_seq_file = os.path.join(out_dir, 'right_barcode_plus%dbp.fasta' % barcode_tail_length)

    barcode_plus_seq_to_barcode_idx_dict = generate_barcode_plus_tail_file(barcode_name_list, barcode_seq_list, amplicon_seq, barcode_tail_length, left_barcode_plus_seq_file, right_barcode_plus_seq_file)
    
    min_read_length = int( (len(amplicon_seq) + barcode_length * 2) * 0.667 )
    max_read_length = int( (len(amplicon_seq) + barcode_length * 2) * 1.5   )

    print('NOTICE: length of amplicon is %d' % len(amplicon_seq))
    print('NOTICE: reads shorter than %d bp would be skipped' % min_read_length)
    print('NOTICE: reads longer  than %d bp would be skipped' % max_read_length)

    read_tail_length = barcode_length + barcode_tail_length + read_tail_add_length

    in_fastq_file_prefix    = os.path.splitext(os.path.split(in_fastq_file)[1])[0]
    left_tail_fastq_file    = os.path.join(out_dir, '%s.left%dbp_tail.fastq'  % (in_fastq_file_prefix, read_tail_length))
    right_tail_fastq_file   = os.path.join(out_dir, '%s.right%dbp_tail.fastq' % (in_fastq_file_prefix, read_tail_length))

    extract_fastq_tail_seq (in_fastq_file, read_tail_length, min_read_length, max_read_length, left_tail_fastq_file, right_tail_fastq_file)

    left_tail_sam_file  = os.path.join(out_dir, '%s.left%dbp_tail.sam'   % (in_fastq_file_prefix, read_tail_length))
    right_tail_sam_file = os.path.join(out_dir, '%s.right%dbp_tail.sam'  % (in_fastq_file_prefix, read_tail_length))

    cmd = 'minimap2 -N 100 --cs -t %d -a -x map-ont %s %s > %s' % (n_threads, left_barcode_plus_seq_file, left_tail_fastq_file, left_tail_sam_file)
    print('running command: %s' % cmd)
    os.system(cmd)

    cmd = 'minimap2 -N 100 --cs -t %d -a -x map-ont %s %s > %s' % (n_threads, right_barcode_plus_seq_file, right_tail_fastq_file, right_tail_sam_file)
    print('running command: %s' % cmd)
    os.system(cmd)


    left_tail_sam_read_barcode_idx_dict  = extract_confident_reads_from_sam (left_tail_sam_file, barcode_length, barcode_plus_seq_to_barcode_idx_dict)

    right_tail_sam_read_barcode_idx_dict = extract_confident_reads_from_sam (right_tail_sam_file, barcode_length, barcode_plus_seq_to_barcode_idx_dict)

    read_to_sample_id_file = os.path.join(out_dir, '%s.%dbp_tail.read_to_sample_id.txt'   % (in_fastq_file_prefix, read_tail_length))
    barcode_key_to_sample_id_dict, sample_id_to_barcode_key_dict, readname_to_sample_id_dict, each_sample_read_count_list = output_read_to_sample_id_file (barcode_name_list, left_tail_sam_read_barcode_idx_dict, right_tail_sam_read_barcode_idx_dict, read_to_sample_id_file)

    
    for sample_id in range(0, len(each_sample_read_count_list)):
        barcode_key = sample_id_to_barcode_key_dict[sample_id]
        left_barcode_idx, right_barcode_idx = barcode_key.split('\t')
        left_barcode_idx = int(left_barcode_idx)
        right_barcode_idx = int(right_barcode_idx)
        left_barcode_name = barcode_name_list[left_barcode_idx]
        right_barcode_name = barcode_name_list[right_barcode_idx]
        print('sample_id = %d\tleft_barcode = %s\tright_barcode = %s\tnum_reads = %d' % (sample_id, left_barcode_name, right_barcode_name, each_sample_read_count_list[sample_id]))

    output_each_sample_fastq(in_fastq_file, readname_to_sample_id_dict, sample_id_to_barcode_key_dict, each_sample_read_count_list, barcode_name_list, out_dir)
    

    return

def output_each_sample_fastq(in_fastq_file, readname_to_sample_id_dict, sample_id_to_barcode_key_dict, each_sample_read_count_list, barcode_name_list, out_dir):

    in_fastq_file_prefix = os.path.splitext(os.path.split(in_fastq_file)[1])[0]
    out_fastq_fp_list = list()
    out_fastq_file_list = list()
    for sample_id in range(0, len(each_sample_read_count_list)):
        barcode_key = sample_id_to_barcode_key_dict[sample_id]
        left_barcode_idx, right_barcode_idx = barcode_key.split('\t')
        left_barcode_idx = int(left_barcode_idx)
        right_barcode_idx = int(right_barcode_idx)
        left_barcode_name = barcode_name_list[left_barcode_idx]
        right_barcode_name = barcode_name_list[right_barcode_idx]
        out_fastq_file = os.path.join(out_dir, '%s.demultiplexed.sample.%d.%s.%s.fastq'   % (in_fastq_file_prefix, sample_id, left_barcode_name, right_barcode_name))

        out_fastq_fp   = open(out_fastq_file, 'w')
        out_fastq_fp_list.append(out_fastq_fp)
        out_fastq_file_list.append(out_fastq_file)

    if '.gz' == in_fastq_file[-3:]:
        in_fastq_fp = gzip.open(in_fastq_file, 'rt')
    else:
        in_fastq_fp = open(in_fastq_file, 'rt')

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
        if readname not in readname_to_sample_id_dict: continue

        sample_id = readname_to_sample_id_dict[readname]
        out_fastq_fp_list[sample_id].write(line1 + line2 + line3 + line4)

    for sample_id in range(0, len(out_fastq_fp_list)):
        out_fastq_fp_list[sample_id].close()

    return

def output_read_to_sample_id_file (barcode_name_list, left_tail_sam_read_barcode_idx_dict, right_tail_sam_read_barcode_idx_dict, read_to_sample_id_file):

    barcode_key_to_sample_id_dict = dict()
    sample_id_to_barcode_key_dict = dict()
    readname_to_sample_id_dict = dict()

    sample_id = 0
    for i in range(0, len(barcode_name_list)):
        for j in range(0, len(barcode_name_list)):
            barcode_key = '%d\t%d' % (i, j)
            barcode_key_to_sample_id_dict[barcode_key] = sample_id
            sample_id_to_barcode_key_dict[sample_id] = barcode_key
            sample_id += 1

    n_samples = sample_id

    read_to_barcode_fp = open(read_to_sample_id_file, 'w')

    for readname in left_tail_sam_read_barcode_idx_dict:
        if readname in right_tail_sam_read_barcode_idx_dict:
            left_barcode_idx  = left_tail_sam_read_barcode_idx_dict[readname]
            right_barcode_idx = right_tail_sam_read_barcode_idx_dict[readname]
            barcode_key = '%d\t%d' % (left_barcode_idx, right_barcode_idx)
            sample_id = barcode_key_to_sample_id_dict[barcode_key]
            readname_to_sample_id_dict[readname] = sample_id
            read_to_barcode_fp.write('%s\t%d\t%d\t%d\n' % (readname, sample_id, left_barcode_idx, right_barcode_idx))

    read_to_barcode_fp.close()

    each_sample_read_count_list = [0] * n_samples

    for readname in readname_to_sample_id_dict:
        sample_id = readname_to_sample_id_dict[readname]
        each_sample_read_count_list[sample_id] += 1


    return barcode_key_to_sample_id_dict, sample_id_to_barcode_key_dict, readname_to_sample_id_dict, each_sample_read_count_list


def generate_barcode_plus_tail_file(barcode_name_list, barcode_seq_list, amplicon_seq, tail_length, left_barcode_plus_seq_file, right_barcode_plus_seq_file):

    barcode_plus_seq_to_barcode_idx_dict = dict()

    revc_amplicon_seq  = rev_comp(amplicon_seq)

    left_amplicon_seq  = amplicon_seq[0:tail_length]
    revc_right_amplicon_seq = revc_amplicon_seq[0:tail_length]

    left_barcode_plus_seq_list  = list()
    right_barcode_plus_seq_list = list()

    for i in range(0, len(barcode_seq_list)):
        barcode_seq  = barcode_seq_list[i]
        barcode_name = barcode_name_list[i]

        left_barcode_plus_seq = barcode_seq + left_amplicon_seq
        left_barcode_plus_seq_list.append(left_barcode_plus_seq)
        right_barcode_plus_seq = barcode_seq + revc_right_amplicon_seq
        right_barcode_plus_seq_list.append(right_barcode_plus_seq)


    left_barcode_plus_seq_fp  = open(left_barcode_plus_seq_file, 'w')
    right_barcode_plus_seq_fp = open(right_barcode_plus_seq_file, 'w')

    for i in range(0, len(barcode_seq_list)):
        left_barcode_plus_seq = left_barcode_plus_seq_list[i]
        left_barcode_plus_seq_name = '%s.plus_%dbp_left_amplicon' % (barcode_name_list[i], tail_length)
        left_barcode_plus_seq_fp.write('>%s\n' % (left_barcode_plus_seq_name))
        left_barcode_plus_seq_fp.write('%s\n' % left_barcode_plus_seq)

        right_barcode_plus_seq = right_barcode_plus_seq_list[i]
        right_barcode_plus_seq_name = '%s.plus_%dbp_right_amplicon' % (barcode_name_list[i], tail_length)
        right_barcode_plus_seq_fp.write('>%s\n' % right_barcode_plus_seq_name)
        right_barcode_plus_seq_fp.write('%s\n' % right_barcode_plus_seq)

        barcode_plus_seq_to_barcode_idx_dict[left_barcode_plus_seq_name] = i
        barcode_plus_seq_to_barcode_idx_dict[right_barcode_plus_seq_name] = i

    left_barcode_plus_seq_fp.close()
    right_barcode_plus_seq_fp.close()

    return barcode_plus_seq_to_barcode_idx_dict

def extract_confident_reads_from_sam (in_sam_file, barcode_length, barcode_plus_seq_to_barcode_idx_dict):

    read_barcode_idx_dict = dict()
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
            print('ERROR!! unknown template name in sam: %s' % contig)
            num_error_alignments += 1
            continue

        read_barcode_idx_dict[readname] = barcode_idx

    in_sam_fp.close()

    print('STATISTICS: sam_file = %s, num_aligned_reads = %d, num_unmapped_reads = %d, num_of_confident_reads = %d' % (in_sam_file, num_aligned_reads, num_unmapped_reads, len(read_barcode_idx_dict) ))

    return read_barcode_idx_dict

def extract_fastq_tail_seq(in_fastq_file, read_tail_length, min_read_length, max_read_length, left_tail_fastq_file, right_tail_fastq_file):

    if '.gz' == in_fastq_file[-3:]:
        in_fastq_fp = gzip.open(in_fastq_file, 'rt')
    else:
        in_fastq_fp = open(in_fastq_file, 'rt')


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

        num_processd_reads += 1

        if num_processd_reads % 100000 == 0:
            print('processed %d reads' % num_processd_reads)

    in_fastq_fp.close()
    fq_left_tail_fp.close()
    fq_right_tail_fp.close()

    print('NOTICE: finished extracting tail sequences from fastq. number of skipped reads = %d' % num_skipped_reads)
    return


def rev_comp(seq):

    complement  = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
    rev_seq = ''.join(reversed(seq))

    rev_comp_seq = ''
    for b in rev_seq:
        rev_comp_seq += complement[b]

    return rev_comp_seq

def read_fasta_file(fasta_file):

    if '.gz' == fasta_file[-3:]:
        fasta_fp = gzip.open(fasta_file, 'rt')
    else:
        fasta_fp = open(fasta_file, 'rt')

    fasta_fp = open(fasta_file, 'r')

    fasta_name_list = list()
    fasta_seq_list = list()
    curr_name = ''
    curr_seq = ''

    while 1:
        line = fasta_fp.readline()
        if not line: break
        line = line.strip()
        if line[0] ==  '>':
            if len(curr_seq) > 0 and len(curr_name) > 0:
                fasta_name_list.append(curr_name)
                fasta_seq_list.append(curr_seq)
            curr_name = line[1:]
            curr_seq = ''
            continue

        curr_seq += line
     
    fasta_fp.close()

    if len(curr_seq) > 0 and len(curr_name) > 0:
        fasta_name_list.append(curr_name)
        fasta_seq_list.append(curr_seq)
    
    return fasta_name_list, fasta_seq_list


if __name__ == '__main__':
    main()
