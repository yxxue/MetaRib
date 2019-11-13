#!/usr/bin/python

# -*- coding: utf-8 -*-
__author__ = "Yaxin Xue"
__license__ = "GPL"
__version__ = "3.0"
__affiliation__ = "CBU, University of Bergen"
__email__ = "xue.ethan@gmail.com, yaxin.xue@uib.no"

# import modules

import ConfigParser, argparse
import os, sys, re, shutil
import random
import pandas as pd

class MyParser(argparse.ArgumentParser):
   def error(self, message):
      sys.stderr.write('error: Check your Parameters!\n')
      sys.stderr.write('error: %s\n' % message)
      self.print_help()
      sys.exit(2)

def parse_arg():
    example_text = '''Example:
    python2 run_MetaRib.py -cfg MetaRib.cfg
    '''
    parser = argparse.ArgumentParser(description='Constructing ribosomal genes from large scale total RNA meta-transcriptomic data\n',epilog=example_text, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser._optionals.title = 'Mandatory Arguments'
    parser.add_argument('-cfg', required=True, help='MetaRib configure file')
    return(parser)

def parse_cfg(config):
    # BASE
    global DATA_DIR, PROJECT_DIR, SAMPLING_NUM, THREAD
    DATA_DIR = config.get('BASE', 'DATA_DIR')
    PROJECT_DIR = config.get('BASE', 'PROJECT_DIR')
    SAMPLING_NUM = config.get('BASE', 'SAMPLING_NUM')
    THREAD = config.getint('BASE','THREAD')
    # EMIRGE
    global EM_PATH, EM_PARA, EM_REF, EM_BT
    EM_PATH = config.get('EMIRGE', 'EM_PATH')
    EM_PARA = config.get('EMIRGE', 'EM_PARA')
    EM_REF = config.get('EMIRGE', 'EM_REF')
    EM_BT = config.get('EMIRGE', 'EM_BT')
    # BBTOOL
    global BM_PATH, MAP_PARA, CLS_PARA
    BM_PATH = config.get('BBTOOL', 'BB_PATH')
    MAP_PARA = config.get('BBTOOL', 'MAP_PARA')
    CLS_PARA = config.get('BBTOOL', 'CLS_PARA')
    # FILTER
    global MIN_COV, MIN_PER
    MIN_COV = config.get('FILTER', 'MIN_COV')
    MIN_PER = config.get('FILTER', 'MIN_PER')
    return(1)

def init(config):
    data_dir = str(config.get('BASE', 'DATA_DIR'))
    samples_list_path = data_dir+'/samples.list.txt'
    samples_list = []
    samples_fq1_path = {}
    samples_fq2_path = {}
    all_fq1 = data_dir+'/all.1.fq'
    all_fq2 = data_dir+'/all.2.fq'
    for i in open(samples_list_path):
        sample_id = i.strip()
        samples_list.append(sample_id)
        fq1_path = data_dir+'/'+sample_id+'.1.fq'
        fq2_path = data_dir+'/'+sample_id+'.2.fq'
        samples_fq1_path[sample_id] = fq1_path
        samples_fq2_path[sample_id] = fq2_path
    return(samples_list, samples_fq1_path, samples_fq2_path, all_fq1, all_fq2)

def cal_fastq_num(fastq):
    fastq = fastq
    total_number_of_file = 0
    with open(fastq) as f:
        total_number_of_file = sum(1 for _ in f)
    number_of_fq = int(total_number_of_file)/4.0
    return (number_of_fq)

def cal_fa_num(in_fa):
    fa_file = open(in_fa, 'r')
    num_fa = 0
    for inp in fa_file:
        inp = inp.strip()
        if re.match('>', inp):
            num_fa +=1
    return(num_fa)

def parse_fa_ids(fa_file):
    fa_h = open(fa_file)
    fa_ids = []
    for inp in fa_h:
        inp = inp.strip()
        if inp.startswith('>'):
            header = inp.split('>')[1]
            fa_ids.append(str(header))
        else:
            continue
    return(fa_ids)

def subsampling_reads(unmap_fq1, unmap_fq2):
    curr_dir = os.getcwd()
    sub_fq1 = curr_dir+'/sub.1.fq'
    sub_fq2 = curr_dir+'/sub.2.fq'
    sampling_num = int(SAMPLING_NUM)
    max_reads = 1000 * sampling_num
    seeds = random.randint(1, 100)
    cmd = ' '.join([BM_PATH+'/reformat.sh','in1='+unmap_fq1, 'in2='+unmap_fq2,'out1='+sub_fq1,'out2='+sub_fq2, 'sample='+str(sampling_num), 'sampleseed='+str(seeds),'ow=t', 'reads='+str(max_reads), '2> subsample.log'])
    os.system(cmd)
    return(sub_fq1, sub_fq2)

def dedup_contig(old_fa, new_fa):
    work_dir = os.getcwd()
    #step1: join fasta
    cur_mg_fa = work_dir+'/current.merged.fasta'
    cmd = 'cat '+old_fa+ ' '+new_fa+' >'+cur_mg_fa
    os.system(cmd)
    #step2: sort by length
    cur_sort_fa = work_dir+'/current.sorted.fasta'
    cmd = ' '.join([BM_PATH+'/sortbyname.sh', 'in='+cur_mg_fa, 'out='+cur_sort_fa, 'length descending', '2> sort.log'])
    os.system(cmd)
    #step3: keep uniq id
    cur_uniq_fa = work_dir+'/current.uniqname.fasta'
    cmd = ' '.join([BM_PATH+'/reformat.sh', 'in='+cur_sort_fa, 'out='+cur_uniq_fa, 'uniquenames', '2> rename.log'])
    os.system(cmd)
    #step4: dedup fasta
    all_dedup_fa = work_dir+'/all.dedup.fasta'
    all_dup_fa = work_dir+'/all.dup.fasta'
    cmd = ' '.join([BM_PATH+'/dedupe.sh', 'in='+cur_uniq_fa, 'out='+all_dedup_fa, 'outd='+all_dup_fa, CLS_PARA, '2> dedupe.log'])
    os.system(cmd)
    return(all_dedup_fa)

def run_align_bbmap(current_iter_fa, unmap_fq1, unmap_fq2):
    ref = current_iter_fa
    # run bbmap alignment, since we may have duplicates, cannot calculate stats
    cmd = ' '.join([BM_PATH+'/bbmap.sh', 'in1='+unmap_fq1, 'in2='+unmap_fq2, 'ref='+ref,
    'threads='+str(THREAD), MAP_PARA, 'outu=bbmap.unmap.fq', 'ow=t', 'statsfile=bbmap.statsfile.txt',
    'sortscafs=t', 'scafstats=bbmap.scafstats.txt', 'covstats=bbmap.covstats.txt','2> bbmap.log'])
    os.system(cmd)
    # reformat to two fastq files
    new_unmap_fq1 = os.getcwd()+'/bbmap.unmaped.1.fq'
    new_unmap_fq2 = os.getcwd()+'/bbmap.unmaped.2.fq'
    cmd = ' '.join([BM_PATH+'/reformat.sh', 'in=bbmap.unmap.fq', 'out1='+new_unmap_fq1, 'out2='+new_unmap_fq2,
    '2> deinterleave.log'])
    os.system(cmd)
    # remove unmapped fq
    cmd = 'rm bbmap.unmap.fq'
    os.system(cmd)
    return(new_unmap_fq1, new_unmap_fq2)

def run_emirge_and_dedup(sub_fq1, sub_fq2, dedup_fa, iter_time):
    cmd = ' '.join(['python2', EM_PATH, 'emirge_amp/', '-1', sub_fq1, '-2', sub_fq2,
    EM_PARA, '-a', str(THREAD), '-f', EM_REF, '-b', EM_BT, '>> iter_'+str(iter_time)+'_emirge.log','2>> iter_'+str(iter_time)+'_emirge.log'])
   # print(cmd)
    os.system(cmd)
    # change to last iteration folder in EMIRGE
    os.chdir('emirge_amp')
    dirs = [d for d in os.listdir('.') if os.path.isdir(d)]
    last_cycle = sorted(dirs, key=lambda x: os.path.getctime(x), reverse=True)[0]
    os.chdir(last_cycle)
    iter_fa = os.getcwd()+'/'+last_cycle+'.cons.fasta'
    # check if it has novel contigs
    all_dedup_fa = dedup_contig(dedup_fa, iter_fa)
    return(all_dedup_fa, iter_fa)

def run_iteration(unmap_fq1, unmap_fq2, dedup_fa, iter_time, keep_running):
    iter_dir = '/'.join([PROJECT_DIR,'MetaRib', 'Iteration', 'iter_'+str(iter_time)])
    if not os.path.isdir(iter_dir):
        os.mkdir(iter_dir)
    os.chdir(iter_dir)
    sub_fq1, sub_fq2 = '', ''
    (sub_fq1, sub_fq2) = subsampling_reads(unmap_fq1, unmap_fq2)
    # run emirge and dedup
    all_dedup_fa, iter_fa = run_emirge_and_dedup(sub_fq1, sub_fq2, dedup_fa, iter_time)
    # EMIRGE stop: no more new contigs
    if os.stat(iter_fa).st_size == 0:
        keep_running = 0
        new_iter_time = iter_time + 1
        return(unmap_fq1, unmap_fq2, all_dedup_fa, new_iter_time, keep_running)
    # print fasta stats
    prev_fa_num = cal_fa_num(dedup_fa)
    cur_fa_num = cal_fa_num(all_dedup_fa)
    new_fa_num = cur_fa_num - prev_fa_num
    fa_stat = 'Iteration: '+str(iter_time)+'\tTotal contigs: '+str(cur_fa_num)+'\tNew contigs: '+str(new_fa_num)
    (new_unmap_fq1, new_unmap_fq2) = run_align_bbmap(all_dedup_fa, unmap_fq1, unmap_fq2)
    # calculate unmapped fq size (MB)
    new_unmap_fq_size = os.stat(new_unmap_fq1).st_size/(1024.0*1024)
    old_unmap_fq_size = os.stat(unmap_fq1).st_size/(1024.0*1024)
    iter_stat = fa_stat+'\tunmapped fastq (MB): '+str(round(new_unmap_fq_size,2))
    print(iter_stat)
    curr_unmap_fq_size = old_unmap_fq_size - new_unmap_fq_size
#    if curr_unmap_fq_size <= (0.01*new_unmap_fq_size):
#        keep_running = 0
    # case2: less than 1% novel contigs
    new_unmap_fq_num = cal_fastq_num(new_unmap_fq1)
    old_unmap_fq_num = cal_fastq_num(unmap_fq1)
    curr_unmap_fq_num = old_unmap_fq_num - new_unmap_fq_num
    if curr_unmap_fq_num <= (0.01*new_unmap_fq_num):
        keep_running = 0
    # case3: maximum iteration
    if iter_time == 10:
        keep_running = 0
    new_iter_time = iter_time + 1
    iter_dir = '/'.join([PROJECT_DIR, '/MetaRib/Iteration'])
    os.chdir(iter_dir)
    return (new_unmap_fq1, new_unmap_fq2, all_dedup_fa, new_iter_time, keep_running)

def run_last_iteration(unmap_fq1, unmap_fq2, dedup_fa, iter_time, keep_running):
    iter_dir = '/'.join([PROJECT_DIR, 'MetaRib','Iteration','iter_'+str(iter_time)+'_L'])
    if not os.path.isdir(iter_dir):
        os.mkdir(iter_dir)
    os.chdir(iter_dir)
    # case1: rest unmaped reads <= subsamping reads
    if (keep_running == 1):
        # use all unmaped reads
        sub_fq1, sub_fq2 = unmap_fq1, unmap_fq2
        # run emrige_amp and dedup
        all_dedup_fa, iter_fa = run_emirge_and_dedup(sub_fq1, sub_fq2, dedup_fa, iter_time)

    # case2 and 3: only a few new contigs or reach the last iteration
    if (keep_running == 0):
        num_unmap_fq = cal_fastq_num(unmap_fq1)
        if (num_unmap_fq <= 2.0*float(SAMPLING_NUM)):
            # use all unmaped reads
            sub_fq1, sub_fq2 = unmap_fq1, unmap_fq2
            # run emrige_amp and dedup
            all_dedup_fa, iter_fa = run_emirge_and_dedup(sub_fq1, sub_fq2, dedup_fa, iter_time)
        else:
            # increase subsampling reads * 2
            work_dir = os.getcwd()
            sub_fq1 = work_dir+'/sub.1.fq'
            sub_fq2 = work_dir+'/sub.2.fq'
            new_sampling_num = 2.0 * float(SAMPLING_NUM)
            max_reads = 100 * new_sampling_num
            cmd = ' '.join([BM_PATH+'/reformat.sh','in1='+unmap_fq1, 'in2='+unmap_fq2,'out1='+sub_fq1,'out2='+sub_fq2, 'sample='+str(new_sampling_num), 'ow=t', 'reads='+str(max_reads), '2> subsample.log'])
            os.system(cmd)
            # run emrige_amp and dedup
            all_dedup_fa, iter_fa = run_emirge_and_dedup(sub_fq1, sub_fq2, dedup_fa, iter_time)
    # print iter fasta stat
    prev_fa_num = cal_fa_num(dedup_fa)
    cur_fa_num = cal_fa_num(all_dedup_fa)
    new_fa_num = cur_fa_num - prev_fa_num
    fa_stat = 'Iteration: '+str(iter_time)+'\tTotal contigs: '+str(cur_fa_num)+'\tNew contigs: '+str(new_fa_num)
    print(fa_stat)
    return(all_dedup_fa)


def cal_mapping_stats(samples_list, samples_fq1_path, samples_fq2_path, all_dedup_fa):
    ab_dir = PROJECT_DIR+'/MetaRib/Abundance/'
    if not os.path.isdir(ab_dir):
        os.mkdir(ab_dir)
    else:
        shutil.rmtree(ab_dir)
        os.mkdir(ab_dir)
    os.chdir(ab_dir)
    dedup_ref = ab_dir+('all.dedup.fasta')
    shutil.copy(all_dedup_fa, dedup_ref)
    all_scafstats_path = {}
    all_covstats_path = {}
    # calculate coverage for each samples
    for idx, val in enumerate(samples_list):
        sample_idx = idx
        sample_name = str(val)
        reads1 = samples_fq1_path[sample_name]
        reads2 = samples_fq2_path[sample_name]
        statsfile = sample_name+'.statsfile.txt'
        scafstats = sample_name+'.scafstats.txt'
        covstats = sample_name+'.covstats.txt'
        # run bbmap alignment, but only we only need statistics file, set ozo=f to print all cov info
        cmd = ' '.join([BM_PATH+'/bbmap.sh', 'in1='+reads1, 'in2='+reads2, 'ref='+dedup_ref,
        'threads='+str(THREAD), MAP_PARA, 'ow=t', 'statsfile='+statsfile, 'nzo=f',
        'sortscafs=t', 'scafstats='+scafstats, 'covstats='+covstats, '2> run.'+sample_name+'.log'])
        scafstats = os.getcwd()+'/'+scafstats
        covstats = os.getcwd()+'/'+covstats
        os.system(cmd)
        # sava scafsats file path
        all_scafstats_path[sample_name] = scafstats
        all_covstats_path[sample_name] = covstats
    return(all_scafstats_path, all_covstats_path, dedup_ref)

def generate_and_filter_abundance_table(samples_list, all_scafstats_path, all_covstats_path, dedup_ref):
    ab_dir = PROJECT_DIR+'/MetaRib/Abundance/'
    os.chdir(ab_dir)
    fa_ids = parse_fa_ids(dedup_ref)
    all_ab_df = pd.DataFrame()
    all_ab_df['Contig_ID'] = fa_ids
    all_ab_df.set_index('Contig_ID')
    all_keeped_ctgs = []
    for sample_name in samples_list:
        sample_df = pd.read_csv(all_scafstats_path[sample_name], sep = '\t')
        saved_df = pd.DataFrame()
        saved_df['Contig_ID'] = sample_df['#name']
        saved_df['Contig_ID'] = saved_df.Contig_ID.astype(str)
        saved_df.set_index('Contig_ID')
        # extract the % of amb and unamb reads, and sum as the abundance
        per_unamb = sample_df['%unambiguousReads']
        per_amb = sample_df['%ambiguousReads']
        per_all = per_unamb + per_amb
        saved_df[sample_name+'_estab'] = per_all
        # merge two df based on the all_ab column keys
        all_ab_df = all_ab_df.merge(saved_df, 'left')
        # filter fasta by coverage information
        # parse coverage info
        cov_df = pd.read_csv(all_covstats_path[sample_name], sep= '\t')
        min_cov = float(MIN_COV)
        min_per = float(MIN_PER)
        # filter low avg fold and low covered percent
        cov_df_filter = cov_df.loc[(cov_df['Avg_fold'] >=min_cov) & (cov_df['Covered_percent'] >=min_per)]
        keeped_ids = cov_df_filter['#ID'].tolist()
        for ids in keeped_ids:
            if ids not in all_keeped_ctgs:
                all_keeped_ctgs.append(ids)
    # save filtered fasta
    filter_dedup_fa = os.getcwd()+'/all.dedup.filtered.fasta'
    keep_id_f = open('all.keeped.ids.txt','w')
    for inp in all_keeped_ctgs:
        inp = inp.strip()
        keep_id_f.write(inp+'\n')
    keep_id_f.close()
    cmd = ' '.join([BM_PATH+'/filterbyname.sh', 'in='+dedup_ref, 'names=all.keeped.ids.txt','out='+filter_dedup_fa,'include=t', '2>ft.log'])
    os.system(cmd)
    # remove id text file
    os.remove('all.keeped.ids.txt')
    # save abundance file
    # save filtered abundance file
    raw_fa_num = cal_fa_num(os.getcwd()+'/all.dedup.fasta')
    ft_fa_num = cal_fa_num(os.getcwd()+'/all.dedup.filtered.fasta')
    fa_stat = 'Raw contig:'+str(raw_fa_num)+'\t'+'Filtered contigs: '+str(ft_fa_num)
    print(fa_stat)
    filter_ab_df = all_ab_df.loc[all_ab_df['Contig_ID'].isin(all_keeped_ctgs)]
    filter_ab_file = os.getcwd()+'/all.dedup.filtered.est.ab.txt'
    filter_ab_df.to_csv(filter_ab_file, sep='\t', header=True, index = False, float_format='%.5f', na_rep='NaN')
    os.chdir(PROJECT_DIR)
    return(filter_ab_file)

def main():
    # parse config file
    parser = parse_arg()
    args = parser.parse_args()
    config_file = args.cfg
    config = ConfigParser.ConfigParser()
    config.read(config_file)
    run_cfg = parse_cfg(config)
    # INIT
    samples_list, samples_fq1_path, samples_fq2_path, all_fq1, all_fq2  = init(config)
    # build work folder
    work_dir = PROJECT_DIR+'/MetaRib'
    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)
    else:
        shutil.rmtree(work_dir)
        os.mkdir(work_dir)
    os.chdir(work_dir)
    dedup_fa = os.getcwd()+'/dedup_contigs.fasta'
    open(dedup_fa,'w').close()
    orig_dedup_fa = dedup_fa
    # iteration parameters
    unmap_fq1, unmap_fq2, iter_time = all_fq1, all_fq2, 1
    keep_running = 1
    max_iter = 10 # set maximum 10 iterations
    keep_running = 1
    iteration_dir = work_dir+'/Iteration/'
    if not os.path.isdir(iteration_dir):
        os.mkdir(iteration_dir)
    # main iterative steps by while loop
    while (keep_running == 1 and iter_time <= max_iter):
        curr_iter_time = iter_time
        print('====START ITERATION '+str(curr_iter_time)+'====')
        (unmap_fq1, unmap_fq2, dedup_fa, iter_time, keep_running) = run_iteration(unmap_fq1, unmap_fq2, dedup_fa, iter_time, keep_running)
        print('====FINISH ITERATION '+str(curr_iter_time)+'====')
    # it reaches maximum iteration
    if (curr_iter_time == max_iter):
        keeep_running = 0
    #  run last iteration
    curr_iter_time = iter_time
    print('====START LAST ITERATION '+str(curr_iter_time)+'====')
    all_dedup_fa = run_last_iteration(unmap_fq1, unmap_fq2, dedup_fa, iter_time, keep_running)
    print('====FINISH ITERATION '+str(curr_iter_time)+'====')
    print('====START POSTPROCESSING====')
    # calculate mapping stats for each sample
    all_scafstats_path, all_covstats_path, dedup_ref = cal_mapping_stats(samples_list, samples_fq1_path, samples_fq2_path, all_dedup_fa)
    # generate abundance table based on scafstats, and filter by coverage info
    (filter_ab_file)  = generate_and_filter_abundance_table(samples_list, all_scafstats_path, all_covstats_path, dedup_ref)
    print('====FINISH POSTPROCESSING====')
    # print final fasta stat
    os.remove(orig_dedup_fa)
    print('====PROGRAM FINISHED!====')
    return()

if __name__ == "__main__":
   main()
