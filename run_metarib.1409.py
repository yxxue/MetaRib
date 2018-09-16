#!/usr/bin/python
# Author: Yaxin Xue, yaxin.xue@uib.no
import ConfigParser
import os, sys, re, shutil
import random
import numpy
import pandas as pd

# source activate py27
# source deactivate
# [BASE]
# PROJECT_DIR : /Users/yaxin/Google_Share/Yaxin_PhD/IARR_project/codes/run_metarib
# DATA_DIR : /Users/yaxin/Google_Share/Yaxin_PhD/IARR_project/codes/run_metarib/data
# SAMPLING_NUM : 1000
# SEED: 100
# THREADS : 30
#
# [EMIRGE]
# EM_PATH : ~/bin/EMIRGE/emirge_amplicon.py
# EM_PARA : --phred33 -l 125 -i 250 -s 50 -a 20 -n 10
# EM_REF : /export/jonassenfs/yaxinx/meta_rna/test_v2/reference/94_otus_16S.fixed.fasta
# EM_BT : /export/jonassenfs/yaxinx/meta_rna/test_v2/reference/94_otus_16S.fixed.fasta.btindex
#
# [BBMAP]
# BM_PATH: /export/jonassenfs/yaxinx/local/software/bbmap
# MAP_PARA : minid=0.96 maxindel=1 minhits=2 idfilter=0.98
# CLS_PARA : fo=t ow=t c=t mcs=1 e=5 mid=99
def parse_cfg(config):
    # BASE
    global DATA_DIR, PROJECT_DIR, SAMPLING_NUM, THREAD
    DATA_DIR = config.get('BASE', 'DATA_DIR')
    PROJECT_DIR = config.get('BASE', 'PROJECT_DIR')
    SAMPLING_NUM = config.get('BASE', 'SAMPLING_NUM')
    #SEED = config.getint('BASE', 'SEED')
    THREAD = config.getint('BASE','THREAD')
    # EMIRGE
    global EM_PATH, EM_PARA, EM_REF, EM_BT
    EM_PATH = config.get('EMIRGE', 'EM_PATH')
    EM_PARA = config.get('EMIRGE', 'EM_PARA')
    EM_REF = config.get('EMIRGE', 'EM_REF')
    EM_BT = config.get('EMIRGE', 'EM_BT')
    # BBMAP
    global BM_PATH, MAP_PARA, CLS_PARA
    BM_PATH = config.get('BBMAP', 'BM_PATH')
    MAP_PARA = config.get('BBMAP', 'MAP_PARA')
    CLS_PARA = config.get('BBMAP', 'CLS_PARA')
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
    return(samples_list, samples_fq1_path, samples_fq1_path, all_fq1, all_fq2)

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

def subsampling_reads(unmap_fq1, unmap_fq2):
    curr_dir = os.getcwd()
    sub_fq1 = curr_dir+'/sub.1.fq'
    sub_fq2 = curr_dir+'/sub.2.fq'
    # sampleseed=SEED
    #reformat.sh in1=sub.1.rf.fq in2=sub.2.rf.fq out1=t1.1.fq out2=t1.2.fq sample=60000 ow=t reads=1000000
    sampling_num = int(SAMPLING_NUM)
    max_reads = 100 * sampling_num
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
    #cmd = ' '.join([BM_PATH+'/rename.sh', 'in='+cur_sort_fa, 'out='+cur_uniq_fa])
    os.system(cmd)
    #step4: dedup fasta
    all_dedup_fa = work_dir+'/all.dedup.fasta'
    all_dup_fa = work_dir+'/all.dup.fasta'
    cmd = ' '.join([BM_PATH+'/dedup.sh', 'in='+cur_uniq_fa, 'out='+all_dedup_fa, 'outd='+all_dup_fa, CLS_PARA, '2> dedup.log'])
    os.system(cmd)
    return(all_dedup_fa)

def run_align_bbmap(current_iter_fa, unmap_fq1, unmap_fq2):
    ref = current_iter_fa
    # build index
    #cmd = ' '.join([BM_PATH+'/bbmap.sh', 'ref='+ref])
    #os.system(cmd)
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
    EM_PARA, '-f', EM_REF, '-b', EM_BT, '2> iter_'+str(iter_time)+'_emirge.log'])
    print(cmd)
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
    # stop iteration if no new contig comes
    if os.stat(iter_fa).st_size == 0:
        keep_running = 0
        new_iter_time = iter_time + 1
        return(unmap_fq1, unmap_fq2, all_dedup_fa, new_iter_time, keep_running)
    # check if it has novel contigs
    #all_dedup_fa = dedup_contig(merged_fa, iter_fa)
    prev_fa_num = cal_fa_num(dedup_fa)
    cur_fa_num = cal_fa_num(all_dedup_fa)
    new_fa_num = cur_fa_num - prev_fa_num
    # case1: less than 20 new contigs, only keep case2 at the moment
    #if new_fa_num <= 20:
    #    keep_running = 0
    # aligning to extract unmaped reads
    (new_unmap_fq1, new_unmap_fq2) = run_align_bbmap(all_dedup_fa, unmap_fq1, unmap_fq2)
    # case2: less than 1% novel contigs
    new_unmap_fq_num = cal_fastq_num(new_unmap_fq1)
    old_unmap_fq_num = cal_fastq_num(unmap_fq1)
    curr_unmap_fq_num = old_unmap_fq_num - new_unmap_fq_num
    if curr_unmap_fq_num <= (0.01*new_unmap_fq_num):
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
        return(all_dedup_fa)

    # case2: no new contig or only a few new contigs
    if (keep_running == 0):
        num_unmap_fq = cal_fastq_num(unmap_fq1)
        if (num_unmap_fq <= 2.0*float(SAMPLING_NUM)):
            # use all unmaped reads
            sub_fq1, sub_fq2 = unmap_fq1, unmap_fq2
            # run emrige_amp and dedup
            all_dedup_fa, iter_fa = run_emirge_and_dedup(sub_fq1, sub_fq2, dedup_fa, iter_time)
            return(all_dedup_fa)
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
            return(all_dedup_fa)



def cal_abundance(samples_names, samples_pe1, samples_pe2, all_dedup_fasta):

    os.chdir(mapping_dir)
    # building the reference
    mapping_ref = mapping_dir+'all_dedup.map.fasta'

    cmd = ' '.join([BBMAP+'/bbmap.sh', 'ref='+mapping_ref])
    os.system(cmd)
    all_rpc_path = []
    for idx, val in enumerate(samples_names):
        os.mkdir(mapping_dir+val)
        os.chdir(mapping_dir+val)
        reads1 = samples_pe1[idx]
        reads2 = samples_pe2[idx]
        sam_file = val+'.derep.sam'
        # mapping
        cmd = ' '.join([BBMAP+'/bbmap.sh', 'in1='+reads1, 'in2='+reads2,'ref='+mapping_ref,'threads='+THREADS, MAP_PARA,'outu='+val+'.bbmap.unmap.fq', 'scafstats='+val+'.ctg.scafs.txt','covstats='+val+'.ctg.covs.txt'])
        os.system(cmd)
        scaf_file = os.getcwd()+'/'+val+'.ctg.scafs.txt'
        all_rpc_path.append(scaf_file)
    return(all_rpc_path)


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
    # calculate coverage for each samples
    for idx, val in enumerate(samples_list):
        sample_idx = idx
        sample_name = str(val)
        reads1 = samples_fq1_path[sample_name]
        reads2 = samples_fq2_path[sample_name]
        statsfile = sample_name+'.statsfile.txt'
        scafstats = sample_name+'.scafstats.txt'
        covstats = sample_name+'.covstats.txt'
        # run bbmap alignment, but only we only need statistics file
        cmd = ' '.join([BM_PATH+'/bbmap.sh', 'in1='+reads1, 'in2='+reads2, 'ref='+dedup_ref,
        'threads='+str(THREAD), MAP_PARA, 'ow=t', 'statsfile='+statsfile,
        'sortscafs=t', 'scafstats='+scafstats, 'covstats='+covstats, '2> run.'+sample_name+'.log'])
        scafstats = os.getcwd()+'/'+scafstats
        os.system(cmd)
        # sava scafsats file path
        all_scafstats_path[sample_name] = scafstats
    return(all_scafstats_path, dedup_ref)

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

def generate_abundance_table(samples_list, all_scafstats_path, dedup_ref):
    ab_dir = PROJECT_DIR+'/MetaRib/Abundance/'
    os.chdir(ab_dir)
    fa_ids = parse_fa_ids(dedup_ref)
    all_ab_df =['Contig_ID'] = fa_ids
    all_ab_df.set_index('Contig_ID')
    for sample_name in samples_list:
        sample_df = pd.read_csv(all_scafstats_path[sample_name], sep = '\t')
        saved_df = pd.DataFrame()
        saved_df['Contig_ID'] = sample_df['#name']
        saved_df['Contig_ID'] = save_df.Contig_ID.astype(str)
        saved_df.set_index('Contig_ID')
        # extract the % of amb and unamb reads, and sum as the abundance
        per_unamb = sample_df['%unambiguousReads']
        per_amb = sample_df['%ambiguousReads']
        per_all = per_unamb + per_amb
        saved_df[sample_name+'_estab'] = per_all
        # merge two df based on the all_ab column keys
        all_ab_df = all_ab_df.merge(saved_df, 'left')
    # save abundance file
    all_ab_file = os.getcwd()+'/all.deup.est.ab.txt'
    all_ab_df.to_csv(all_ab_file, sep='\t', header=True, index = False, float_format='%.5f', na_rep='NaN')
    os.chdir(PROJECT_DIR)
    return(all_ab_file)

def main():
    # parse config file
    config = ConfigParser.ConfigParser()
    config.read(sys.argv[1])
    run_cfg = parse_cfg(config)
    # INIT
    samples_list, samples_fq1_path, samples_fq2_path, all_fq1, all_fq2  = init(config)

    # build work folder
    #print(str(SEED),str(THREAD))
    work_dir = PROJECT_DIR+'/MetaRib'
    if not os.path.isdir(work_dir):
        os.mkdir(work_dir)
    else:
        shutil.rmtree(work_dir)
        os.mkdir(work_dir)
    os.chdir(work_dir)
    dedup_fa = os.getcwd()+'/dedup_contigs.fasta'
    open(dedup_fa,'w').close()
    global LOG
    log_file = os.getcwd()+'/run.log'
    LOG = open(log_file, 'a')
    # iteration parameters
    unmap_fq1, unmap_fq2, iter_time = all_fq1, all_fq2, 0
    keep_running = 1
    max_iter = 10 # set maximum 20 iterations
    keep_running = 1
    iteration_dir = work_dir+'/Iteration/'
    if not os.path.isdir(iteration_dir):
        os.mkdir(iteration_dir)
    reads_num = cal_fastq_num(unmap_fq1)
    #print(sub_fq1, sub_fq2, all_fq2, BM_PATH)
    # main iterative steps by while loop
    while (reads_num >= 1.1 * float(SAMPLING_NUM) and keep_running == 1 and iter_time <= max_iter):
        print('====START RUN ITERATION '+str(iter_time)+'====')
        (unmap_fq1, unmap_fq2, dedup_fa, iter_time, keep_running) = run_iteration(unmap_fq1, unmap_fq2, dedup_fa, iter_time, keep_running)
    # run last iteration, increase the size
    print('====START RUN ITERATION '+str(iter_time)+'====')
    all_dedup_fa = run_last_iteration(unmap_fq1, unmap_fq2, dedup_fa, iter_time, keep_running)
    # calculate mapping stats for each sample
    all_scafstats_path, dedup_ref = cal_mapping_stats(samples_list, samples_fq1_path, samples_fq2_path, all_dedup_fa)
    # generate abundance table based on scafstats
    all_ab_file = generate_abundance_table(samples_list, all_scafstats_path, dedup_ref)
    os.remove(dedup_fa)
    print('====PROGRAM FINISHED!====')
    return()

if __name__ == "__main__":
   main()

# python final.metarib.v1.py MetaRib.configure.ini
