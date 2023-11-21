#!/usr/bin/env python3
##!/opt/gensoft/exe/Python/3.8.1/scripts/python3

import logging
import argparse
import subprocess
import io
import os
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
import pyranges as pr

def detect_ring(args):
    sample = args.cram.split('/')[-1][:-5] if args.sample_name == '$CRAM' else args.sample_name
    os.makedirs(os.path.dirname(args.output_dir+"/"), exist_ok=True)
    logging.basicConfig(filename=args.output_dir+'/'+sample+'.log', encoding='utf-8', level=logging.DEBUG, filemode='w', format="[%(asctime)s] {%(filename)s:%(lineno)d} %(levelname)s - \n%(message)s\n", datefmt='%H:%M:%S')
    s_params=''
    for arg in vars(args):
        s_params+=arg+': '+str(getattr(args, arg))+'\n'
    logging.info(s_params)
    samtools_command = ["samtools", "view", "-F", "2", "-q", "40", "--reference", args.fasta, args.cram, "chr22:15550986-51324926"]
    logging.debug('samtools command for discordant extraction:\n'+' '.join(samtools_command))
    disc_extraction = subprocess.run(samtools_command, stdout=subprocess.PIPE, text=True) #"chr22:15550986-51324926","chr22:47096982-47098982"

    df_disc_pre = pd.read_csv(io.StringIO(disc_extraction.stdout), sep='\t', usecols=[0,2,3,4,5,6,7], names=['read_id', 'chr', 'pos', 'mq', 'map', 'chr_mate', 'pos_mate'])
    
    logging.info('nb discordant extracted: '+str(len(df_disc_pre)))

    df_bed_blacklist = pd.read_csv(args.blacklist_bed, sep='\t', names=['Chromosome', 'Start', 'End'])
    gr_bed_blacklist = pr.PyRanges(df_bed_blacklist[['Chromosome', 'Start', 'End']])
    def isBlacklisted(row):
        start = int(row['pos'])
        end = int(row['pos']+1)
        if len( gr_bed_blacklist[row['chr'], start:end] ) > 0:
            return True
        return False
    df_disc_pre['isblacklisted'] = df_disc_pre.apply(isBlacklisted, axis=1)
    df_disc_pre = df_disc_pre[~df_disc_pre['isblacklisted']]

    df_disc_pre['chr_mate'] = df_disc_pre.apply(lambda x: x['chr'] if x['chr_mate']=='=' else x['chr_mate'], axis=1)
    df_disc_pre['clustered'] = False
    df_disc = df_disc_pre[df_disc_pre['mq']>=60]
    df_acro = pd.read_csv(args.acro_shortarm_file, sep='\t', names=['Chromosome', 'Start', 'End'])
    df_disc_acro = df_disc[df_disc['chr_mate'].isin(df_acro['Chromosome'])]
    df_disc_acro = df_disc_acro.merge(df_acro, left_on=['chr_mate'], right_on=['Chromosome']).drop(['Chromosome', 'Start'], axis=1).rename({'End':'centromere'}, axis=1)
    df_disc_acro = df_disc_acro[df_disc_acro['pos_mate']<=df_disc_acro['centromere']]
    
    if len(df_disc_acro)<=0:
        logging.info('no ring detected.')
        return
    
    X = np.array(df_disc_acro['pos'])
    clustering = DBSCAN(eps=args.dbscan_eps, min_samples=args.dbscan_min).fit(X.reshape(-1, 1))
    df_disc_acro['cluster'] = clustering.labels_
    logging.info('nb reads clustered: '+str(len(df_disc_acro[df_disc_acro['cluster']!=-1])))
    
    df_temp = df_disc_acro[df_disc_acro['cluster']!=-1].groupby(['cluster']).agg({'pos':['mean', 'min', 'max'], 'read_id':'count'}).reset_index()
    df_temp.columns = ['cluster','pos_mean', 'pos_min', 'pos_max', 'nb_reads']
    
    if len(df_temp)<=0:
        logging.info('no ring detected.')
        return
    
    for row in df_temp.itertuples():
        df_disc_pre.loc[ (df_disc_pre['pos']>=getattr(row, 'pos_min')) & (df_disc_pre['pos']<=getattr(row, 'pos_max')) , 'clustered' ] = True
    
    logging.info('Extracting mate from the cram file')
    with open(args.output_dir+"/"+sample+'.read2_for_blat.fasta', 'w') as f:
        for row in df_disc_pre[df_disc_pre['clustered']][['chr_mate', 'pos_mate', 'pos_mate', 'read_id']].itertuples():
            samtools_command = ["samtools","view", "-F", "2", "--reference", args.fasta, 
                                args.cram, getattr(row, 'chr_mate')+":"+str(getattr(row, 'pos_mate'))+"-"+str(getattr(row, 'pos_mate'))]
            mate_extraction = subprocess.run(samtools_command, stdout=subprocess.PIPE, text=True)
            df_samtools_temp = pd.read_csv(io.StringIO(mate_extraction.stdout), sep='\t', usecols=[0,9], names=['read_id', 'seq'])
            # print(df_samtools_temp[df_samtools_temp['read_id']==getattr(row, 'read_id')].to_markdown())
            for row in df_samtools_temp[df_samtools_temp['read_id']==getattr(row, 'read_id')].itertuples():
                print(f'>{row.read_id}\n{row.seq}', file=f)
        
    blat_command = ['blat', args.fasta, args.output_dir+"/"+sample+'.read2_for_blat.fasta', '-out=blast8', args.output_dir+"/"+sample+'.read2.blat.tsv']
    logging.info("running blat on the mates:\n"+' '.join(blat_command))
    blat_run = subprocess.run(blat_command, stdout=subprocess.PIPE, text=True)
    
    df_blat = pd.read_csv(args.output_dir+"/"+sample+'.read2.blat.tsv', sep='\t', names=['read_id', 'chr', 'identity', 'alignment length', 'mismatches', "gap openings", "q.start","q.end","s.start","s.end","e-value","bit score"])
    df_blat = df_blat[(df_blat['e-value']<=args.e_value) & (df_blat['bit score']>=args.bit_score)]

    # specificity score for acrocentric
    df_blat_tmp = df_blat.groupby(['read_id', 'chr']).agg({'e-value':'min','bit score':'max', 'alignment length':'max'}).reset_index()
    df_blat_tmp = df_blat_tmp.groupby(['chr']).agg({'read_id':'count','e-value':'mean','bit score':'sum'}).reset_index().sort_values(['bit score', 'e-value'], ascending=[False,True])
    acrosum = df_blat_tmp[df_blat_tmp['chr'].isin(['chr22','chr21', 'chr15', 'chr14', 'chr13'])]['bit score'].sum()
    nonacrosum = df_blat_tmp[~df_blat_tmp['chr'].isin(['chr22','chr21', 'chr15', 'chr14', 'chr13'])]['bit score'].sum()
    logging.info('specificity score for acrocentric: '+str(acrosum/nonacrosum))
    
    if len(df_blat)<=0:
        logging.info('no ring detected.')
        return

    df_bed_guarracino = pd.read_csv(args.phr_bed, sep='\t', names=['Chromosome', 'Start', 'End', 'avg Shannon Diversity Index', 'avg number of contigs'])
    if args.chr22_only:
        df_bed_guarracino = df_bed_guarracino[df_bed_guarracino['Chromosome']=='chr22']
    gr_bed_guarracino = pr.PyRanges(df_bed_guarracino[['Chromosome', 'Start', 'End']])
    def isPHR(row):
        start = int(row['s.start'])
        end = int(row['s.end'])
        if len( gr_bed_guarracino[row['chr'], start:end] ) > 0:
            return True
        return False
    df_blat['isPHR'] = df_blat.apply(isPHR, axis=1)
    logging.info('nb read mapped within PHR: '+str(len(df_blat[df_blat['isPHR']].groupby(['read_id']).count())))
    
    df_disc_blat = df_disc_pre[df_disc_pre['clustered']].merge(df_blat[df_blat['isPHR']].groupby(['read_id']).agg({'identity':'max', 'alignment length':'max', 'bit score':'max'}).reset_index(), on='read_id')
    
    if len(df_disc_blat)<=0:
        logging.info('no ring detected.')
        return
    
    X = np.array(df_disc_blat['pos'])
    clustering = DBSCAN(eps=args.dbscan_eps, min_samples=args.dbscan_min).fit(X.reshape(-1, 1))
    df_disc_blat['cluster'] = clustering.labels_
    
    df_disc_blat[df_disc_blat['cluster']!=-1].to_csv(args.output_dir+'/'+sample+'.supporting_reads.tsv', index=False, sep='\t')
    s_final = ''
    if len(set(df_disc_blat[df_disc_blat['cluster']!=-1]['cluster']))>=1 :
        s_final = 'sample\tcluster\tbreakpoint\tsupporting_reads\tacrocentric_specificity'
        for cluster in set(df_disc_blat[df_disc_blat['cluster']!=-1]['cluster']):
            pos_mean = df_disc_blat[df_disc_blat['cluster']==cluster]['pos'].mean()
            read_count = df_disc_blat[df_disc_blat['cluster']==cluster]['read_id'].count()

            df_blat_tmp = df_blat.groupby(['read_id', 'chr']).agg({'e-value':'min','bit score':'max', 'alignment length':'max'}).reset_index()
            df_blat_tmp = df_blat_tmp[df_blat_tmp['read_id'].isin(list(df_disc_blat[df_disc_blat['cluster']==cluster]['read_id']))].groupby(['chr']).agg({'read_id':'count','e-value':'mean','bit score':'sum'}).reset_index().sort_values(['bit score', 'e-value'], ascending=[False,True])
            acrosum = df_blat_tmp[df_blat_tmp['chr'].isin(['chr22','chr21', 'chr15', 'chr14', 'chr13'])]['bit score'].sum()
            nonacrosum = df_blat_tmp[~df_blat_tmp['chr'].isin(['chr22','chr21', 'chr15', 'chr14', 'chr13'])]['bit score'].sum()

            if nonacrosum==0 and acrosum==0:
                ratio = 'NA'
            elif nonacrosum==0:
                ratio = 'inf'
            else:
                ratio = str(acrosum/nonacrosum)

            s_final+= '\n'+sample+'\tcluster_'+str(cluster)+'\tchr22:'+str(int(pos_mean))+'\t'+str(read_count)+'\t'+ratio
    else:
        s_final = 'no ring detected.'
    logging.info(s_final)
    print(s_final)
    
parser = argparse.ArgumentParser()
parser.add_argument('cram', help='cram file to search for ring')
parser.add_argument('--chr', default='22', help='chromosome for which to search for a ring (need to be acrocentric)')
parser.add_argument('--sample-name', default='$CRAM', help='sample name, used for output files (default is the same as the cram name)')
parser.add_argument('--e-value', default=0.01, help='e-value threshold for blat (keeping <=)')
parser.add_argument('--bit-score', default=150, help='bit score threshold for blat (keeping >=)')
parser.add_argument('--fasta', help='path to the fasta used for the cram alignment')
parser.add_argument('--acro-shortarm-file', help='path to the bed file having the coordinates of all acrocentric chromosome short arms.')
parser.add_argument('--blacklist-bed', help='bed of blacklisted regions')
parser.add_argument('--output-dir', help='path to the output directory')
parser.add_argument('--dbscan-eps', default=100, help='EPS value used by DBSCAN inside when running the clustering of the positions (default: 100)')
parser.add_argument('--dbscan-min', default=5, help='Min number of reads inside a DBSCAN cluster (default: 5)')
parser.add_argument('--phr-bed', help='bed file for the PHR region from Guarracino et al.')
parser.add_argument('--chr22-only', action='store_true', help='if we test only for the PHRs of the 22p')

detect_ring(parser.parse_args())