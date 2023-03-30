#!/usr/bin/env python3

import pandas as pd
import os, sys
import json
from Bio.PDB import *
import gzip

import argparse
import logging

from statistics import mean
import fnmatch


def parse_predicted2_to_dict(dir, organism, uniprot_file):
    dict_predicted = {}
    dict_regions = {}

    logging.debug('start uniprot dict')


    with gzip.open(dir + organism + '.csv.gz') as f:
        for line in f:
            uniprot_id = str(line).split(',')[0][2:-1]
            if uniprot_file == uniprot_id:
                # print(str(line).split(',')[1])
                if str(line).split(',')[1] == 'no regions':
                    continue

                # units
                units = json.loads(str(line).split('"')[1])
                for unit in units:
                    unit.extend([str(line).split(',')[-6], str(line).split(',')[-5]]) # add class, topology
                if uniprot_id in dict_predicted:
                    dict_predicted[uniprot_id].extend(units)
                if uniprot_id not in dict_predicted:
                    dict_predicted[uniprot_id] = []
                    dict_predicted[uniprot_id].extend(units)

                # regions
                rg = str(line).split(',')[1]
                regions = [int(rg.split('-')[0]), int(rg.split('-')[1])]
                if uniprot_id in dict_regions:
                    dict_regions[uniprot_id].append(regions)
                if uniprot_id not in dict_regions:
                    dict_regions[uniprot_id] = []
                    dict_regions[uniprot_id].append(regions)

    logging.debug('end uniprot dict')
    return [dict_predicted, dict_regions]
def dict_to_binary(d_pred2, d_reg, filename, uniprot_file):
    logging.debug(f'dict_to_binary: {filename} ')

    # processed = 0
    # totals = len(d_pred2)
    # chrstot = len(str(totals))

    uniprot_id = uniprot_file.split('-')[1]
    # processed += 1
    # progress = processed / totals * 100
    # if processed % 100 == 0:
    #     logging.debug(f'  {processed: {chrstot}}/{totals} [{progress:5.1f}%] {uniprot_id} ...')

    parser = MMCIFParser()
    logging.debug(filename)
    f = args.in_dataset

    structure = parser.get_structure(uniprot_file, gzip.open(f, "rt"))

    # Using CA-CA
    ppb = CaPPBuilder()
    seq = ''
    for pp in ppb.build_peptides(structure):
        seq = pp.get_sequence()

    logging.debug(f'sequence: {seq} ')
    df_residues = pd.DataFrame({'uniprot_sequence': [*seq]})

    logging.debug(f'-----------------------------')

    df_residues['residue_id'] = range(1, 1 + len(df_residues))
    df_residues['RDB2'] = 0
    df_residues['REGION'] = '0'
    df_residues['UNIT'] = '0'
    df_residues['uniprot_id'] = uniprot_id

    logging.debug(f'Df: {df_residues} ')

    if uniprot_id in d_pred2:
        for interval2 in d_pred2[uniprot_id]:
            fill_repeated_interval(uniprot_id, df_residues, interval2, 'RDB2')
    if uniprot_id in d_reg:
        for interval3 in d_reg[uniprot_id]:
            fill_repeated_interval(uniprot_id, df_residues, interval3, 'REGION')


    # logging.debug(f'  {processed: {chrstot}}/{totals} [{progress:5.1f}%] {uniprot_id} ...')

    logging.debug(f'making dir {args.out} {filename}...')

    os.makedirs(args.out + filename, exist_ok=True)
    logging.debug(f'writing to {args.out} {filename}...')

    df_residues.to_csv(args.out + filename + '/' + uniprot_file + '.csv.gz', index=False)


def fill_repeated_interval(uniprot, residues, interval, type):

    start_seqres = residues.loc[residues['residue_id'] == interval[0]]
    end_seqres = residues.loc[residues['residue_id'] == interval[1]]
    if not start_seqres.empty:
        start_seqres = start_seqres['residue_id'].iloc[0]
        if not end_seqres.empty:
            end_seqres = end_seqres['residue_id'].iloc[0]

            # scrive la regione in una colonna
            if type == 'REGION':
                residues.loc[(residues['residue_id'] >= start_seqres) &
                         (residues['residue_id'] <= end_seqres),'REGION'] = uniprot + '_' + str(
                    interval[0]) + '_' + str(
                    interval[1])
            else:
                # scrive l'unità in una colonna
                residues.loc[(residues['residue_id'] >= start_seqres) &
                             (residues['residue_id'] <= end_seqres), 'UNIT'] = uniprot + '_' + str(
                    interval[0]) + '_' + str(
                    interval[1])
                # mette 1 in una colonna quando ripetuto
                residues.loc[(residues['residue_id'] >= start_seqres) &
                             (residues['residue_id'] <= end_seqres), 'RDB2'] = 1
                # mette classe
                residues.loc[(residues['residue_id'] >= start_seqres) &
                             (residues['residue_id'] <= end_seqres), 'CLASS'] = interval[2]
                # mette topology
                residues.loc[(residues['residue_id'] >= start_seqres) &
                             (residues['residue_id'] <= end_seqres), 'TOPOLOGY'] = interval[3]

    return residues



def find_seq(db, uniprot_id):
    print('Current uniprot is %s' % str(uniprot_id))
    l = list(db.uniprot_20220823.find({'uniprot_id': uniprot_id}))
    return l


def count_delta(x, type):
    average = 0
    x = x.loc[x[type] != '0'].drop_duplicates(subset=[type])

    if not x.empty:
        average = []
        for region in x.iterrows():
            end = int(region[1][type].split('_')[2])
            st = int(region[1][type].split('_')[1])
            delta = end - st
            average.append(delta)
        average = mean(average)
    return average

def get_data_table(filename, uniprot_file):
    data_binary = pd.read_csv(args.out + filename + '/' + uniprot_file + '.csv.gz', sep=',', dtype={'UNIT': str, 'REGION': str})


    trp = 0
    trp_res = len(data_binary.loc[data_binary['RDB2'] != 0])
    if trp_res > 0: trp = 1
    prt_len = len(data_binary)
    units = len(data_binary.loc[data_binary['UNIT'] != '0'].drop_duplicates(subset=['UNIT']))


    avg_reg = count_delta(data_binary, 'REGION')
    avg_unit = count_delta(data_binary, 'UNIT')
    # df_uniprots.rename({0: 'delta'}, axis='columns', inplace=True)
    # organism_average = df_uniprots['delta'].mean()
    df = pd.DataFrame(columns=['trp', 'trp_residues','protein_len', 'units', 'unit_avg_len', 'region_avg_len'])
    df.loc[0] = [trp, str(trp_res), str(prt_len), str(units), round(avg_unit, 2), round(avg_reg, 2)]

    df.to_csv(args.out + filename + '/table_' + uniprot_file + '.csv', index=False)


if __name__ == '__main__':
    # Define argument parser
    parser = argparse.ArgumentParser(description='Generate files with AF info')
    parser.add_argument('--in_prediction', '-i', type=str, help='Path to input predictions')
    parser.add_argument('--in_dataset', '-db', type=str, help='Path to input AF dataset')
    parser.add_argument('--out', '-o', type=str, help='Path to output')
    parser.add_argument("-ll", "--log-level", type=str, choices=["debug", "info"], default="info",
                        help="set log level for logging")

    # args.in_prediction
    # /mnt/projects/repeatsdb/prediction/rdbl_2/srul_20220829/af_4/results/

    # args.in_dataset
    # /mnt/db/af/UP000002716_300267_SHIDS_v4/AF-Q32DX3-F1-model_v4.cif.gz


    # args.out
    # /home/martina/PycharmProjects/rdb2_AF/data/af_4/

    # Parse arguments
    args = parser.parse_args()

    # Initialize logger
    logging.basicConfig(format='%(asctime)s %(levelname)-5.5s %(message)s',
                        level=logging.getLevelName(args.log_level.upper()), stream=sys.stdout)
    logging.info(f'{os.path.basename(__file__)} started')
    logging.debug(f'Arguments: {vars(args)}')

    filename = os.path.basename(args.in_dataset).split('.')[0]

    organism = args.in_dataset.split('/')
    del organism[-1]

    organism = organism[-1]
    uniprot = filename.split('AF-')[1]
    uniprot = uniprot.split('-F1')[0]

    logging.info(f'Working with {uniprot} [{organism}]')
    logging.debug('Processing dict ...')
    dict_predicted2 = parse_predicted2_to_dict(args.in_prediction, organism, uniprot)
    logging.debug(f'Processing residues ...')
    dict_to_binary(dict_predicted2[0], dict_predicted2[1], organism, filename)  # get only file name from args.in_prediction
    logging.debug(f'Calculating stats ...')
    get_data_table(organism, filename)

    logging.info(f'{os.path.basename(__file__)} on {uniprot} [{organism}] finished with 0')
    sys.exit(0)
