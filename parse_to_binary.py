import config as cfg
import pandas as pd
import os
import json
from Bio.PDB import *
import gzip
from statistics import mean
import fnmatch
import argparse

def parse_predicted2_to_dict(f):
    print('predicted dict')
    dict_predicted = {}
    data = pd.read_csv(f, sep=',')

    for i in data.iterrows():
        if i[1][2] == 'no regions':
            continue
        units = json.loads(i[1][2])
        for unit in units:
            unit.extend([i[1][5], i[1][6]])
        if i[1][0] in dict_predicted:
            dict_predicted[i[1][0]].extend(units)
        if i[1][0] not in dict_predicted:
            dict_predicted[i[1][0]] = []
            dict_predicted[i[1][0]].extend(units)
    return dict_predicted

def dict_to_binary(d_pred2, filename):
    df = pd.DataFrame()

    for uniprot in d_pred2:

        df_residues = pd.DataFrame()
        uniprot_id = uniprot[:-1]
        parser = MMCIFParser()
        f = cfg.data['dataset'] + filename.split('.')[0] + '/AF-' + uniprot_id + '-F1-model_v4.cif.gz'
        structure = parser.get_structure(uniprot, gzip.open(f,"rt"))

        # Using CA-CA
        ppb = CaPPBuilder()
        seq = ''
        for pp in ppb.build_peptides(structure):
            seq = pp.get_sequence()

        df_residues['uniprot_sequence'] = seq
        df_residues['residue_id'] = range(1, 1+len(df_residues))
        df_residues['RDB2'] = 0
        df_residues['REGION'] = 0
        df_residues['uniprot_id'] = uniprot


        for interval2 in d_pred2[uniprot]:
            fill_repeated_interval(uniprot, df_residues, interval2)

        df = pd.concat([df, df_residues])
        print(filename)
        print(uniprot_id)

    df.to_csv(cfg.data['binary'] + filename, index=False)

def fill_repeated_interval(uniprot, residues, interval):
    start_seqres = residues.loc[residues['residue_id'] == interval[0]]
    end_seqres = residues.loc[residues['residue_id'] == interval[1]]
    if not start_seqres.empty:
        start_seqres = start_seqres['residue_id'].iloc[0]
        if not end_seqres.empty:
            end_seqres = end_seqres['residue_id'].iloc[0]

            # scrive la regione in una colonna
            residues.loc[(residues['residue_id'] >= start_seqres) &
                         (residues['residue_id'] <= end_seqres), 'REGION'] = uniprot + '_' + str(interval[0]) + '_' + str(
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


def count_delta(x):
    uniprot_average = []
    x = x.drop_duplicates(subset=['uniprot_id', 'REGION'])
    for region in x.iterrows():
        end = int(region[1]['REGION'].split('_')[2])
        st = int(region[1]['REGION'].split('_')[1])
        delta = end - st
        uniprot_average.append(delta)
    uniprot_average = mean(uniprot_average)
    return uniprot_average


def get_data_table(f_binary,f_dataset):
    data_binary = pd.read_csv(f_binary, sep=',')
    uniprot_count_total = len(fnmatch.filter(os.listdir(f_dataset), '*.*'))

    uniprot_count_repeated = len(data_binary['uniprot_id'].drop_duplicates())
    regions = data_binary.loc[data_binary['REGION'] != '0']
    df_uniprots = regions.groupby(['uniprot_id']).apply(lambda x: count_delta(x)).reset_index()
    df_uniprots.rename({0: 'delta'}, axis='columns', inplace=True)
    organism_average = df_uniprots['delta'].mean()
    df = pd.DataFrame(columns=['uniprot_count_repeated', 'uniprot_count_total','organism_average'])
    df.loc[0] = [str(uniprot_count_repeated), uniprot_count_total, round(organism_average, 2)]

    df.to_csv(cfg.data['dataset_table'] + filename.split('.')[0] + '.csv', index=False)



if __name__ == '__main__':



    # Define argument parser
    parser = argparse.ArgumentParser(description='Generate alignments between units in the same PDB structure')
    parser.add_argument('--in_prediction', '-i', type=str, help='Path to input predictions')
    parser.add_argument('--in_dataset', '-db', type=str, help='Path to input AF dataset')
    parser.add_argument('--out', '-o', type=str, help='Path to output')

    # args.in_prediction
    # /mnt/projects/repeatsdb/prediction/rdbl_2/srul_20220829/af_4/results/UP000000429_85962_HELPY_v4.csv.gz
    # /mnt/projects/repeatsdb/prediction/rdbl_2/srul_20220829/af_4/results/UP000000579_71421_HAEIN_v4.csv.gz

    # args.in_dataset
    # /mnt/db/af/UP000000429_85962_HELPY_v4
    # /mnt/db/af/UP000000579_71421_HAEIN_v4

    # args.out
    # absolute + 'data/af_4/binary/

    # Parse arguments
    args = parser.parse_args()
    print(args.in_dataset)



    f = os.path.join(args.in_prediction)
    dict_predicted2 = parse_predicted2_to_dict(f)
    dict_to_binary(dict_predicted2, args.in_prediction) # get only file name from args.in_prediction


    f_binary = os.path.join(args.out)
    f_dataset = os.path.join(args.in_dataset)
    get_data_table(f_binary, f_dataset)

