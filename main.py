import config as cfg
import pandas as pd
import os
import json
from pymongo import MongoClient
from src.blast import Blast
import argparse

def parse_predicted2_to_dict(f):
    print('predicted dict')
    dict_predicted = {}
    data = pd.read_csv(f, sep=',')

    for i in data.iterrows():
        if i[1][2] == 'no regions':
            continue
        units = json.loads(i[1][2])
        if i[1][0] in dict_predicted:
            dict_predicted[i[1][0]].extend(units)
        if i[1][0] not in dict_predicted:
            dict_predicted[i[1][0]] = []
            dict_predicted[i[1][0]].extend(units)
    return dict_predicted


def dict_to_binary(d_pred2, filename):

    df = pd.DataFrame()

    for uniprot in d_pred2:


        residues = df_residues[['uniprot_id', 'uniprot_sequence']]
        residues['RDB2'] = 0
        residues['REGION'] = 0

        for interval2 in d_pred2[uniprot]:
            fill_repeated_interval(uniprot, residues, interval2)


        df = pd.concat([df, residues])

    df.to_csv(cfg.data['binary'] + filename, index=False)

def fill_repeated_interval(pdb, residues, interval):
    start_seqres = residues.loc[residues['residue_id'] == str(interval[0])]
    end_seqres = residues.loc[residues['residue_id'] == str(interval[1])]
    if not start_seqres.empty:
        start_seqres = start_seqres['residue_id'].iloc[0]
        if not end_seqres.empty:
            end_seqres = end_seqres['residue_id'].iloc[0]

            # scrivere la regione in una colonna
            residues.loc[(residues['residue_id'] >= start_seqres) &
                         (residues['residue_id'] <= end_seqres), 'REGION'] = pdb + '_' + str(interval[0]) + '_' + str(
                interval[1])
            # mettere 1 in una colonna quando ripetuto
            residues.loc[(residues['residue_id'] >= start_seqres) &
                         (residues['residue_id'] <= end_seqres), 'RDB2'] = 1
    return residues

def find_seq(db, uniprot_id):

    print('Current uniprot is %s' % str(uniprot_id))
    l = list(db.uniprot_20220823.find({'uniprot_id': uniprot_id}))
    return l


if __name__ == '__main__':

    for filename in os.listdir(cfg.data['af']):
        f = os.path.join(cfg.data['af'], filename)
        dict_predicted2 = parse_predicted2_to_dict(f)
        dict_to_binary(dict_predicted2, filename)
