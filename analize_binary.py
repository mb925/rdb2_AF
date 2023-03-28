import os
import pandas as pd
import numpy as np
import argparse


def count(x):
    file_names = list(x['files'])
    return file_names

def datasets(org):
    mapping = pd.read_csv(args.data)
    list_files = mapping.groupby('organism').apply(lambda x: count(x)).reset_index() # dare in input un organismo
    # each organism
    # list_files = list_files.loc[list_files['organism'] == args.organism]
    list_files = list_files.loc[list_files['organism'] == org]
    total_trp_res = 0
    total_res = 0
    for el in list_files.iterrows():
        count_proteins = 0
        count_trp_proteins = 0
        prt_length = []
        unit_length = []
        region_length = []
        avg_prt_len = 0
        files = el[1][0]
        avg_unit_len = 0
        avg_region_len = 0

        # per organism
        organism = el[1]['organism']
        for file in files:
            possible_paths = mapping.loc[mapping['files'] == file]
            for path in possible_paths.iterrows():

                try_path = args.in_files + path[1]['organism'] + '/table_' + path[1]['files'].split('.cif')[0] + '.csv'
                print(try_path)

                if os.path.isfile(try_path):
                    af_prt = pd.read_csv(try_path)
                    print(af_prt)
                    count_proteins += 1
                    trp = af_prt['trp'].iloc[0]

                    if trp == 1:
                        count_trp_proteins += 1
                        total_trp_res += af_prt['trp_residues'].iloc[0]
                        total_res += af_prt['protein_len'].iloc[0]
                        prt_length.append(af_prt['protein_len'].iloc[0])
                        unit_length.append(af_prt['unit_avg_len'].iloc[0])
                        region_length.append(af_prt['region_avg_len'].iloc[0])

        if len(prt_length) != 0: avg_prt_len = np.mean(prt_length)


        if len(prt_length) != 0: avg_unit_len = np.mean(unit_length)
        if len(prt_length) != 0: avg_region_len = np.mean(region_length)

        df_organism = pd.DataFrame({'organism': [organism],  'total prt': [count_proteins], 'total trp': [count_trp_proteins],
                                    'not_trp': [count_proteins - count_trp_proteins],
                                    'trp_residues': [total_trp_res], 'total_residues': [total_res],
                                     'avg_prt_len': [round(avg_prt_len,2)],
                                    'avg_unit_lenght': [round(avg_unit_len,2)], 'avg_region_len':[round(avg_region_len,2)]})

        df_organism.to_csv(args.out + '/summary_' + organism + '.csv', index=False)

def covered_residues_per_prot(org):
    mapping = pd.read_csv(args.data)
    list_files = mapping.groupby('organism').apply(lambda x: count(x)).reset_index()  # dare in input un organismo
    # each organism
    # list_files = list_files.loc[list_files['organism'] == args.organism]
    list_files = list_files.loc[list_files['organism'] == org]
    total_trp_res = 0
    total_res = 0
    for el in list_files.iterrows():
        files = el[1][0]

        # per organism
        organism = el[1]['organism']
        for file in files:
            possible_paths = mapping.loc[mapping['files'] == file]
            for path in possible_paths.iterrows():

                try_path = args.in_files + path[1]['organism'] + '/' + path[1]['files'].split('.cif')[0] + '.csv.gz'
                print(try_path)

                if os.path.isfile(try_path):
                    df = pd.read_csv(try_path)
                    print(df)


if __name__ == '__main__':
    # Define argument parser
    parser = argparse.ArgumentParser(description='Generate statistics')
    parser.add_argument('--in_files', '-i', type=str, help='Path to input files')
    parser.add_argument('--data', '-d', type=str, help='Path to mapping files')
    parser.add_argument('--organism', '-og', type=str, help='Organism folder')
    parser.add_argument('--out', '-o', type=str, help='Path to output')
    parser.add_argument("-ll", "--log-level", type=str, choices=["debug", "info"], default="info",
                        help="set log level for logging")

    # -i /mnt/ldap/marbev/rdb-af/output/
    # -og UP000020681_1299332_MYCUL_v4
    # -d data/mapping.csv
    # -o data/af_4/stats


    # Parse arguments
    args = parser.parse_args()
    total_trp_res = 0
    total_res = 0
    # organisms = ['UP000000579_71421_HAEIN_v4', 'UP000002716_300267_SHIDS_v4', 'UP000020681_1299332_MYCUL_v4',
    #              'UP000008827_3847_SOYBN_v4', 'UP000094526_86049_9EURO1_v4','UP000000589_10090_MOUSE_v4',
    #              'UP000274756_318479_DRAME_v4','UP000024404_6282_ONCVO_v4',
    #              'UP000006304_1133849_9NOCA1_v4', 'UP000030665_36087_TRITR_v4']

    organisms = ['UP000000579_71421_HAEIN_v4']
    for org in organisms:
        # datasets(org)
        covered_residues_per_prot(org)

