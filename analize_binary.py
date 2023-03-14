import os
import pandas as pd
from statistics import mean
import argparse


def count(x):
    file_names = list(x['files'])
    return file_names

def datasets():
    mapping = pd.read_csv(args.data)
    list_files = mapping.groupby('organism').apply(lambda x: count(x)).reset_index()
    df = pd.DataFrame()
    df_organisms = []
    df_count_proteins = []
    df_count_trp_proteins = []
    df_prt_average_length = []
    # each organism
    for el in list_files.iterrows():
        count_proteins = 0
        count_trp_proteins = 0
        prt_length = []
        prt_average_length = 0
        files = el[1][0]
        organism = el[1]['organism']
        for file in files:
            possible_paths = mapping.loc[mapping['files'] == file]
            for path in possible_paths.iterrows():

                try_path = args.in_files + '/' + path[1]['organism'] + '/table_' + path[1]['files'].split('.cif')[0] + '.csv'
                print(try_path)

                if os.path.isfile(try_path):
                    af_prt = pd.read_csv(try_path)
                    print(af_prt)
                    count_proteins += 1
                    trp = af_prt['trp'].iloc[0]
                    print(trp)
                    if trp == 1:
                        count_trp_proteins += 1
                        prt_length.append(af_prt['prot_len'])
                        print(count_trp_proteins)

    # for organism in os.listdir(args.in_files):
    #     print(organism)
        # for file in os.listdir(args.in_files + '/' + organism):
        #     trp_residues = pd.read_csv(args.in_files + '/' + organism + '/table_' + file.split('.gz')[0])['trp_residues'].iloc[0]
        #     print(trp_residues)




if __name__ == '__main__':
    # Define argument parser
    parser = argparse.ArgumentParser(description='Generate statistics')
    parser.add_argument('--in_files', '-i', type=str, help='Path to input files')
    parser.add_argument('--data', '-d', type=str, help='Path to mapping files')
    parser.add_argument('--out', '-o', type=str, help='Path to output')
    parser.add_argument("-ll", "--log-level", type=str, choices=["debug", "info"], default="info",
                        help="set log level for logging")

    # args.in_files
    # /mnt/ldap/marbev/rdb-af/output/binary

    # args.out
    # data/af_4/stats


    # Parse arguments
    args = parser.parse_args()

    datasets()
