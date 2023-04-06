import os
import pandas as pd
import numpy as np
import argparse
import gzip


def count(x):
    file_names = list(x['files'])
    return file_names


def datasets(org):
    mapping = pd.read_csv(args.data)
    list_files = mapping.groupby('organism').apply(lambda x: count(x)).reset_index()  # dare in input un organismo
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

        df_organism = pd.DataFrame(
            {'organism': [organism], 'total prt': [count_proteins], 'total trp': [count_trp_proteins],
             'not_trp': [count_proteins - count_trp_proteins],
             'trp_residues': [total_trp_res], 'total_residues': [total_res],
             'avg_prt_len': [round(avg_prt_len, 2)],
             'avg_unit_lenght': [round(avg_unit_len, 2)], 'avg_region_len': [round(avg_region_len, 2)]})

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


def trp_by_topology(folder, organism_classes_trp, organism_classes_res):

    organism = folder.split('_')[2]
    print(folder)

    organism_classes_trp[organism] = {}
    organism_classes_res[organism] = {}
    for file in os.listdir(args.in_files + folder):


        if file.endswith('.gz'):
            topologies = []
            with gzip.open(args.in_files + folder + '/' + file, 'rb') as f:
                next(f)
                for line in f:
                    l = str(line).split(',')
                    if len(l) > 6:
                        topology = l[6] + ', ' + l[7].replace("\\n'", '')
                        if len(topology) > 5:
                            topologies.append(topology)

                if 0 < len(set(topologies)) < 2: # not considering those with multiple topologies
                    res = len(topologies)
                    if topologies[0] not in organism_classes_res[organism]:
                        organism_classes_res[organism][topologies[0]] = 0
                        organism_classes_res[organism][topologies[0]] += res

                    else:
                        organism_classes_res[organism][topologies[0]] += res

                    if topologies[0] not in organism_classes_trp[organism]:
                        organism_classes_trp[organism][topologies[0]] = 0
                        organism_classes_trp[organism][topologies[0]] += 1
                    else:
                        organism_classes_trp[organism][topologies[0]] += 1
    df_trp = pd.DataFrame(organism_classes_trp[organism].items(), columns=['topology', 'trp'])
    df_res = pd.DataFrame(organism_classes_res[organism].items(), columns=['topology', 'trp_res'])

    df_trp.to_csv(args.out + '/topology/trp_by_topology_' + organism + '.csv',  index=False)
    df_res.to_csv(args.out + '/topology/res_by_topology_' + organism + '.csv', index=False)


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

    # organisms = ['mane_overlap_v4', 'UP000000579_71421_HAEIN_v4']
    # for org in organisms:

    organism_classes_res = {}
    organism_classes_trp = {}
    for folder in os.listdir(args.in_files):
        # datasets(folder)
        # covered_residues_per_prot(folder)


        trp_by_topology(folder, organism_classes_trp, organism_classes_res) # fill dicts organism_classes_res