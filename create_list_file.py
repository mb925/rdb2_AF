#!/usr/bin/env python3

import os
import sys
import pandas as pd


def create_list(listfile, dbdir):
    df = pd.DataFrame()
    paths = []
    files = []
    organisms = []
    for organism in os.listdir('/mnt' + dbdir):
        print(organism)

        if organism == 'download_metadata.json':
            continue
        for file in os.listdir('/mnt' + dbdir + '/' + organism):
            if file.endswith('.gz'):
                paths.append(dbdir + '/' + organism + '/' + file)
                files.append(file)
                organisms.append(organism)

    df['paths'] = paths
    df['files'] = files
    df['organism'] = organisms
    df.to_csv('data/mapping.csv', index=False)


    df = df.drop_duplicates(subset=['files'])

    with open(listfile, 'w') as the_file:
       for row in df.iterrows():
           the_file.write(row[1]['paths'] + '\n')




if __name__ == '__main__':


    if len(sys.argv) != 3:
        print(sys.argv)

        print("Usage: {} path/to/list.dat /path/to/af".format(os.path.basename(__file__)))
        exit(1)
    listfile = sys.argv[1]
    dbdir = sys.argv[2]
    create_list(listfile, dbdir)
