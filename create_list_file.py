#!/usr/bin/env python3

import os
import sys


def create_list(listfile, dbdir):
    with open(listfile, 'w') as the_file:
        for organism in os.listdir(dbdir):
            if organism == 'mane_overlap_v4' or organism == 'swissprot_cif_v4' or organism == 'download_metadata.json':
                continue

            the_file.write(
                '/projects/repeatsdb/prediction/rdbl_2/srul_20220829/af_4/results/' + organism + '.csv.gz' '\n')


if __name__ == '__main__':


    if len(sys.argv) != 3:
        print(sys.argv)

        print("Usage: {} path/to/list.dat /path/to/af".format(os.path.basename(__file__)))
        exit(1)
    listfile = sys.argv[1]
    dbdir = sys.argv[2]
    create_list(listfile, dbdir)
