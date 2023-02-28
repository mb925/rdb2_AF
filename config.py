import os

absolute = os.path.abspath(os.getcwd())
absolute = absolute + '/'

data = {'af': '/mnt/projects/repeatsdb/prediction/rdbl_2/srul_20220829/af_4/results',
        'binary': absolute + 'data/af_4/binary',
        'dataset': '/mnt/db/af'
}

# Define Mongo configuration
db_host = '172.21.2.89'
db_port = 27017
db_name = 'biodbs'
db_coll = 'uniprot_20220823'
