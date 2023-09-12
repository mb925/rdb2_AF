import os
import pandas as pd
import plotly.express as px
import argparse
import matplotlib.pyplot as plt
import numpy as np
import glob

def plot_table_dataset():
    df_all = pd.DataFrame()
    for file in os.listdir(args.in_files):
        if 'summary' in file:
            df = pd.read_csv(args.in_files + '/' + file)
            df_all = pd.concat([df_all, df])
    df_all.to_csv(args.out + '/filtered_neg_datasets.csv', index=False)

def covered_residues():
    # residui ripetuti
    # residui ripetuti su totale
    # è una info che si può ottenere dai summary files
    total_res = 0
    trp_res = 0

    for file in os.listdir(args.in_files):
        if os.path.isfile(args.in_files + '/' + file):
            df = pd.read_csv(args.in_files + '/' + file)
            #for organism sum residues:
            trp_res += df['trp_residues'].iloc[0]
            total_res += df['total_residues'].iloc[0]

    perc_trp_residues = round(trp_res/total_res * 100, 2)
    # print(perc_trp_residues)
    df_residues = pd.DataFrame({'total trp residues': [trp_res], '% trp residues': [perc_trp_residues]})
    df_residues.to_csv(args.in_files + '/covered_residues.csv', index=False)



def piecharts_proteins():
    df = pd.read_csv(args.out + '/filtered_neg_datasets.csv')
    total_trp = df['total trp'].sum()
    total_prt = df['total prt'].sum()
    perc = round((total_trp * 100) / total_prt)
    diff_perc = 100 - perc
    print(perc)
    print(diff_perc)

    # trp_curated = 3480
    # total_curated = trp_curated + 1821
    # perc_cur = round((trp_curated * 100) / total_curated)
    # diff_perc = 100 - perc_cur

    fig, ax = plt.subplots()

    group_names = ['TRP organisms ' + str(perc) + '%', 'not TRP organisms ' + str(100 - perc) + '%']
    size = 0.3
    vals = np.array([[perc, perc], [diff_perc, diff_perc]]) # percentuali

    cmap = plt.colormaps["tab20c"]
    outer_colors = cmap(np.arange(4) * 4)

    ax.pie(vals.sum(axis=1), radius=1, colors=outer_colors,labels=group_names,
           wedgeprops=dict(width=size, edgecolor='w'), labeldistance=0.5)

    ax.set(aspect="equal", title='Pie plot proteins')
    plt.savefig(args.out + "/piechart_neg_prt.png")

def piecharts_residues():
    df = pd.read_csv(args.out + '/filtered_neg_datasets.csv')
    total_trp = df['trp_residues'].sum()
    total_prt = df['total_residues'].sum()
    perc = round((total_trp * 100) / total_prt)
    diff_perc = 100 - perc
    print(perc)
    print(diff_perc)

    trp_curated = 3480
    total_curated = trp_curated + 1821
    perc_cur = round((trp_curated * 100) / total_curated)
    diff_perc = 100 - perc_cur

    fig, ax = plt.subplots()

    group_names = ['TRP organisms ' + str(perc) + '%', 'not TRP organisms ' + str(100 - perc) + '%']
    size = 0.3
    vals = np.array([[perc, perc], [diff_perc, diff_perc]])  # percentuali

    cmap = plt.colormaps["tab20c"]
    outer_colors = cmap(np.arange(4) * 4)

    ax.pie(vals.sum(axis=1), radius=1, colors=outer_colors, labels=group_names,
           wedgeprops=dict(width=size, edgecolor='w'), labeldistance=0.5)


    ax.set(aspect="equal", title='Pie plot residues`')
    plt.savefig(args.out + "/filtered_neg_piechart_residues.png")



def boxplots():

    df_reg = pd.read_csv(args.out + '/datasets.csv')
    fig = px.box(df_reg, y="avg_region_len")
    fig.write_image(args.out + "/avg_region_len.png")

    df_uni = pd.read_csv(args.out + '/datasets.csv')
    fig = px.box(df_uni, y="avg_unit_lenght")
    fig.write_image(args.out + "/avg_unit_lenght.png")

def topology_to_number():
    ontology = pd.read_csv(args.ontology + '/ontology.csv', sep=',')
    df_trp = pd.DataFrame()
    df_res = pd.DataFrame()
    for file in os.listdir(args.in_files + '/topology_neg/'):
        if 'v4' in file:
            continue
        organism = file.split('_')[3].split('.')[0]
        matching_files = glob.glob(args.in_files +f"/*{organism}*")
        total_trp = pd.read_csv(matching_files[0])['total prt'][0]
        total_res = pd.read_csv(matching_files[0])['total_residues'][0]
        df = pd.read_csv(args.in_files + '/topology_neg/' + file)
        for idx, row in df.iterrows():
            tp = row['topology'].split(',')[1].lstrip()
            if tp == 'Beta hairpin repeat': tp = 'Beta hairpins'

            tp_name = ontology.loc[ontology['Name'].str.contains(tp)]['Code'].iloc[0]
            df['organism'] = organism
            if 'res' in file:
                df.loc[idx, 'trp_res'] = round(row['trp_res'] / total_res, 2)
            else:
                df.loc[idx, 'trp'] = round(row['trp'] / total_trp,2)

            df.loc[idx, 'topology'] = tp_name

        if 'res' in file:
            df_res = pd.concat([df_res, df])
        else:
            df_trp = pd.concat([df_trp, df])

    df_trp[['c', 't']] = df_trp['topology'].str.split('.', expand=True)
    df_trp["t"] = pd.to_numeric(df_trp["t"])
    df_trp = df_trp.sort_values(['organism', 'c', 't'])
    df_trp = df_trp.pivot_table(values='trp', index='organism', columns='topology', aggfunc='first')
    df_trp.to_csv(args.out + '/all_trp_by_topology_neg.csv')
    df_trp.plot.bar(stacked=True, figsize=(8,10))
    plt.legend(bbox_to_anchor=(0, 1.01, 1, 0.2), loc="lower left", mode='expand', ncol=6)
    plt.savefig(args.out + "/all_trp_by_topology_neg.png")

    df_res[['c', 't']] = df_res['topology'].str.split('.', expand=True)
    df_res["t"] = pd.to_numeric(df_res["t"])
    df_res = df_res.sort_values(['organism', 'c', 't'])
    df_res = df_res.pivot_table(values='trp_res', index='organism', columns='topology', aggfunc='first')
    df_res.to_csv(args.out + '/all_res_by_topology_neg.csv')
    df_res.plot.bar(stacked=True, figsize=(8,10))
    plt.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left", mode='expand', ncol=6)
    plt.savefig(args.out + "/all_res_by_topology_neg.png")


def amino_frequency():
    amino = pd.read_csv(args.in_files + '/amino_comp/residues.csv', sep=',')
    amino = amino.rename(columns={'0': 'residues', '1': 'quantity'})
    fig = px.bar(amino, x='residues', y='quantity')
    fig.write_image(args.out + "/residues.png")

def pie_classes():

    classes_rdb = pd.DataFrame()
    classes_rdb['class'] = ['2','3','4','5']
    classes_rdb['total'] = ['21','864','1160','120']
    fig = px.pie(classes_rdb, values='total', names='class', title='PBD (RepeatsDB)', width=400, height=400)

    ontology = pd.read_csv(args.ontology + '/ontology.csv', sep=',')
    df_trp = pd.DataFrame()
    for file in os.listdir(args.in_files + '/topology_neg/'):
        if 'v4' in file:
            continue
        organism = file.split('_')[3].split('.')[0]
        df = pd.read_csv(args.in_files + '/topology_neg/' + file)
        for idx, row in df.iterrows():
            tp = row['topology'].split(',')[1].lstrip()
            if tp == 'Beta hairpin repeat': tp = 'Beta hairpins'
            tp_name = ontology.loc[ontology['Name'].str.contains(tp)]['Code'].iloc[0]
            df['organism'] = organism
            if 'res' not in file:
                df.loc[idx, 'trp'] += row['trp']
            df.loc[idx, 'topology'] = tp_name

        if 'res' not in file:
            df_trp = pd.concat([df_trp, df])


    df_trp[['c', 't']] = df_trp['topology'].str.split('.', expand=True)
    df_trp["t"] = pd.to_numeric(df_trp["t"])
    df_trp = df_trp.sort_values(['organism', 'c', 't'])
    df_trp = df_trp.groupby('c')['trp'].sum()
    classes_af = pd.DataFrame()
    classes_af['class'] = ['2','3','4','5']
    classes_af['total'] = ['0','510086','104092','806']
    fig_af = px.pie(classes_af, values='total', names='class', title='AlphaFold', width=400, height=400)
    fig_af.show()



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate plots')
    parser.add_argument('--ontology', '-ont', type=str, help='Path to ontology file')
    parser.add_argument('--in_files', '-i', type=str, help='Path to input files')
    parser.add_argument('--out', '-o', type=str, help='Path to output')


    # args.ontology
    # -ont data/af_4/

    # args.in_files
    # -i data/af_4/stats/filtered

    # args.out
    # -o data/af_4/plots


    # Parse arguments
    args = parser.parse_args()

    # plot_table_dataset()
    # piecharts_proteins()
    # piecharts_residues()
    # covered_residues()
    # boxplots()
    # topology table parsing and plotting histogram
    # topology_to_number()
    # amino_frequency()
    pie_classes()