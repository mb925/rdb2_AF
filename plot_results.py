import os
import pandas as pd
import plotly.express as px
import argparse
import matplotlib.pyplot as plt
import numpy as np

def plot_table_dataset():
    df_all = pd.DataFrame()
    for file in os.listdir(args.in_files):
        df = pd.read_csv(args.in_files + '/' + file)
        df_all = pd.concat([df_all, df])
    df_all.to_csv(args.out + '/datasets.csv', index=False)

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
    print(perc_trp_residues)
    df_residues = pd.DataFrame({'total trp residues': [trp_res], '% trp residues': [perc_trp_residues]})
    df_residues.to_csv(args.in_files + '/covered_residues.csv', index=False)



def piecharts():
    df = pd.read_csv(args.out + '/datasets.csv')
    total_trp = df['total trp'].sum()
    total_prt = df['total prt'].sum()
    perc = round((total_trp * 100) / total_prt)
    diff_perc = 100 - perc
    print(perc)
    print(diff_perc)

    trp_curated = 3480
    total_curated = trp_curated + 1821
    perc_cur = round((trp_curated * 100) / total_curated)
    diff_perc = 100 - perc_cur

    fig, ax = plt.subplots()

    group_names = ['TRP organisms', 'not TRP organisms']
    subgroup_names = ['TRP curated', 'not TRP curated']
    size = 0.3
    vals = np.array([[perc, perc], [diff_perc, diff_perc]]) # percentuali
    vals2 = np.array([[perc_cur, perc_cur], [diff_perc, diff_perc]])

    cmap = plt.colormaps["tab20c"]
    outer_colors = cmap(np.arange(4) * 4)

    ax.pie(vals.sum(axis=1), radius=1, colors=outer_colors,labels=group_names,
           wedgeprops=dict(width=size, edgecolor='w'), labeldistance=0.5)


    ax.pie(vals2.sum(axis=1), radius=1 - size, colors=outer_colors, labels=subgroup_names,
           wedgeprops=dict(width=size, edgecolor='w'),  labeldistance=0.5)

    ax.set(aspect="equal", title='Pie plot with `ax.pie`')
    plt.savefig(args.out + "/piechart.png")
def boxplots():

    df_reg = pd.read_csv(args.out + '/datasets.csv')
    fig = px.box(df_reg, y="avg_region_len")
    fig.write_image(args.out + "/avg_region_len.png")

    df_uni = pd.read_csv(args.out + '/datasets.csv')
    fig = px.box(df_uni, y="avg_unit_lenght")
    fig.write_image(args.out + "/avg_unit_lenght.png")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate plots')
    parser.add_argument('--in_files', '-i', type=str, help='Path to input files')
    parser.add_argument('--out', '-o', type=str, help='Path to output')

    # args.in_files
    # data/af_4/stats

    # args.out
    # data/af_4/plots


    # Parse arguments
    args = parser.parse_args()

    # plot_table_dataset()
    # covered_residues()
    # piecharts()
    # boxplots()