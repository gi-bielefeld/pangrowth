#!/usr/bin/env python3

import seaborn as sns
import argparse
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
import pandas as pd
import sys
import os

def get_name_file(filename):
    name = os.path.basename(filename)
    if os.path.dirname(filename) != "":
        name = os.path.dirname(filename)+"/"+name
    return name

def load_data(filename, norm_x, norm_y):
    with open(filename, 'r') as f:
        data = [int(x.strip()) for x in f.readlines()]
        n = len(data)
        y = data
        if norm_y == "multiplicity":
            y = [y[i]*(i+1) for i in range(len(y))]
        if norm_y == "percentage":
            y = [y[i]/sum(data) for i in range(len(y))]
        if norm_y == "both":
            y = [y[i]*(i+1) for i in range(len(y))]
            y = [y[i]/sum(y) for i in range(len(y))]

        x = [(i+1) for i in range(n)]

        name = get_name_file(filename)

        df = pd.DataFrame({"x":x, "y":y, "file":name})

        if norm_x:
            ranges = np.arange(0, 1.05, 0.05)
            x = [(i+1)/n for i in range(n)]

            y_norm = df["y"].groupby(pd.cut(x, bins=ranges)).sum()
            y_norm /= y_norm.sum()

            df_norm = {}
            df_norm["x"] = y_norm.index
            df_norm["y"] = y_norm.values
            df_norm["file"] = [name]*len(df_norm["y"])
            df = pd.DataFrame(df_norm)

        return df

def main():
    parser = argparse.ArgumentParser(description="Plot a multiple histograms from histogram files.")
    parser.add_argument("input_files", nargs='+', help="Input histogram file(s).")
    parser.add_argument("output_file", help="Output figure file in pdf.")
    parser.add_argument("--norm_x", action="store_true", default=False, help="Normalize the x-axis between (0, 1.0] to compare histograms with different number of genomes")
    parser.add_argument('--norm_y', choices=['multiplicity', 'percentage', 'both', 'none'], 
                        default="none", 
                        help='Normalize the y-axis values: '
                         '"multiplicity" multiply each bar by the number of genomes it appears (unique items are multiplied by one, items appearing in two genomes are multiplied by two, and so on..), '
                         '"percentage" normalizes by converting counts to a percentage of the total, '
                         '"both" applies both multiplicity and percentage normalizations in sequence, '
                         '"none" applies no normalization.')
    parser.add_argument("--max_tot_bars", action="store_true", default=190, help="Maximum number of bars shown. It considers the number of files passed.")
    args = parser.parse_args()

    sns.set_context("paper")
    sns.set_style("whitegrid")
    _, ax = plt.subplots(figsize=(7, 3))

    #dd = {"x":load_data(args.input_files[0],args.norm).index}
    dataframes = []
    x_max = 0
    n_files = len(args.input_files)
    for filename in args.input_files:
        dataframe = load_data(filename, args.norm_x, args.norm_y)
        dataframes.append(dataframe)
        if not args.norm_x:
            x_max = max(x_max, max(dataframe["x"]))

    max_tot_bars = args.max_tot_bars
    max_bars = np.ceil(max_tot_bars/n_files)

    if not args.norm_x and x_max*n_files > max_tot_bars:
        print("WARNING: The number of bars exceeds the maximum allowed ('max_tot_bars'). Bars are now grouped into intervals, representing the sum of values.")
        bin_width = np.ceil(x_max / max_bars)
        ranges = np.arange(0, x_max + bin_width, bin_width)

        dataframes_norm = []
        for i in range(len(dataframes)):
            dataframe = dataframes[i]
            name = get_name_file(args.input_files[i])
            y_norm = dataframe["y"].groupby(pd.cut(dataframe["x"], bins=ranges)).sum()
            y_norm /= y_norm.sum()
            df_norm = {}
            df_norm["x"] = y_norm.index
            df_norm["y"] = y_norm.values
            df_norm["file"] = [name]*len(df_norm["y"])
            dataframes_norm.append(pd.DataFrame(df_norm))

        dataframes = dataframes_norm

    df = pd.concat(dataframes)

    sns.barplot(x='x', y='y', hue='file', data=df, edgecolor='none')
    sns.set_style("whitegrid")

    plt.gca().spines["top"].set_visible(True)
    plt.gca().spines["right"].set_visible(True)
    plt.gca().spines["bottom"].set_visible(True)
    plt.gca().spines["left"].set_visible(True)
    plt.gca().spines["left"].set_visible(True)
    line_width = 0.28
    plt.gca().spines["left"].set_linewidth(line_width)
    plt.gca().spines["right"].set_linewidth(line_width)
    plt.gca().spines["bottom"].set_linewidth(line_width)
    plt.gca().spines["top"].set_linewidth(line_width)
    plt.grid(True, linewidth=line_width)
    # Axis labels
    if not args.norm_x and x_max*n_files > max_tot_bars:
        plt.xticks(rotation=315,fontsize=6,ha="left", rotation_mode="anchor")
    elif not args.norm_x:
        plt.xticks(fontsize=max(7-np.ceil(n_files/26),5),ha="left", rotation_mode="anchor")
    else:
        plt.xticks(rotation=315,fontsize=8,ha="left", rotation_mode="anchor")
    # Set the x and y labels and the title
    plt.xlabel('Percentage of Genomes') 
    if args.norm_y != "none":
        plt.ylabel('Normalized $k$-mers')
    else:
        plt.ylabel('$k$-mers')

    if args.norm_y == "percentage" or args.norm_y == "both":
        plt.ylim(0,1.0)

    # Legend
    plt.tight_layout()
    font = FontProperties()
    font.set_style('italic')
    plt.legend(loc="upper center",prop=font)

    root, ext = os.path.splitext(args.output_file)
    if ext.lower() != '.pdf':
        args.output_file = args.output_file + '.pdf'
    plt.savefig(args.output_file, bbox_inches='tight')

if __name__ == "__main__":
    main()
