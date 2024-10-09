#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
from scipy.optimize import curve_fit

def deg_color(color, deg):
    import matplotlib.colors as mc
    import colorsys
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    #return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])
    return colorsys.hls_to_rgb(((c[0]*360)-deg)/360, c[1], c[2])

def load_data(filename):
    with open(filename, 'r') as f:
        return  [float(x) for x in f.read().split()]

def get_max_x(input_files):
    max_x = 0
    for input_file in input_files:
        y = load_data(input_file)
        max_x = max(len(y),max_x)
    return max_x

def exp_dec(x, a, b, c):
    return a * np.exp(-b * x) + c

def sort_by_length(input_files):
    files = []
    for input_file in input_files:
        y = load_data(input_file)
        files.append((input_file, len(y)))
    files = sorted(files, key=lambda x: -x[1])
    return [f for f,_ in files]

def plot_scatter_from_file(input_files, output_file):
    colors = sns.color_palette("tab10", len(input_files))

    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1)
    _, ax = plt.subplots(figsize=(6, 3))

    max_x = get_max_x(input_files)

    for idx, input_file in enumerate(input_files):
        y = load_data(input_file)
        x = [i for i in range(1,len(y)+1)]

        name = os.path.basename(input_file)
        if os.path.dirname(input_file) != "":
            name = os.path.dirname(input_file)+"/"+name

        sns.scatterplot(x=x, y=y, s=20,
                color=colors[idx], linewidth=1.2,
                edgecolor=colors[idx], alpha=0.7, label=name)
        sns.lineplot(x=x, y=y, color=colors[idx], alpha=0.7, linewidth=1.3)

        x_fit = np.linspace(1, max_x, num=1000, dtype='float')

        p0 = (max(y), 1, min(y))
        popt, _ = curve_fit(exp_dec, range(1, len(y) + 1), y, p0=p0)
        print(rf'{name} core: ({popt[0]:.2e}) e^(-{popt[1]:.2f}x) + {popt[2]:.2f}')
        print(rf"{name} percentage of core in a genome:")
        print("    (predicted core)/(average items in a single genome)")
        percentage_core = popt[2]/y[0]
        print(rf"    {popt[2]:.2f}/{y[0]:.2f} = {percentage_core}")
        sns.lineplot(x=x_fit, y=exp_dec(x_fit, *popt),
                     color=deg_color(colors[idx], 10), alpha=0.8, linewidth=1.8,
                     label=rf'({popt[0]:.2e})$\cdot e^{{ -{popt[1]:.2f}x }} + {popt[2]:.2f}$',
                     linestyle='--')

    plt.xlabel('No. of Genomes')
    plt.ylabel('$k$-mers')

    # Save the figure
    plt.tight_layout()
    plt.savefig(output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot pangenome growth from histogram(s).")
    parser.add_argument("input_files", nargs='+', help="Input histogram file(s).")
    parser.add_argument("output_file", help="Output figure file in pdf.")

    args = parser.parse_args()

    # Check if the output file has a .pdf extension; if not, add it
    root, ext = os.path.splitext(args.output_file)
    if ext.lower() != '.pdf':
        args.output_file = args.output_file + '.pdf'

    plot_scatter_from_file(sort_by_length(args.input_files), args.output_file)
