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

def power_law_growth(x, a, b):
    return a * x **b

def power_law_new(x, y):
    log_x = np.log(x)
    log_y = np.log(y)
    b, log_a = np.polyfit(log_x, log_y, 1)
    a = np.exp(log_a)
    return [a, b]

def load_data(filename):
    with open(filename, 'r') as f:
        return  [float(x) for x in f.read().split()]

def get_max_x(input_files):
    max_x = 0
    for input_file in input_files:
        y = load_data(input_file)
        max_x = max(len(y),max_x)
    return max_x

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
    sns.set_context("paper", font_scale=2)
    _, (ax1,ax2) = plt.subplots(1,2, figsize=(16, 7))

    max_x = get_max_x(input_files)

    for idx, input_file in enumerate(input_files):
        y = load_data(input_file)
        x = [i for i in range(1,len(y)+1)]

        name = os.path.basename(input_file)
        if os.path.dirname(input_file) != "":
            name = os.path.dirname(input_file)+"/"+name

        sns.scatterplot(x=x, y=y, s=110,
                color=colors[idx], linewidth=3,
                edgecolor=colors[idx], alpha=0.7, label=name, ax=ax1)
        x_fit = np.linspace(1, max_x, num=1000, dtype='float')

        p0 = (max(y), 1)
        popt, _ = curve_fit(power_law_growth, range(1, len(y) + 1), y, p0=p0)
        print(rf'{name} growth: ({popt[0]:.2e}) x^{popt[1]:.2f}')
        sns.lineplot(x=x_fit, y=power_law_growth(x_fit, *popt),
                     color=deg_color(colors[idx], 10), alpha=0.9, linewidth=3.3,
                     label=rf'({popt[0]:.2e})$\cdot x^{{ {popt[1]:.2f} }}$', ax=ax1,
                     linestyle='--')

    for idx, input_file in enumerate(input_files):
        y = load_data(input_file)
        y_new = y.copy()
        for i in range(1,len(y)):
            y_new[i] = y[i]-y[i-1]
        y_new = y_new[1:]
        x = [i for i in range(2, len(y)+1)]
        y = y_new

        name = os.path.basename(input_file)
        if os.path.dirname(input_file) != "":
            name = os.path.dirname(input_file)+"/"+name

        sns.scatterplot(x=x, y=y, s=110,
                color=colors[idx], linewidth=3,
                edgecolor=colors[idx], alpha=0.7, label=name, ax=ax2)

        x_fit = np.linspace(2, max_x, num=1000, dtype='float')
        [a,b] = power_law_new(x,y)
        fit_y = a*x_fit**b
        print(rf'{name} average new: ({a:.2e}) x^{b:.2f}')
        sns.lineplot(x=x_fit, y=fit_y, color=deg_color(colors[idx], 10),
                     alpha=0.9, linewidth=3.3, 
                     label=rf'({a:.2e})$\cdot x^{{ {b:.2f} }}$',ax=ax2, linestyle='--')


    ax1.set_xlabel('No. of Genomes')
    ax1.set_ylabel('$k$-mers')
    #ax1.legend()
    ax2.set_xlabel('No. of Genomes')
    ax2.set_ylabel('$k$-mers')
    #ax2.legend()

    # Save the figure
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
