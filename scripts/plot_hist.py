#!/usr/bin/env python3

"""Plot one or more pangrowth frequency histograms."""

import argparse
import math
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


PERCENTAGE_BINS = list(range(5, 101, 5))


def load_counts(filename, norm_y):
    with open(filename, "r", encoding="utf-8") as stream:
        try:
            counts = [int(line.strip()) for line in stream if line.strip()]
        except ValueError as error:
            raise ValueError(f"{filename}: histogram values must be integers") from error

    if not counts:
        raise ValueError(f"{filename}: the histogram is empty")
    if any(value < 0 for value in counts):
        raise ValueError(f"{filename}: histogram values must be non-negative")

    values = [float(value) for value in counts]
    if norm_y in ("multiplicity", "both"):
        values = [value * frequency for frequency, value in enumerate(values, 1)]
    if norm_y in ("percentage", "both"):
        total = sum(values)
        if total == 0:
            raise ValueError(f"{filename}: cannot normalize an all-zero histogram")
        values = [value / total for value in values]
    return values


def make_dataframe(filename, label, norm_x, norm_y):
    values = load_counts(filename, norm_y)
    genome_count = len(values)

    if not norm_x:
        return pd.DataFrame(
            {
                "x": range(1, genome_count + 1),
                "y": values,
                "file": label,
            }
        )

    # Assign every frequency to the upper edge of a five-percentage-point bin.
    # Only x is normalized; the k-mer counts retain their original scale.
    binned = {percentage: 0.0 for percentage in PERCENTAGE_BINS}
    for frequency, value in enumerate(values, 1):
        percentage = min(100, 5 * math.ceil(20 * frequency / genome_count))
        binned[percentage] += value
    return pd.DataFrame(
        {
            "x": PERCENTAGE_BINS,
            "y": [binned[percentage] for percentage in PERCENTAGE_BINS],
            "file": label,
        }
    )


def readable_tick_indices(number_of_ticks, maximum=12):
    """Return evenly spaced tick indices while preserving both endpoints."""
    if number_of_ticks <= maximum:
        return list(range(number_of_ticks))
    step = math.ceil((number_of_ticks - 1) / (maximum - 1))
    indices = list(range(0, number_of_ticks - 1, step))
    if indices[-1] != number_of_ticks - 1:
        indices.append(number_of_ticks - 1)
    return indices


def plot_histograms(input_files, output_file, norm_x=False, norm_y="none", labels=None):
    labels = labels or [os.path.basename(filename) for filename in input_files]
    if len(labels) != len(input_files):
        raise ValueError(
            f"received {len(labels)} labels for {len(input_files)} input files"
        )

    dataframes = [
        make_dataframe(filename, label, norm_x, norm_y)
        for filename, label in zip(input_files, labels)
    ]
    data = pd.concat(dataframes, ignore_index=True)
    if norm_x:
        order = PERCENTAGE_BINS
    else:
        maximum_frequency = max(int(dataframe["x"].max()) for dataframe in dataframes)
        order = list(range(1, maximum_frequency + 1))

    with sns.axes_style("whitegrid"), sns.plotting_context("paper", font_scale=1.0):
        figure, axis = plt.subplots(figsize=(7, 3))
        sns.barplot(
            x="x",
            y="y",
            hue="file",
            data=data,
            order=order,
            edgecolor="none",
            ci=None,
            ax=axis,
        )

        if norm_x:
            tick_indices = list(range(len(order)))
            tick_labels = [f"{percentage}%" for percentage in order]
            rotation = 45
            font_size = 7
        else:
            tick_indices = readable_tick_indices(len(order))
            tick_labels = [str(order[index]) for index in tick_indices]
            rotation = 0
            font_size = 8
        axis.set_xticks(tick_indices)
        axis.set_xticklabels(
            tick_labels, rotation=rotation, ha="right" if rotation else "center"
        )
        axis.tick_params(
            axis="x", which="both", top=False, labelsize=font_size, width=0.45
        )
        axis.tick_params(axis="y", which="both", right=False, width=0.45)
        axis.grid(which="major", axis="y", color="#d9d9d9", linewidth=0.42)
        axis.grid(which="major", axis="x", visible=False)
        axis.grid(which="minor", visible=False)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.set_xlabel("Percentage of genomes" if norm_x else "No. of genomes")
        axis.set_ylabel(
            "Normalized $k$-mers" if norm_y != "none" else "$k$-mers"
        )
        if norm_y in ("percentage", "both"):
            axis.set_ylim(0, 1.0)
        axis.legend(loc="upper center", frameon=False)

        figure.tight_layout()
        figure.savefig(output_file, bbox_inches="tight")
        plt.close(figure)


def main():
    parser = argparse.ArgumentParser(
        description="Plot one or more pangrowth frequency histograms."
    )
    parser.add_argument("input_files", nargs="+", help="Input histogram file(s).")
    parser.add_argument("output_file", help="Output figure (PDF by default).")
    parser.add_argument(
        "--norm_x",
        action="store_true",
        help="Group genome frequencies into five-percentage-point bins.",
    )
    parser.add_argument(
        "--norm_y",
        choices=("multiplicity", "percentage", "both", "none"),
        default="none",
        help="Optional y-axis normalization (not used by the CloWM workflow).",
    )
    parser.add_argument(
        "--label",
        action="append",
        dest="labels",
        help="Legend label for an input file; repeat once per input file.",
    )
    args = parser.parse_args()

    root, extension = os.path.splitext(args.output_file)
    if not extension:
        args.output_file = root + ".pdf"

    try:
        plot_histograms(
            args.input_files,
            args.output_file,
            norm_x=args.norm_x,
            norm_y=args.norm_y,
            labels=args.labels,
        )
    except (OSError, ValueError) as error:
        parser.exit(1, f"plot_hist.py: error: {error}\n")


if __name__ == "__main__":
    main()
