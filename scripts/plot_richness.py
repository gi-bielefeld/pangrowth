#!/usr/bin/env python3

"""Plot interpolated, observed, and extrapolated pangenome richness."""

import argparse
import os

import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D
import pandas as pd
import seaborn as sns


FIT_VALUES = {"int", "obs", "ext"}


def load_data(filename):
    """Load and validate a richness table produced by ``pangrowth hill``."""
    data = pd.read_csv(filename, sep="\t")
    required = {"fit", "m", "richness"}
    missing = sorted(required.difference(data.columns))
    if missing:
        raise ValueError("missing required column(s): " + ", ".join(missing))

    data = data.loc[:, ["fit", "m", "richness"]].copy()
    data["fit"] = data["fit"].astype(str).str.strip().str.lower()
    unknown = sorted(set(data["fit"]).difference(FIT_VALUES))
    if unknown:
        raise ValueError("unknown fit value(s): " + ", ".join(unknown))

    for column in ("m", "richness"):
        data[column] = pd.to_numeric(data[column], errors="raise")
    if data.empty:
        raise ValueError("the input table has no data rows")
    if (data["m"] <= 0).any():
        raise ValueError("column 'm' must contain positive sample sizes")
    if data["m"].duplicated().any():
        duplicates = sorted(data.loc[data["m"].duplicated(keep=False), "m"].unique())
        raise ValueError(
            "duplicate sample size(s) in column 'm': "
            + ", ".join(map(str, duplicates))
        )
    if (data["richness"] < 0).any():
        raise ValueError("richness must be non-negative")

    observed = data[data["fit"] == "obs"]
    if len(observed) != 1:
        raise ValueError(
            f"expected exactly one observed row ('obs'); found {len(observed)}"
        )
    return data.sort_values("m"), observed.iloc[0]


def segment(data, fit):
    """Return an interpolation or extrapolation segment joined to observation."""
    return data[data["fit"].isin((fit, "obs"))].sort_values("m")


def plot_richness(input_files, output_file, labels=None):
    labels = labels or [os.path.splitext(os.path.basename(path))[0] for path in input_files]
    if len(labels) != len(input_files):
        raise ValueError(
            f"received {len(labels)} labels for {len(input_files)} input files"
        )

    datasets = []
    for filename, label in sorted(
        zip(input_files, labels), key=lambda item: item[1].casefold()
    ):
        data, observed = load_data(filename)
        datasets.append((data, observed, label))

    colors = sns.color_palette("tab10", n_colors=len(datasets))
    richness_max = max(data["richness"].max() for data, _, _ in datasets)
    y_max = richness_max * 1.05 if richness_max > 0 else 1.0
    legend_columns = min(3, len(datasets))
    legend_rows = (len(datasets) + legend_columns - 1) // legend_columns

    with (
        sns.axes_style("whitegrid"),
        sns.plotting_context("paper", font_scale=1.15),
        plt.rc_context(
            {
                "font.family": "serif",
                "font.size": 9,
                "axes.titlesize": 10,
                "axes.labelsize": 10,
                "axes.linewidth": 0.55,
                "xtick.labelsize": 8,
                "ytick.labelsize": 8,
                "legend.fontsize": 8,
                "lines.solid_capstyle": "round",
                "lines.dash_capstyle": "butt",
            }
        ),
    ):
        figure_height = 3.25 + 0.22 * (legend_rows - 1)
        figure, axis = plt.subplots(figsize=(5.0, figure_height))

        for (data, observed, _), color in zip(datasets, colors):
            interpolation = segment(data, "int")
            extrapolation = segment(data, "ext")

            if len(interpolation) > 1:
                sns.lineplot(
                    x=interpolation["m"],
                    y=interpolation["richness"],
                    color=color,
                    linewidth=1.8,
                    ax=axis,
                )
            if len(extrapolation) > 1:
                sns.lineplot(
                    x=extrapolation["m"],
                    y=extrapolation["richness"],
                    color=color,
                    linewidth=1.8,
                    linestyle=(0, (2.2, 2.0)),
                    ax=axis,
                )
            sns.scatterplot(
                x=[observed["m"]],
                y=[observed["richness"]],
                s=30,
                color=color,
                edgecolor=color,
                linewidth=0.5,
                zorder=3,
                ax=axis,
            )

        axis.set_xlabel("No. of genomes")
        axis.set_ylabel("Richness")
        axis.tick_params(
            which="both", direction="in", top=False, right=False,
            width=0.45, length=3.5,
        )
        axis.minorticks_on()
        axis.grid(which="major", color="#d9d9d9", linewidth=0.42)
        axis.grid(which="minor", visible=False)
        axis.ticklabel_format(axis="y", style="sci", scilimits=(0, 0), useMathText=True)
        axis.margins(x=0.03)
        axis.set_ylim(0, y_max)

        legend_symbols = [
            (
                Line2D([], [], color=color, linewidth=1.8),
                Line2D(
                    [], [], color=color, marker="o", linestyle="none", markersize=4.8
                ),
                Line2D(
                    [], [], color=color, linewidth=1.8,
                    linestyle=(0, (2.2, 2.0)),
                ),
            )
            for color in colors
        ]
        figure.legend(
            legend_symbols,
            [label for _, _, label in datasets],
            handler_map={tuple: HandlerTuple(ndivide=None, pad=0.12)},
            prop={"style": "italic", "size": 8},
            loc="lower center",
            bbox_to_anchor=(0.5, -0.015),
            frameon=False,
            handlelength=4.8,
            ncol=legend_columns,
        )
        bottom = 0.16 + 0.06 * (legend_rows - 1)
        figure.tight_layout(rect=(0, bottom, 1, 1))
        figure.savefig(output_file, bbox_inches="tight", dpi=300)
        plt.close(figure)


def main():
    parser = argparse.ArgumentParser(
        description="Plot richness interpolation, observation, and extrapolation."
    )
    parser.add_argument(
        "input_files", nargs="+", help="TSV file(s) produced by 'pangrowth hill'."
    )
    parser.add_argument("output_file", help="Output figure (PDF by default).")
    parser.add_argument(
        "--label", action="append", dest="labels",
        help="Legend label for an input file; repeat once per input file.",
    )
    args = parser.parse_args()

    root, extension = os.path.splitext(args.output_file)
    if not extension:
        args.output_file = root + ".pdf"

    try:
        plot_richness(args.input_files, args.output_file, args.labels)
    except (OSError, ValueError, pd.errors.ParserError) as error:
        parser.exit(1, f"plot_richness.py: error: {error}\n")


if __name__ == "__main__":
    main()
