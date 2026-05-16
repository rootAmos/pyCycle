"""Plot the reference mission profile from sizing/mission_ref.json."""

from argparse import ArgumentParser
import json
from pathlib import Path

import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MISSION_REF = REPO_ROOT / "sizing" / "mission_ref.json"


def load_mission(path):
    with Path(path).open("r", encoding="utf-8") as stream:
        return json.load(stream)


def build_profile(segments):
    distance = [0.0]
    altitude = []
    mach = []
    labels = []
    segment_edges = []

    if not segments:
        return distance, altitude, mach, labels, segment_edges

    altitude.append(float(segments[0]["start_altitude_ft"]))
    mach.append(float(segments[0]["start_mach"]))

    cumulative_distance = 0.0
    for segment in segments:
        segment_distance = float(segment.get("distance_nmi", 0.0))
        cumulative_distance += segment_distance

        distance.append(cumulative_distance)
        altitude.append(float(segment["end_altitude_ft"]))
        mach.append(float(segment["end_mach"]))
        labels.append(segment["name"])
        segment_edges.append(cumulative_distance)

    return distance, altitude, mach, labels, segment_edges


def plot_mission(mission, output_path=None, show=False):
    distance, altitude, mach, labels, segment_edges = build_profile(mission["segments"])

    fig, altitude_axis = plt.subplots(figsize=(11, 6))
    mach_axis = altitude_axis.twinx()

    altitude_axis.plot(distance, altitude, marker="o", color="tab:blue", label="Altitude")
    mach_axis.plot(distance, mach, marker="s", color="tab:red", label="Mach")

    for edge, label in zip(segment_edges, labels):
        altitude_axis.axvline(edge, color="0.85", linewidth=0.8, zorder=0)
        altitude_axis.text(
            edge,
            0.02,
            label.replace("_", " "),
            rotation=90,
            va="bottom",
            ha="right",
            fontsize=8,
            color="0.35",
            transform=altitude_axis.get_xaxis_transform(),
        )

    units = mission.get("units", {})
    altitude_axis.set_xlabel(f"Distance ({units.get('distance', 'nmi')})")
    altitude_axis.set_ylabel(f"Altitude ({units.get('altitude', 'ft')})", color="tab:blue")
    mach_axis.set_ylabel(f"Speed ({units.get('speed', 'Mach')})", color="tab:red")
    altitude_axis.tick_params(axis="y", labelcolor="tab:blue")
    mach_axis.tick_params(axis="y", labelcolor="tab:red")
    altitude_axis.grid(True, color="0.9")

    title = mission.get("mission_name", "Mission reference").replace("_", " ")
    altitude_axis.set_title(title)

    lines = altitude_axis.get_lines() + mach_axis.get_lines()
    altitude_axis.legend(lines, [line.get_label() for line in lines], loc="upper left")
    fig.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=200)
        print(f"Saved {output_path}")

    if show:
        plt.show()

    return fig


def parse_args():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "mission_ref",
        nargs="?",
        default=DEFAULT_MISSION_REF,
        type=Path,
        help="Path to the mission reference JSON file.",
    )
    parser.add_argument(
        "--output",
        "-o",
        default=REPO_ROOT / "sizing" / "mission_ref_profile.png",
        type=Path,
        help="Path for the output plot image. Use --no-save to skip writing a file.",
    )
    parser.add_argument("--no-save", action="store_true", help="Do not save the plot image.")
    parser.add_argument("--show", action="store_true", help="Display the plot window.")
    return parser.parse_args()


def main():
    args = parse_args()
    mission = load_mission(args.mission_ref)
    output_path = None if args.no_save else args.output
    plot_mission(mission, output_path=output_path, show=args.show)


if __name__ == "__main__":
    main()
