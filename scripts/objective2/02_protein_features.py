from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO
from collections import Counter
import argparse
import json
import math
import os
import time

import pandas as pd


# ============================================================
# PROJECT CONFIGURATION
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DEFAULT_INPUT = (
    PROJECT_ROOT
    / "data"
    / "raw"
    / "alternaria"
    / "proteome"
    / "alternaria_alternata_GCF_001642055.1_protein.faa"
)

DEFAULT_OUTPUT = (
    PROJECT_ROOT
    / "data"
    / "processed"
    / "annotation"
    / "alternaria"
)

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

HYDROPHOBIC = set("AVILMFWY")
POSITIVELY_CHARGED = set("KRH")
NEGATIVELY_CHARGED = set("DE")
AROMATIC = set("FWY")


# ============================================================
# FEATURE EXTRACTION
# ============================================================

def extract_features(item):
    protein_id, sequence = item

    sequence = str(sequence).upper()

    # Keep only the 20 standard amino acids for composition.
    # Record the original sequence length separately.
    original_length = len(sequence)
    clean_sequence = "".join(
        aa for aa in sequence if aa in AMINO_ACIDS
    )

    length = len(clean_sequence)
    counts = Counter(clean_sequence)

    row = {
        "protein_id": protein_id,
        "sequence_length": original_length,
        "standard_aa_length": length,
        "hydrophobic_fraction": (
            sum(counts[aa] for aa in HYDROPHOBIC) / length
            if length else 0.0
        ),
        "positive_fraction": (
            sum(counts[aa] for aa in POSITIVELY_CHARGED) / length
            if length else 0.0
        ),
        "negative_fraction": (
            sum(counts[aa] for aa in NEGATIVELY_CHARGED) / length
            if length else 0.0
        ),
        "aromatic_fraction": (
            sum(counts[aa] for aa in AROMATIC) / length
            if length else 0.0
        ),
        "sequence_entropy": (
            -sum(
                (count / length) * math.log2(count / length)
                for count in counts.values()
                if count
            )
            if length else 0.0
        ),
    }

    for aa in AMINO_ACIDS:
        row[f"frac_{aa}"] = (
            counts[aa] / length if length else 0.0
        )

    return row


# ============================================================
# INPUT DISCOVERY
# ============================================================

def find_protein_fasta(explicit_path=None):
    if explicit_path:
        path = Path(explicit_path).expanduser()

        if not path.is_absolute():
            path = PROJECT_ROOT / path

        if not path.is_file():
            raise FileNotFoundError(
                f"Specified protein FASTA not found: {path}"
            )

        return path.resolve()

    # Prefer the known proteome file.
    if DEFAULT_INPUT.is_file():
        return DEFAULT_INPUT

    # Fall back to another .faa in the Alternaria directory.
    search_root = PROJECT_ROOT / "data" / "raw" / "alternaria"

    candidates = sorted(
        search_root.rglob("*.faa"),
        key=lambda p: (
            "proteome" not in str(p).lower(),
            len(str(p))
        )
    )

    if not candidates:
        raise FileNotFoundError(
            "No .faa file found under data/raw/alternaria. "
            "Check the project directory and protein filename."
        )

    if len(candidates) > 1:
        print("Multiple protein FASTA files found:")
        for candidate in candidates:
            print(f"  {candidate}")
        print(f"Using: {candidates[0]}")

    return candidates[0].resolve()


# ============================================================
# MAIN PIPELINE
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract protein sequence features and benchmark "
            "serial versus parallel CPU processing."
        )
    )

    parser.add_argument(
        "--input",
        default=None,
        help="Optional path to the fungal protein FASTA."
    )

    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT),
        help="Directory for feature and benchmark outputs."
    )

    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of parallel worker processes."
    )

    parser.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="Number of timing repetitions."
    )

    args = parser.parse_args()

    if args.workers < 1 or args.repeats < 1:
        raise ValueError(
            "--workers and --repeats must both be at least 1."
        )

    input_path = find_protein_fasta(args.input)
    output_dir = Path(args.output_dir).expanduser()

    if not output_dir.is_absolute():
        output_dir = PROJECT_ROOT / output_dir

    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("ALTERNARIA PROTEIN FEATURE EXTRACTION")
    print("=" * 60)
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Input FASTA:  {input_path}")
    print(f"Output dir:   {output_dir}")
    print(f"CPU cores:    {os.cpu_count()}")
    print(f"Workers:      {args.workers}")
    print(f"Repeats:      {args.repeats}")

    # Read actual protein sequences.
    records = [
        (record.id, str(record.seq))
        for record in SeqIO.parse(input_path, "fasta")
    ]

    if not records:
        raise ValueError(
            f"No protein sequences were read from {input_path}"
        )

    print(f"\nProtein records loaded: {len(records)}")

    # Serial benchmark.
    serial_times = []

    for run in range(args.repeats):
        start = time.perf_counter()

        serial_rows = [
            extract_features(record)
            for record in records
        ]

        elapsed = time.perf_counter() - start
        serial_times.append(elapsed)

        print(
            f"Serial run {run + 1}: {elapsed:.4f} seconds"
        )

    # Parallel benchmark.
    parallel_times = []
    parallel_rows = None

    for run in range(args.repeats):
        start = time.perf_counter()

        with ProcessPoolExecutor(
            max_workers=args.workers
        ) as executor:
            parallel_rows = list(
                executor.map(
                    extract_features,
                    records,
                    chunksize=100
                )
            )

        elapsed = time.perf_counter() - start
        parallel_times.append(elapsed)

        print(
            f"Parallel run {run + 1}: {elapsed:.4f} seconds"
        )

    # Verify both methods generated matching IDs and row counts.
    serial_ids = [row["protein_id"] for row in serial_rows]
    parallel_ids = [row["protein_id"] for row in parallel_rows]

    if serial_ids != parallel_ids:
        raise RuntimeError(
            "Serial and parallel outputs do not match."
        )

    # Save features.
    feature_path = output_dir / "protein_features.csv"

    pd.DataFrame(parallel_rows).to_csv(
        feature_path,
        index=False
    )

    serial_mean = sum(serial_times) / len(serial_times)
    parallel_mean = sum(parallel_times) / len(parallel_times)

    speedup = (
        serial_mean / parallel_mean
        if parallel_mean > 0 else None
    )

    efficiency = (
        speedup / args.workers
        if speedup is not None else None
    )

    benchmark = {
        "input_file": str(input_path.relative_to(PROJECT_ROOT)),
        "protein_count": len(records),
        "cpu_cores_reported": os.cpu_count(),
        "workers": args.workers,
        "repeats": args.repeats,
        "serial_times_seconds": serial_times,
        "parallel_times_seconds": parallel_times,
        "serial_mean_seconds": serial_mean,
        "parallel_mean_seconds": parallel_mean,
        "speedup": speedup,
        "parallel_efficiency": efficiency,
        "environment": "Local CPU multiprocessing",
        "note": (
            "This is a local benchmark, not a distributed HPC "
            "cluster benchmark. Results depend on system load "
            "and multiprocessing overhead."
        )
    }

    benchmark_path = output_dir / "hpc_benchmark.json"

    benchmark_path.write_text(
        json.dumps(benchmark, indent=4),
        encoding="utf-8"
    )

    # Save a compact text summary.
    summary_path = output_dir / "hpc_benchmark_summary.txt"

    summary = [
        "ALTERNARIA PROTEIN PROCESSING BENCHMARK",
        "=" * 42,
        f"Input: {benchmark['input_file']}",
        f"Protein count: {len(records)}",
        f"Workers: {args.workers}",
        f"Repeats: {args.repeats}",
        f"Mean serial time: {serial_mean:.4f} seconds",
        f"Mean parallel time: {parallel_mean:.4f} seconds",
        f"Speedup: {speedup:.4f}x" if speedup is not None
        else "Speedup: unavailable",
        f"Parallel efficiency: {efficiency:.4f}"
        if efficiency is not None
        else "Parallel efficiency: unavailable",
        "",
        "Note: local CPU benchmark; not a distributed HPC run."
    ]

    summary_path.write_text(
        "\n".join(summary),
        encoding="utf-8"
    )

    print("\n" + "=" * 60)
    print("COMPLETED")
    print("=" * 60)
    print(f"Feature table: {feature_path}")
    print(f"Benchmark JSON: {benchmark_path}")
    print(f"Benchmark summary: {summary_path}")
    print(f"Mean serial time: {serial_mean:.4f} seconds")
    print(f"Mean parallel time: {parallel_mean:.4f} seconds")

    if speedup is not None:
        print(f"Speedup: {speedup:.4f}x")
        print(f"Parallel efficiency: {efficiency:.4f}")


if __name__ == "__main__":
    main()