#!/usr/bin/env python3
"""Concatenate aligned FASTA files and write an IQ-TREE partition file."""

from __future__ import annotations

import argparse
from pathlib import Path


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []

    def finish() -> None:
        nonlocal name, chunks
        if name is None:
            return
        if name in records:
            raise ValueError(f"Duplicate FASTA header {name!r} in {path}")
        sequence = "".join(chunks).replace(" ", "").replace("\t", "")
        if not sequence:
            raise ValueError(f"Empty sequence for {name!r} in {path}")
        records[name] = sequence
        name = None
        chunks = []

    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                name = line[1:].split()[0]
            else:
                if name is None:
                    raise ValueError(f"Sequence data before FASTA header in {path}")
                chunks.append(line)
    finish()
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    lengths = {len(sequence) for sequence in records.values()}
    if len(lengths) != 1:
        raise ValueError(f"Unequal sequence lengths in {path}: {sorted(lengths)}")
    return records


def write_fasta(path: Path, sequences: dict[str, str]) -> None:
    with path.open("w") as handle:
        for name, sequence in sequences.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment-dir", type=Path, required=True)
    parser.add_argument("--genes", nargs="+", required=True)
    parser.add_argument("--suffix", default=".aln.fa")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--partitions", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()

    alignments: list[tuple[str, dict[str, str], int]] = []
    taxa: set[str] = set()
    for gene in args.genes:
        path = args.alignment_dir / f"{gene}{args.suffix}"
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Alignment not found or empty: {path}")
        records = read_fasta(path)
        length = len(next(iter(records.values())))
        alignments.append((gene, records, length))
        taxa.update(records)

    ordered_taxa = sorted(taxa)
    concatenated = {taxon: "" for taxon in ordered_taxa}
    start = 1
    manifest_rows = ["gene\tlength\tstart\tend\tn_taxa\tn_missing"]
    partition_rows: list[str] = []
    for gene, records, length in alignments:
        end = start + length - 1
        partition_rows.append(f"DNA, {gene} = {start}-{end}")
        missing = 0
        for taxon in ordered_taxa:
            sequence = records.get(taxon)
            if sequence is None:
                sequence = "N" * length
                missing += 1
            concatenated[taxon] += sequence
        manifest_rows.append(
            f"{gene}\t{length}\t{start}\t{end}\t{len(records)}\t{missing}"
        )
        start = end + 1

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.partitions.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(args.output, concatenated)
    args.partitions.write_text("\n".join(partition_rows) + "\n")
    args.manifest.write_text("\n".join(manifest_rows) + "\n")


if __name__ == "__main__":
    main()
