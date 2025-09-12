#!/usr/bin/env python3
import argparse
import gzip
import io
import os
import sys

def open_maybe_gzip(path):
    if path in ("-", "", None):
        return sys.stdin
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")

def detect_style(attrs_field):
    # crude but effective: GTF uses key "value"; GFF3 uses key=value;
    if "=" in attrs_field and '"' not in attrs_field:
        return "gff3"
    return "gtf"

def parse_attrs_gtf(field):
    # GTF attributes: key "value"; key2 "value2"; ...
    out = {}
    # Split by ; then parse key "value"
    for part in field.strip().strip(";").split(";"):
        if not part.strip():
            continue
        # parts may be: key "value" or key "value";
        try:
            key, val = part.strip().split(" ", 1)
        except ValueError:
            # malformed; skip
            continue
        val = val.strip().strip('"')
        out[key] = val
    return out

def parse_attrs_gff3(field):
    # GFF3 attributes: key=value;key2=value2
    out = {}
    for part in field.strip().split(";"):
        if not part:
            continue
        if "=" not in part:
            # Some files have bare flags; ignore
            continue
        key, val = part.split("=", 1)
        # GFF3 values can be URL-encoded; do a minimal decode for common chars
        val = val.replace("%3B", ";").replace("%2C", ",").replace("%3D", "=").replace("%25", "%")
        out[key] = val
    return out

def main():
    ap = argparse.ArgumentParser(
        description="Create tx2gene.tsv (TXNAME, GENEID) from a GTF or GFF3 file."
    )
    ap.add_argument("input", help="Input GTF/GFF3 file (use - for stdin)")
    ap.add_argument("-o", "--output", default="tx2gene.tsv", help="Output TSV file (default: tx2gene.tsv)")
    ap.add_argument("--tx-field", default="transcript_id",
                    help="Attribute holding transcript ID (default: transcript_id). "
                         "For GFF3 this is often ID (for mRNA) or transcript_id if provided.")
    ap.add_argument("--gene-field", default="gene_id",
                    help="Attribute holding gene ID (default: gene_id). "
                         "For GFF3 this is often Parent (for mRNA) or gene_id if provided.")
    ap.add_argument("--format", choices=["auto", "gtf", "gff3"], default="auto",
                    help="Force parse style (auto tries to infer per-file; default: auto)")
    args = ap.parse_args()

    # We’ll gather pairs; last one wins if duplicates; we’ll deduplicate at the end.
    tx2gene = {}

    with open_maybe_gzip(args.input) as fh:
        inferred_style = None
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue  # not a valid GTF/GFF row
            attrs_field = cols[8]

            style = args.format
            if style == "auto":
                if inferred_style is None:
                    inferred_style = detect_style(attrs_field)
                style = inferred_style

            if style == "gtf":
                attrs = parse_attrs_gtf(attrs_field)
            else:  # gff3
                attrs = parse_attrs_gff3(attrs_field)

            tx = attrs.get(args.tx_field)
            gene = attrs.get(args.gene_field)

            # Helpful fallbacks for common GFF3 conventions if user left defaults:
            if tx is None and style == "gff3":
                # Transcript often on mRNA rows: ID=<tx>; gene might be Parent=<gene>
                if args.tx_field == "transcript_id":
                    tx = attrs.get("ID")
            if gene is None and style == "gff3":
                if args.gene_field == "gene_id":
                    # Most GFF3s have Parent linking mRNA->gene
                    gene = attrs.get("Parent")

            if tx is None or gene is None:
                continue  # skip rows that lack either field

            # Some GFF3 Parent can be a comma-separated list; take first (typical for mRNA->gene)
            if "," in gene:
                gene = gene.split(",", 1)[0]

            tx2gene[tx] = gene

    # Write output
    out_path = args.output
    with (sys.stdout if out_path in ("-", "", None) else open(out_path, "w", encoding="utf-8")) as out:
        out.write("TXNAME\tGENEID\n")
        for tx, gene in tx2gene.items():
            out.write(f"{tx}\t{gene}\n")

    if out_path not in ("-", "", None):
        sys.stderr.write(f"Wrote {len(tx2gene)} pairs to {out_path}\n")

if __name__ == "__main__":
    main()
