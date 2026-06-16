import re

# ─── CONFIG ──────────────────────────────────────────────────────────────────
INPUT_FILE  = "//ggd01/Seismik/2026_LVL_Erfurt/Bearbeitung/PICKS/job.output"
OUTPUT_FILE = "//ggd01/Seismik/2026_LVL_Erfurt/Bearbeitung/PICKS/715_picks.txt"

# Column header that marks data lines (regex, case-sensitive)
HEADER_PATTERN = r"Ensemble\s+#\s+SOURCE\s+CHAN\s+OFFSET\s+FB_PICK"

# A data line starts with optional spaces then an integer (the ensemble number)
DATA_LINE_PATTERN = r"^\s+\d"
# ─────────────────────────────────────────────────────────────────────────────


def parse_picks(path):
    """Return the header string and a list of (source_id, raw_line) tuples."""
    header = None
    rows   = []          # list of (source_id, line)

    with open(path, "r") as fh:
        for raw in fh:
            line = raw.rstrip("\n")

            # Capture the first header line encountered
            if re.search(HEADER_PATTERN, line):
                if header is None:
                    header = line
                continue  # skip repeated headers

            # Keep only data lines
            if re.match(DATA_LINE_PATTERN, line):
                cols      = line.split()
                source_id = int(cols[1])          # SOURCE is column index 1
                rows.append((source_id, line))

    return header, rows


def write_clean(header, rows, out_path):
    """Write one header followed by data rows, with a blank line between sources."""
    with open(out_path, "w") as fh:
        fh.write(header + "\n")

        prev_source = None
        for source_id, line in rows:
            # Insert blank separator when the SOURCE value changes
            if prev_source is not None and source_id != prev_source:
                fh.write("\n")
            fh.write(line + "\n")
            prev_source = source_id


def main():
    header, rows = parse_picks(INPUT_FILE)

    if header is None:
        raise ValueError("No header line found – check HEADER_PATTERN.")
    if not rows:
        raise ValueError("No data lines found – check DATA_LINE_PATTERN.")

    write_clean(header, rows, OUTPUT_FILE)
    print(f"Done. {len(rows)} rows written to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()