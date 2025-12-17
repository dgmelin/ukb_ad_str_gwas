import gzip
import os
import csv
import argparse

addon_lines = [
    # "##INFO=<ID=END,Number=1,Type=Integer,Description=\"ending position\">",
    "##INFO=<ID=START,Number=1,Type=Integer,Description=\"starting position\">",
    "##INFO=<ID=PERIOD,Number=1,Type=Integer,Description=\"period\">"
]

def log(message, verbose = True):
    if verbose:
        print(message)

def run():
    parser = argparse.ArgumentParser(description="Annotate imputed VCF files with STR information.")
    parser.add_argument('input', type=str, help='Input VCF file (gzipped)')
    parser.add_argument('output', type=str, help='Output file path')
    args = parser.parse_args()

    infile = args.input
    outfile = args.output

    if not os.path.exists(infile):
        log(f"Input file {infile} does not exist.", True)
        exit(1)

    log(f"Processing file: {infile}")
    outdir = os.path.dirname(outfile)
    os.makedirs(outdir, exist_ok=True)
    
    run_linewise(infile, outfile)
    log('Finished annotating imputed VCF file.')

def run_linewise(in_path, out_path, verbose=True):
    with gzip.open(in_path, 'rt') as infile, gzip.open(out_path, 'wt') as outfile:
        for ln in infile:
            if ln.startswith('#'):
                outfile.write(ln)
                if 'ID=END' in ln:
                    for add_ln in addon_lines:
                        outfile.write(add_ln + "\n")
                continue

            row = ln.strip().split('\t')
            start = int(row[1])
            ref = row[3]
            end = start + len(ref) - 1
            period, start_offset = find_pattern(ref)
            alt_pers = [find_pattern(alt)[0] for alt in row[4].split(',')]
            min_period = min(alt_pers + [period])
            start = start + start_offset
            info = row[7]
            if 'END=' in info:
                row[7] = f"START={start};PERIOD={min_period};{info}"
            else:
                row[7] = f"START={start};END={end};PERIOD={min_period};{info}"

            outfile.write('\t'.join(row) + '\n')


def find_pattern(seq: str) -> int:
    s = seq.upper()
    n = len(s)
    best_len = 0
    best_k = 0

    # Scan all start positions
    for i in range(n):
        # For a tandem repeat with period k, we need at least two copies:
        # i + 2*k <= n
        max_k = (n - i) // 2
        for k in range(1, max_k + 1):
            block = s[i:i + k]
            # Count how many times block repeats contiguously
            r = 0
            while i + (r + 1) * k <= n and s[i + r * k:i + (r + 1) * k] == block:
                r += 1
            if r >= 2:
                rep_len = r * k
                # Prefer longer repeated region; break ties by smaller period
                if rep_len > best_len or (rep_len == best_len and (best_k == 0 or k < best_k)):
                    best_len = rep_len
                    best_k = k
                    best_i = i

    best_k = best_k if best_len > 1 else 1
    best_i = best_i if best_len > 1 else 0
    return best_k, best_i



if __name__ == "__main__":
    run()