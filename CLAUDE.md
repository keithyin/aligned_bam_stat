# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`bam_stat` is a Rust CLI tool for computing per-read alignment concordance metrics from BAM files. The main use case is evaluating alignment quality by counting matches, mismatches, homopolymer indels, and other CIGAR-derived statistics per BAM record.

## Commands

```bash
# Build
cargo build

# Run tests
cargo test

# Run a specific test
cargo test test_name

# Run the CLI (single subcommand: bam-concordance)
cargo run -- bam-concordance <reffasta> <aligned_bam> [--hcregions <bed>] [--hcvariants <vcf>] [--chrom <chrom>]
```

## Architecture

### Entry point (`src/main.rs`)

Routes the single subcommand `BamConcordance` to `bam_concordance::bam_concordance()`. Uses `clap` derive for CLI parsing.

### CLI (`src/cli.rs`)

Defines the `BamConcordanceArgs` struct:
- Positional: `reffasta` (reference FASTA), `aligned_bam` (input BAM)
- Optional: `--hcregions` (BED file to restrict concordance regions), `--hcvariants` (VCF file to exclude variant sites), `--chrom` (unused filter)

### Core logic (`src/bam_concordance.rs`)

The main pipeline uses a threaded producer-consumer pattern via `crossbeam`:

1. **Preload phase** (outside scope): Loads HC regions (BED) and HC variants (VCF) in separate threads, then constructs a `FastaFile` (full reference genome loaded into memory as `HashMap<refname, seq>`).

2. **Producer thread** (inside `thread_scope`): Reads the BAM file using `rust_htslib::bam` with 10 threads, wraps each `bam::Record` into a `RecordReplica`, and sends it through a bounded channel.

3. **Consumer threads** (8 workers): Each pulls records from the channel and runs `stat_record_core()` which iterates the CIGAR string, comparing against the reference sequence to count:
   - `m` (match/equal bases)
   - `mm` (mismatch/diff bases)
   - `h_ins` / `non_h_ins` (homopolymer vs non-homopolymer insertions)
   - `h_del` / `non_h_del` (homopolymer vs non-homopolymer deletions)
   - `ignore_bps` (bases in excluded regions or at variant sites)

4. **Output**: Writes per-read stats to `{bam_basename}.metric.csv` as TSV with columns for all metrics plus `mmRefPositions`, `insRefPositions`, `delRefPositions` (semicolon-separated positions).

Key data structures:
- `RecordReplica` - deserialized BAM record (extracted fields, no borrow from htslib)
- `Stat` - computed alignment statistics for a single read

### Common modules (`src/common/`)

- **`bam_ext.rs`** - `BamRecordExt` wraps `bam::Record` to extract custom tags (`ec`, `iy`, `rq`, `ch`, `np`, `dw`, `cr`) via `aux()` parsing. Also `AlignedRecord` with helper methods `compute_effective_coverage()` and `compute_identity()`.

- **`file_reader/mod.rs`** - Re-exports four file parsers, each with a dedicated iterator and a compiled-index struct:
  - `fasta_reader.rs` - `FastaFile` loads full FASTA into `HashMap<name, seq>` in memory
  - `bed_reader.rs` - `BedInfo` stores intervals per chromosome as sorted `Vec<(start, end)>` for binary search lookups
  - `vcf_reader.rs` - `VcfInfo` stores variant positions per chromosome as sorted `Vec<usize>` for binary search lookups
  - `fastq_reader.rs` - `FastqReaderIter` iterator (unused by main pipeline)

- **`pb_tools.rs`** - Progress bar helpers (`get_spin_pb`, `get_bar_pb`) using `indicatif`.

## Key Dependencies

- `rust-htslib` - BAM/FASTA/HTS file I/O
- `bio` - CIGAR/alignment utilities
- `polars` - CSV/data frame processing (imported but currently unused in the pipeline; outputs raw CSV via `BufWriter`)
- `rayon` + `indicatif[rayon]` - parallelism and progress bars
- `clap` - CLI argument parsing
- `crossbeam` - thread channels for the producer-consumer pipeline
