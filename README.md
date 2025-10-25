# Huffman Compression Project

This repository implements a **Huffman encoding/decoding system** in C, with support for **serial**, **OpenMP**, and **Pthreads** parallelization. It includes performance tests on multiple sample files and visualization scripts.

## Directory Structure

- `src/` – Source code for serial and parallel Huffman compressors.
- `samples/` – Sample input files of different sizes and types.
- `scripts/` – Helper scripts:
  - `run_performance.sh` – Runs compression tests and records timings.
  - `plot_compression.py` – Plots compression timings using Matplotlib.
- `results/` – Stores output files and timing results.
- `.gitignore` – Ignore compiled binaries and temporary files.

## Compilation

```bash
cd src
gcc -O2 -Wall huffman_serial.c -o huffman_serial
gcc -O2 -Wall huffman_parallel.c -fopenmp -pthread -o huffman_parallel
