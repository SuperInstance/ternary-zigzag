# ternary-zigzag

Zigzag scanning for ternary matrix compression.

[![crates.io](https://img.shields.io/crates/v/ternary-zigzag.svg)](https://crates.io/crates/ternary-zigzag)

## Overview

Zigzag scanning reorders 2D matrix elements in a diagonal serpentine pattern that
groups low-frequency coefficients together — the same technique used in JPEG image
compression. After transform coding (DCT, wavelets, Walsh), most signal energy
concentrates in the top-left corner of the coefficient matrix. Zigzag scanning
exploits this by traversing low-frequency coefficients first, producing runs of
zeros that compress extremely well with run-length encoding (RLE).

This crate provides:

- **Zigzag scanning** — Convert 2D NxN matrices to 1D sequences in zigzag order
- **Inverse zigzag** — Reconstruct the original matrix from a scanned sequence
- **Block-based scanning** — Scan individual blocks from a larger matrix
- **Run-length encoding** — Compress zigzag output with RLE for efficient storage
- **Full matrix block scanning** — Break a matrix into blocks and scan each independently
- **Index generation** — Get zigzag traversal order as (row, col) pairs

## Installation

Add to your `Cargo.toml`:

```toml
[dependencies]
ternary-zigzag = "0.1.0"
```

## Quick Start

```rust
use ternary_zigzag::*;

// A 4x4 ternary coefficient matrix (e.g., after wavelet transform)
let matrix = vec![
    vec![7, 0, 0, 0],
    vec![0, 0, -3, 0],
    vec![0, 0, 0, 0],
    vec![0, 2, 0, 0],
];

// Scan in zigzag order
let scanned = zigzag_scan(&matrix);
// Result: [7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0, 2, 0, 0]

// Run-length encode for compression
let encoded = rle_encode(&scanned);
// Many consecutive zeros compress to a single (value, count) pair

// Decode and reconstruct
let decoded = rle_decode(&encoded);
let recovered = inverse_zigzag(&decoded, 4);
assert_eq!(recovered, matrix);

// Block-based scanning (e.g., 8x8 JPEG-style blocks)
let block = block_zigzag(&matrix, 0, 0, 2);
```

## Mathematical Background

### Zigzag Scanning

For an NxN matrix, zigzag scanning traverses elements along anti-diagonals,
alternating between upward and downward directions:

```text
4x4 traversal order:
 0  1  5  6
 2  4  7 12
 3  8 11 13
 9 10 14 15
```

This order ensures that coefficients are visited in approximately increasing
"frequency" order (top-left = low frequency, bottom-right = high frequency).

### Why Zigzag for Compression?

After applying a decorrelating transform (DCT, wavelets, Walsh), coefficient
matrices have a characteristic energy distribution:

1. **Most energy** is concentrated in the top-left (low-frequency) corner
2. **Most coefficients** in the bottom-right are zero or near-zero
3. **Zigzag ordering** groups non-zero coefficients early and zeros late
4. **RLE** then efficiently encodes the long runs of trailing zeros

For ternary signals specifically, the coefficient matrices tend to be even more
sparse because the original signal only takes three values, producing fewer
distinct transform coefficients.

### Block-Based Processing

In practice (e.g., JPEG), images are divided into 8×8 blocks. Each block is:

1. Transformed (DCT/wavelet)
2. Zigzag scanned
3. Run-length encoded
4. Entropy coded (Huffman/arithmetic)

This crate supports this workflow through `block_zigzag` and `full_block_zigzag`.

## API Reference

### Functions

| Function | Description |
|----------|-------------|
| `zigzag_scan(matrix)` | Scan NxN matrix in zigzag order → Vec |
| `inverse_zigzag(flat, n)` | Reconstruct NxN matrix from zigzag sequence |
| `block_zigzag(matrix, row, col, size)` | Scan a block from a larger matrix |
| `inverse_block_zigzag(flat, n, dest, row, col, size)` | Place reconstructed block into matrix |
| `full_block_zigzag(matrix, block_size)` | Scan entire matrix block-by-block |
| `zigzag_indices(n)` | Get zigzag order as (row, col) index pairs |
| `rle_encode(data)` | Run-length encode any sequence |
| `rle_decode(pairs)` | Decode RLE back to original sequence |
| `rle_compression_ratio(pairs, original_len)` | Compute (encoded_size, original_size, ratio) |

### Types

| Type | Description |
|------|-------------|
| `RLEPair<T>` | Run-length pair with `value: T` and `count: usize` |

## Properties Verified by Tests

- **Round-trip**: zigzag_scan → inverse_zigzag recovers the original matrix
- **All positions covered**: Every (row, col) pair appears exactly once
- **No duplicates**: No position is visited twice
- **Known patterns**: Correct ordering verified against standard JPEG zigzag
- **RLE round-trip**: encode → decode recovers the original sequence
- **RLE compression**: All-same data compresses to 1 pair; all-different stays same size
- **Block scanning**: Blocks correctly extracted and scanned from larger matrices
- **Full workflow**: Matrix → zigzag → RLE → decode → inverse zigzag = original
- **Sizes**: Works for 1×1, 2×2, 3×3, 4×4, 8×8, 16×16 matrices
- **Out-of-bounds safety**: Block scanning panics on out-of-bounds access

## Use Cases

- **JPEG-like compression** for ternary images/matrices
- **Wavelet coefficient reordering** before entropy coding
- **DCT coefficient serialization** for transform coding pipelines
- **Sparse matrix compression** for ternary-valued grids
- **Pre-processing for Huffman/arithmetic coding**

## License

MIT
