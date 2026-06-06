//! # ternary-zigzag
//!
//! Zigzag scanning for ternary matrix compression.
//!
//! Zigzag scanning reorders 2D matrix elements in a diagonal zigzag pattern,
//! which groups low-frequency coefficients together — critical for efficient
//! compression of transform-domain data (e.g., after DCT or wavelet transforms).
//! This crate implements zigzag scanning, inverse reconstruction, block-based
//! scanning, and run-length encoding for ternary-valued matrices.

/// Scan a 2D matrix in zigzag (diagonal) order.
///
/// The zigzag pattern traverses the matrix diagonally, alternating between
/// upward-right and downward-left directions. This is the same scanning order
/// used in JPEG compression.
///
/// For an NxN matrix, produces an N*N element vector.
///
/// # Arguments
/// * `matrix` - A square matrix as Vec<Vec<T>> (must be non-empty and square)
///
/// # Example order for 4x4:
/// ```text
/// 0  1  5  6
/// 2  4  7  12
/// 3  8  11 13
/// 9  10 14 15
/// ```
pub fn zigzag_scan<T: Clone>(matrix: &[Vec<T>]) -> Vec<T> {
    let n = matrix.len();
    assert!(!matrix.is_empty(), "matrix must not be empty");
    assert!(
        matrix.iter().all(|row| row.len() == n),
        "matrix must be square"
    );

    let mut result = Vec::with_capacity(n * n);
    let mut row = 0usize;
    let mut col = 0usize;
    let mut going_up = true;

    for _ in 0..(n * n) {
        result.push(matrix[row][col].clone());

        if going_up {
            if col == n - 1 {
                row += 1;
                going_up = false;
            } else if row == 0 {
                col += 1;
                going_up = false;
            } else {
                row -= 1;
                col += 1;
            }
        } else {
            if row == n - 1 {
                col += 1;
                going_up = true;
            } else if col == 0 {
                row += 1;
                going_up = true;
            } else {
                row += 1;
                col -= 1;
            }
        }
    }

    result
}

/// Reconstruct a 2D matrix from a zigzag-scanned 1D sequence.
///
/// This is the inverse of `zigzag_scan`. Given a flat vector of N*N elements,
/// reconstructs the original NxN matrix by placing elements in zigzag order.
pub fn inverse_zigzag<T: Clone + Default>(flat: &[T], n: usize) -> Vec<Vec<T>> {
    assert_eq!(flat.len(), n * n, "flat length must be n*n");

    let mut matrix = vec![vec![T::default(); n]; n];
    let mut row = 0usize;
    let mut col = 0usize;
    let mut going_up = true;

    for i in 0..(n * n) {
        matrix[row][col] = flat[i].clone();

        if going_up {
            if col == n - 1 {
                row += 1;
                going_up = false;
            } else if row == 0 {
                col += 1;
                going_up = false;
            } else {
                row -= 1;
                col += 1;
            }
        } else {
            if row == n - 1 {
                col += 1;
                going_up = true;
            } else if col == 0 {
                row += 1;
                going_up = true;
            } else {
                row += 1;
                col -= 1;
            }
        }
    }

    matrix
}

/// Scan an NxN block from a larger matrix at the given offset.
///
/// Extracts a block of size `block_size` starting at (row_offset, col_offset)
/// from the matrix, then scans it in zigzag order.
///
/// # Panics
/// Panics if the block would extend beyond the matrix boundaries.
pub fn block_zigzag<T: Clone>(
    matrix: &[Vec<T>],
    row_offset: usize,
    col_offset: usize,
    block_size: usize,
) -> Vec<T> {
    let n = matrix.len();
    assert!(row_offset + block_size <= n, "block extends beyond matrix rows");
    assert!(
        col_offset + block_size <= matrix.get(0).map_or(0, |r| r.len()),
        "block extends beyond matrix columns"
    );

    let block: Vec<Vec<T>> = (row_offset..row_offset + block_size)
        .map(|i| matrix[i][col_offset..col_offset + block_size].to_vec())
        .collect();

    zigzag_scan(&block)
}

/// Reconstruct a block from zigzag-scanned data and place it in a matrix.
///
/// Creates a block of size `block_size` from the zigzag-scanned data and places
/// it into the destination matrix at (row_offset, col_offset).
pub fn inverse_block_zigzag<T: Clone + Default>(
    flat: &[T],
    n: usize,
    dest: &mut [Vec<T>],
    row_offset: usize,
    col_offset: usize,
    block_size: usize,
) {
    let block = inverse_zigzag(flat, block_size);
    for i in 0..block_size {
        for j in 0..block_size {
            dest[row_offset + i][col_offset + j] = block[i][j].clone();
        }
    }
}

/// A run-length encoded pair: (value, count).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RLEPair<T> {
    pub value: T,
    pub count: usize,
}

/// Run-length encode a sequence.
///
/// Consecutive equal values are grouped into (value, count) pairs.
/// This is particularly effective on zigzag-scanned ternary data,
/// where many zeros (or other values) appear in long runs.
pub fn rle_encode<T: Clone + PartialEq>(data: &[T]) -> Vec<RLEPair<T>> {
    if data.is_empty() {
        return Vec::new();
    }

    let mut result = Vec::new();
    let mut current = data[0].clone();
    let mut count = 1;

    for i in 1..data.len() {
        if data[i] == current {
            count += 1;
        } else {
            result.push(RLEPair {
                value: current,
                count,
            });
            current = data[i].clone();
            count = 1;
        }
    }

    result.push(RLEPair {
        value: current,
        count,
    });
    result
}

/// Decode run-length encoded data back to a flat sequence.
pub fn rle_decode<T: Clone>(pairs: &[RLEPair<T>]) -> Vec<T> {
    let mut result = Vec::new();
    for pair in pairs {
        for _ in 0..pair.count {
            result.push(pair.value.clone());
        }
    }
    result
}

/// Compute the compression ratio of RLE encoding.
///
/// Returns (encoded_size, original_size, ratio) where ratio < 1.0 means compression.
pub fn rle_compression_ratio<T>(pairs: &[RLEPair<T>], original_len: usize) -> (usize, usize, f64) {
    let encoded_size = pairs.len();
    let ratio = encoded_size as f64 / original_len.max(1) as f64;
    (encoded_size, original_len, ratio)
}

/// Zigzag scan a full matrix by breaking it into blocks.
///
/// Divides the matrix into NxN blocks (where N = block_size) and scans each
/// block in zigzag order. Blocks are processed in row-major order.
/// If the matrix dimensions are not divisible by block_size, the remaining
/// elements are handled with a smaller final block (or padding could be applied externally).
///
/// Returns a flat vector of all scanned elements in block order.
pub fn full_block_zigzag<T: Clone + Default>(
    matrix: &[Vec<T>],
    block_size: usize,
) -> Vec<T> {
    let rows = matrix.len();
    let cols = matrix.get(0).map_or(0, |r| r.len());
    let mut result = Vec::new();

    let mut r = 0;
    while r < rows {
        let mut c = 0;
        let actual_block_r = block_size.min(rows - r);
        while c < cols {
            let actual_block_c = block_size.min(cols - c);

            if actual_block_r == block_size && actual_block_c == block_size {
                // Full block
                let scanned = block_zigzag(matrix, r, c, block_size);
                result.extend(scanned);
            } else {
                // Partial block — extract and scan
                let block: Vec<Vec<T>> = (r..r + actual_block_r)
                    .map(|i| matrix[i][c..c + actual_block_c].to_vec())
                    .collect();
                result.extend(zigzag_scan(&block));
            }
            c += block_size;
        }
        r += block_size;
    }

    result
}

/// Generate the zigzag scan order as indices for an NxN matrix.
///
/// Returns a vector of (row, col) pairs in zigzag order.
pub fn zigzag_indices(n: usize) -> Vec<(usize, usize)> {
    let mut indices = Vec::with_capacity(n * n);
    let mut row = 0usize;
    let mut col = 0usize;
    let mut going_up = true;

    for _ in 0..(n * n) {
        indices.push((row, col));

        if going_up {
            if col == n - 1 {
                row += 1;
                going_up = false;
            } else if row == 0 {
                col += 1;
                going_up = false;
            } else {
                row -= 1;
                col += 1;
            }
        } else {
            if row == n - 1 {
                col += 1;
                going_up = true;
            } else if col == 0 {
                row += 1;
                going_up = true;
            } else {
                row += 1;
                col -= 1;
            }
        }
    }

    indices
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_zigzag_1x1() {
        let matrix = vec![vec![42]];
        let result = zigzag_scan(&matrix);
        assert_eq!(result, vec![42]);
    }

    #[test]
    fn test_zigzag_2x2() {
        let matrix = vec![vec![1, 2], vec![3, 4]];
        let result = zigzag_scan(&matrix);
        // Order: (0,0), (0,1), (1,0), (1,1)
        assert_eq!(result, vec![1, 2, 3, 4]);
    }

    #[test]
    fn test_zigzag_3x3() {
        let matrix = vec![
            vec![1, 2, 3],
            vec![4, 5, 6],
            vec![7, 8, 9],
        ];
        let result = zigzag_scan(&matrix);
        // Zigzag order: 1,2,4,7,5,3,6,8,9
        assert_eq!(result, vec![1, 2, 4, 7, 5, 3, 6, 8, 9]);
    }

    #[test]
    fn test_zigzag_4x4() {
        let matrix: Vec<Vec<i32>> = (0..4)
            .map(|i| (0..4).map(|j| i * 4 + j).collect())
            .collect();
        let result = zigzag_scan(&matrix);
        // Standard JPEG zigzag for 4x4:
        // 0 1 4 5 / 2 3 6 7 / 8 9 10 11 / 12 13 14 15
        // should map to: 0,1,4,8,5,2,3,6,9,12,13,10,7,11,14,15
        assert_eq!(result, vec![0, 1, 4, 8, 5, 2, 3, 6, 9, 12, 13, 10, 7, 11, 14, 15]);
    }

    #[test]
    fn test_zigzag_8x8_known_pattern() {
        // Fill with sequential values and verify specific positions
        let matrix: Vec<Vec<i32>> = (0..8)
            .map(|i| (0..8).map(|j| i * 8 + j).collect())
            .collect();
        let result = zigzag_scan(&matrix);

        assert_eq!(result.len(), 64);
        // First element should be top-left (0,0)
        assert_eq!(result[0], 0);
        // Last element should be bottom-right (7,7)
        assert_eq!(result[63], 63);
        // Second element should be (0,1)
        assert_eq!(result[1], 1);
    }

    #[test]
    fn test_inverse_zigzag_roundtrip() {
        let matrices = vec![
            vec![vec![1]],
            vec![vec![1, 2], vec![3, 4]],
            vec![vec![1, 2, 3], vec![4, 5, 6], vec![7, 8, 9]],
            (0..4).map(|i| (0..4).map(|j| i * 4 + j + 1).collect()).collect(),
            (0..8).map(|i| (0..8).map(|j| i * 8 + j).collect()).collect(),
        ];

        for matrix in &matrices {
            let n = matrix.len();
            let scanned = zigzag_scan(matrix);
            let recovered = inverse_zigzag(&scanned, n);

            assert_eq!(*matrix, recovered, "Roundtrip failed for {}x{}", n, n);
        }
    }

    #[test]
    fn test_inverse_zigzag_ternary() {
        let matrix = vec![
            vec![1, -1, 0, 1],
            vec![0, 0, 1, -1],
            vec![-1, 1, 0, 0],
            vec![1, 0, -1, 1],
        ];
        let scanned = zigzag_scan(&matrix);
        let recovered = inverse_zigzag(&scanned, 4);
        assert_eq!(matrix, recovered);
    }

    #[test]
    fn test_block_zigzag() {
        let matrix = vec![
            vec![1, 2, 3, 4],
            vec![5, 6, 7, 8],
            vec![9, 10, 11, 12],
            vec![13, 14, 15, 16],
        ];

        // Top-left 2x2 block
        let block = block_zigzag(&matrix, 0, 0, 2);
        assert_eq!(block, vec![1, 2, 5, 6]);

        // Bottom-right 2x2 block
        let block = block_zigzag(&matrix, 2, 2, 2);
        assert_eq!(block, vec![11, 12, 15, 16]);

        // Center-ish block
        let block = block_zigzag(&matrix, 1, 1, 2);
        assert_eq!(block, vec![6, 7, 10, 11]);
    }

    #[test]
    fn test_inverse_block_zigzag() {
        let scanned = vec![1, 2, 5, 6];
        let mut dest = vec![vec![0; 4]; 4];
        inverse_block_zigzag(&scanned, 2, &mut dest, 0, 0, 2);

        assert_eq!(dest[0][0], 1);
        assert_eq!(dest[0][1], 2);
        assert_eq!(dest[1][0], 5);
        assert_eq!(dest[1][1], 6);
        // Rest should be 0
        assert_eq!(dest[2][2], 0);
    }

    #[test]
    fn test_rle_encode_basic() {
        let data = vec![1, 1, 1, 2, 2, 3];
        let encoded = rle_encode(&data);
        assert_eq!(encoded.len(), 3);
        assert_eq!(encoded[0], RLEPair { value: 1, count: 3 });
        assert_eq!(encoded[1], RLEPair { value: 2, count: 2 });
        assert_eq!(encoded[2], RLEPair { value: 3, count: 1 });
    }

    #[test]
    fn test_rle_encode_ternary() {
        let data = vec![0, 0, 0, 0, 1, -1, -1, 0, 0];
        let encoded = rle_encode(&data);
        assert_eq!(encoded.len(), 4);
        assert_eq!(encoded[0], RLEPair { value: 0, count: 4 });
        assert_eq!(encoded[1], RLEPair { value: 1, count: 1 });
        assert_eq!(encoded[2], RLEPair { value: -1, count: 2 });
        assert_eq!(encoded[3], RLEPair { value: 0, count: 2 });
    }

    #[test]
    fn test_rle_decode_roundtrip() {
        let data_sets = vec![
            vec![1, 1, 1, 2, 2, 3],
            vec![0, 0, 0, 0, 0],
            vec![1],
            vec![1, -1, 1, -1],
            vec![0i32; 100],
        ];

        for data in &data_sets {
            let encoded = rle_encode(data);
            let decoded = rle_decode(&encoded);
            assert_eq!(*data, decoded);
        }
    }

    #[test]
    fn test_rle_empty() {
        let data: Vec<i32> = vec![];
        let encoded = rle_encode(&data);
        assert!(encoded.is_empty());
        let decoded = rle_decode(&encoded);
        assert!(decoded.is_empty());
    }

    #[test]
    fn test_rle_compression_ratio() {
        // All same -> excellent compression
        let data = vec![0i32; 100];
        let encoded = rle_encode(&data);
        let (enc_size, orig_size, ratio) = rle_compression_ratio(&encoded, data.len());
        assert_eq!(enc_size, 1);
        assert_eq!(orig_size, 100);
        assert!((ratio - 0.01).abs() < 1e-10);

        // All different -> no compression
        let data: Vec<i32> = (0..50).collect();
        let encoded = rle_encode(&data);
        let (enc_size, _, ratio) = rle_compression_ratio(&encoded, data.len());
        assert_eq!(enc_size, 50);
        assert!((ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_zigzag_ternary_compression() {
        // Simulate a transformed ternary matrix (mostly zeros, a few non-zero)
        let matrix = vec![
            vec![7, 0, 0, 0],
            vec![0, 0, -3, 0],
            vec![0, 0, 0, 0],
            vec![0, 2, 0, 0],
        ];

        let scanned = zigzag_scan(&matrix);
        let encoded = rle_encode(&scanned);
        let decoded = rle_decode(&encoded);

        assert_eq!(scanned, decoded);
        // Should compress well due to many zeros
        assert!(encoded.len() < scanned.len());
    }

    #[test]
    fn test_full_block_zigzag() {
        let matrix = vec![
            vec![1, 2, 3, 4],
            vec![5, 6, 7, 8],
            vec![9, 10, 11, 12],
            vec![13, 14, 15, 16],
        ];

        let result = full_block_zigzag(&matrix, 2);
        // Four 2x2 blocks, each scanned in zigzag
        assert_eq!(result.len(), 16);

        // First block (top-left): 1,2,5,6
        assert_eq!(result[0..4], [1, 2, 5, 6]);
        // Second block (top-right): 3,4,7,8
        assert_eq!(result[4..8], [3, 4, 7, 8]);
        // Third block (bottom-left): 9,10,13,14
        assert_eq!(result[8..12], [9, 10, 13, 14]);
        // Fourth block (bottom-right): 11,12,15,16
        assert_eq!(result[12..16], [11, 12, 15, 16]);
    }

    #[test]
    fn test_zigzag_indices() {
        let indices = zigzag_indices(3);
        assert_eq!(indices.len(), 9);
        assert_eq!(indices[0], (0, 0));
        assert_eq!(indices[1], (0, 1));
        assert_eq!(indices[2], (1, 0));
        assert_eq!(indices[8], (2, 2));

        // Verify all indices are present
        let mut all_pairs: Vec<(usize, usize)> = (0..3).flat_map(|i| (0..3).map(move |j| (i, j))).collect();
        all_pairs.sort();
        let mut sorted_indices = indices.clone();
        sorted_indices.sort();
        assert_eq!(sorted_indices, all_pairs);
    }

    #[test]
    fn test_zigzag_indices_no_duplicates() {
        for &n in &[1, 2, 3, 4, 8, 16] {
            let indices = zigzag_indices(n);
            assert_eq!(indices.len(), n * n);
            let mut seen = std::collections::HashSet::new();
            for &(r, c) in &indices {
                assert!(
                    seen.insert((r, c)),
                    "Duplicate index ({}, {}) for n={}",
                    r,
                    c,
                    n
                );
            }
        }
    }

    #[test]
    fn test_zigzag_all_positions_covered() {
        for &n in &[1, 2, 3, 4, 5, 8] {
            let matrix: Vec<Vec<usize>> = (0..n)
                .map(|i| (0..n).map(|j| i * n + j).collect())
                .collect();
            let scanned = zigzag_scan(&matrix);

            // All n*n values should be present exactly once
            let mut sorted = scanned.clone();
            sorted.sort();
            let expected: Vec<usize> = (0..(n * n)).collect();
            assert_eq!(sorted, expected, "Missing values for n={}", n);
        }
    }

    #[test]
    fn test_rle_with_zigzag_workflow() {
        // Full workflow: matrix -> zigzag -> RLE -> decode -> inverse zigzag
        let original = vec![
            vec![5, 0, 0, 0],
            vec![0, 0, 0, 0],
            vec![0, 0, -3, 0],
            vec![0, 0, 0, 0],
        ];

        let scanned = zigzag_scan(&original);
        let encoded = rle_encode(&scanned);
        let decoded = rle_decode(&encoded);
        assert_eq!(decoded, scanned);

        let recovered = inverse_zigzag(&decoded, 4);
        assert_eq!(recovered, original);
    }

    #[test]
    fn test_known_jpeg_zigzag_8x8() {
        // Verify the first few and last few indices of the standard 8x8 zigzag
        let indices = zigzag_indices(8);
        // Standard JPEG zigzag order:
        assert_eq!(indices[0], (0, 0));
        assert_eq!(indices[1], (0, 1));
        assert_eq!(indices[2], (1, 0));
        assert_eq!(indices[3], (2, 0));
        assert_eq!(indices[4], (1, 1));
        assert_eq!(indices[5], (0, 2));
        assert_eq!(indices[6], (0, 3));
        assert_eq!(indices[7], (1, 2));
        // Last
        assert_eq!(indices[63], (7, 7));
        assert_eq!(indices[62], (7, 6));
    }

    #[test]
    fn test_large_matrix_zigzag() {
        let n = 16;
        let matrix: Vec<Vec<i32>> = (0..n)
            .map(|i| (0..n).map(|j| (i * n + j) as i32).collect())
            .collect();

        let scanned = zigzag_scan(&matrix);
        assert_eq!(scanned.len(), n * n);

        let recovered = inverse_zigzag(&scanned, n);
        assert_eq!(recovered, matrix);
    }

    #[test]
    fn test_rle_single_element() {
        let data = vec![42];
        let encoded = rle_encode(&data);
        assert_eq!(encoded.len(), 1);
        assert_eq!(encoded[0], RLEPair { value: 42, count: 1 });
        let decoded = rle_decode(&encoded);
        assert_eq!(decoded, data);
    }

    #[test]
    fn test_block_zigzag_out_of_bounds() {
        let matrix = vec![vec![1, 2], vec![3, 4]];
        // This should panic because block extends beyond matrix
        let result = std::panic::catch_unwind(|| {
            block_zigzag(&matrix, 1, 1, 2);
        });
        assert!(result.is_err());
    }
}
