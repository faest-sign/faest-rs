use crate::field::{GF2View, GF2p8};
use ndarray::{Array1, Array2, ArrayView1, ArrayView2, Axis};

type GF2p8Vector = Array1<GF2p8>;
type GF2p8VectorView<'a> = ArrayView1<'a, GF2p8>;
type GF2p8Matrix = Array2<GF2p8>;
type GF2p8MatrixView<'a> = ArrayView2<'a, GF2p8>;

pub fn hash_bitvector(r: GF2p8VectorView, v: &GF2View) -> GF2p8Vector {
    let tau = r.len();
    let mut r_powers = r.to_owned();
    let mut h = if v[0] {
        GF2p8Vector::ones(tau)
    } else {
        GF2p8Vector::zeros(tau)
    };
    if v[1] {
        h += &r;
    }
    for v_i in v.iter().by_vals().skip(2) {
        r_powers *= &r;
        if v_i {
            h += &r_powers;
        }
    }

    h
}

#[allow(non_snake_case)]
pub fn hash_matrix(r: GF2p8VectorView, M: GF2p8MatrixView) -> GF2p8Matrix {
    let tau = r.len();
    let n = M.shape()[0];
    let m = M.shape()[1];
    // compute the matrix product H := RM where R is the (tau x n) Vandermonde matrix which has
    // successive powers of r as columns

    let mut r_powers = r.to_owned();
    let mut H = M.row(0).broadcast((tau, m)).unwrap().to_owned();
    for (&rp_i, mut h_i) in r_powers.iter().zip(H.axis_iter_mut(Axis(0))) {
        h_i.scaled_add(rp_i, &M.row(1));
    }
    for i in 2..n {
        r_powers *= &r;
        for (&rp_i, mut h_i) in r_powers.iter().zip(H.axis_iter_mut(Axis(0))) {
            h_i.scaled_add(rp_i, &M.row(i));
        }
    }
    H
}

#[allow(non_snake_case)]
pub fn hash_bitvector_and_matrix(
    r: GF2p8VectorView,
    v: &GF2View,
    M: GF2p8MatrixView,
) -> (GF2p8Vector, GF2p8Matrix) {
    // TODO: optimize
    (hash_bitvector(r, v), hash_matrix(r, M))
}

const fn gen_mask_table() -> [u64; 256] {
    const fn f(x: u8) -> u64 {
        let mut m = 0;
        let mut i = 0;
        while i < 8 {
            if x & (1 << i) != 0 {
                m |= 0xff << (i * 8);
            }
            i += 1;
        }
        m
    }
    let mut table = [0; 256];
    let mut x = 0;
    while x < 256 {
        table[x] = f(x as u8);
        x += 1;
    }
    table
}

const MASK_TABLE: [u64; 256] = gen_mask_table();

/// Compute y_i += x * b_i for i = 0..
pub fn bitmul_accumulate(ys: &mut [GF2p8], x: GF2p8, bs: &[u8]) {
    let n = ys.len();
    debug_assert!(
        (bs.len() * 8 - 7 <= n) && (n <= bs.len() * 8),
        "ys.len() = {}, but bs.len() * 8 = {}",
        ys.len(),
        bs.len() * 8
    );

    let n_blocks = n / 16;

    let xs = u128::from_le_bytes([x.0; 16]);
    // This cast is fine, since GF2p8 is wrapping a u8 and we have computed the number of 16B
    // blocks in the array.
    let xmm_slice =
        unsafe { std::slice::from_raw_parts_mut(ys.as_mut_ptr() as *mut u128, n_blocks) };
    for j in 0..n_blocks {
        let mask_lo = MASK_TABLE[bs[2 * j] as usize] as u128;
        let mask_hi = MASK_TABLE[bs[2 * j + 1] as usize] as u128;
        let mask = mask_hi << 64 | mask_lo;
        xmm_slice[j] ^= mask & xs;
    }

    if n & 8 != 0 {
        let xs = u64::from_le_bytes([x.0; 8]);
        // This cast is also fine, since the remainder of the array is more than 8B big.
        let q_slice = unsafe {
            std::slice::from_raw_parts_mut(ys.as_mut_ptr() as *mut u64, 2 * n_blocks + 1)
        };
        let mask = MASK_TABLE[bs[2 * n_blocks] as usize];
        q_slice[2 * n_blocks] ^= mask & xs;
    }

    if n & 7 != 0 {
        let start_idx = n & !7;
        eprintln!("s = {start_idx}");
        for i in 0..(n & 7) {
            let last_bs = bs.last().unwrap();
            if last_bs & (1 << i) != 0 {
                ys[start_idx + i] += x;
            }
        }
    }
}

pub fn bitmul_accumulate_naive(ys: &mut [GF2p8], x: GF2p8, bs: &[u8]) {
    let n = ys.len();
    debug_assert!(
        (bs.len() * 8 - 7 <= n) && (n <= bs.len() * 8),
        "ys.len() = {}, but bs.len() * 8 = {}",
        ys.len(),
        bs.len() * 8
    );
    let bs = GF2View::from_slice(bs);
    for (i, b) in bs.iter().take(n).enumerate() {
        if *b {
            ys[i] += x;
        }
    }
}

#[cfg(test)]
#[allow(non_snake_case)]
mod tests {
    use super::*;
    use crate::field::GF2Vector;
    use ff::Field;
    use rand::{thread_rng, Rng};

    #[test]
    fn test_hash() {
        let tau = 7;
        let n = 5;
        let m = 3;

        let r = GF2p8Vector::from(vec![
            GF2p8(0xf5),
            GF2p8(0xe8),
            GF2p8(0x13),
            GF2p8(0xb9),
            GF2p8(0xe3),
            GF2p8(0x48),
            GF2p8(0x86),
        ]);
        let v = GF2Vector {
            bits: [true, true, false, false, true].iter().collect(),
        };
        let M = GF2p8Matrix::from_shape_vec(
            (n, m),
            vec![
                GF2p8(0x25),
                GF2p8(0xad),
                GF2p8(0x90),
                GF2p8(0x9b),
                GF2p8(0xbe),
                GF2p8(0x57),
                GF2p8(0xec),
                GF2p8(0x7b),
                GF2p8(0xa0),
                GF2p8(0xf7),
                GF2p8(0xee),
                GF2p8(0x48),
                GF2p8(0x6f),
                GF2p8(0xc5),
                GF2p8(0xde),
            ],
        )
        .unwrap();
        let Rv = GF2p8Vector::from(vec![
            GF2p8(0x51),
            GF2p8(0xa3),
            GF2p8(0x5d),
            GF2p8(0x1e),
            GF2p8(0x12),
            GF2p8(0x51),
            GF2p8(0x49),
        ]);
        let RM = GF2p8Matrix::from_shape_vec(
            (tau, m),
            vec![
                GF2p8(0x8f),
                GF2p8(0xde),
                GF2p8(0xea),
                GF2p8(0x56),
                GF2p8(0x2f),
                GF2p8(0x64),
                GF2p8(0x67),
                GF2p8(0xd5),
                GF2p8(0xd5),
                GF2p8(0x10),
                GF2p8(0x15),
                GF2p8(0x1f),
                GF2p8(0xd4),
                GF2p8(0xbd),
                GF2p8(0x68),
                GF2p8(0x98),
                GF2p8(0x6b),
                GF2p8(0x10),
                GF2p8(0xb7),
                GF2p8(0x9c),
                GF2p8(0x1b),
            ],
        )
        .unwrap();
        assert_eq!(hash_matrix((&r).into(), (&M).into()), RM);
        assert_eq!(hash_bitvector((&r).into(), v.as_ref()), Rv);
        assert_eq!(
            hash_bitvector_and_matrix((&r).into(), v.as_ref(), (&M).into()),
            (Rv, RM)
        );
    }

    #[test]
    fn test_bitmul_accumulate() {
        let n = 16 + 8 + 3;
        let x = GF2p8::random(thread_rng());
        let mut ys_1: Vec<_> = (0..n).map(|_| GF2p8::random(thread_rng())).collect();
        let mut ys_2 = ys_1.clone();
        let mut bs = GF2Vector::with_capacity(n);
        bs.bits.resize(n, false);
        thread_rng().fill(bs.as_raw_mut_slice());

        bitmul_accumulate_naive(&mut ys_1, x, bs.as_raw_slice());
        bitmul_accumulate(&mut ys_2, x, bs.as_raw_slice());
        assert_eq!(ys_1, ys_2);
    }
}
