use crate::field::GF2p8;
use bitvec::slice::BitSlice;
use ndarray::{Array1, Array2, ArrayView1, ArrayView2, Axis};

type GF2View = BitSlice<u8>;
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
    (hash_bitvector(r, v), hash_matrix(r, M))
}

#[cfg(test)]
#[allow(non_snake_case)]
mod tests {
    use super::*;
    use bitvec::vec::BitVec;
    type GF2Vector = BitVec<u8>;

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
        let v: GF2Vector = [true, true, false, false, true].iter().collect();
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
        assert_eq!(hash_bitvector((&r).into(), &v), Rv);
        assert_eq!(
            hash_bitvector_and_matrix((&r).into(), &v, (&M).into()),
            (Rv, RM)
        );
    }
}
