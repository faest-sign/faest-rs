use criterion::{criterion_group, criterion_main, Criterion};
use ff::Field;
use homcomzk::arithmetic::{
    bit_xor_assign, bit_xor_assign_naive, bitmul_accumulate, bitmul_accumulate_naive,
};
use homcomzk::field::{GF2Vector, GF2p8};
use rand::{thread_rng, Rng};

pub fn bench_bitmul_accumulate(c: &mut Criterion) {
    let n = 2000;
    let x = GF2p8::random(thread_rng());
    let mut ys: Vec<_> = (0..n).map(|_| GF2p8::random(thread_rng())).collect();
    let mut bs = GF2Vector::with_capacity(n);
    bs.bits.resize(n, false);
    thread_rng().fill(bs.as_raw_mut_slice());

    let mut g = c.benchmark_group("bitmul_accumulate");

    g.bench_function("clever", |b| {
        b.iter(|| bitmul_accumulate(&mut ys, x, bs.as_raw_slice()));
    });
    g.bench_function("naive", |b| {
        b.iter(|| bitmul_accumulate_naive(&mut ys, x, bs.as_raw_slice()));
    });
}

pub fn bench_bit_xor_assign(c: &mut Criterion) {
    let n = 10000;
    let mut xs = GF2Vector::with_capacity(n);
    xs.bits.resize(n, false);
    thread_rng().fill(xs.as_raw_mut_slice());
    let mut ys = GF2Vector::with_capacity(n);
    ys.bits.resize(n, false);
    thread_rng().fill(ys.as_raw_mut_slice());

    let mut g = c.benchmark_group("bit_xor_assign");

    g.bench_function("clever", |b| {
        b.iter(|| bit_xor_assign(&mut ys, &xs));
    });
    g.bench_function("naive", |b| {
        b.iter(|| bit_xor_assign_naive(&mut ys, &xs));
    });
}

criterion_group!(benches, bench_bitmul_accumulate, bench_bit_xor_assign);
criterion_main!(benches);
