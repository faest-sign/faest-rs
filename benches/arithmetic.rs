use criterion::{criterion_group, criterion_main, Criterion};
use ff::Field;
use homcomzk::arithmetic::{bitmul_accumulate, bitmul_accumulate_naive};
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

criterion_group!(benches, bench_bitmul_accumulate);
criterion_main!(benches);
