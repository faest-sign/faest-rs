use criterion::{black_box, criterion_group, criterion_main, Criterion};
use faest::arithmetic::{
    bit_xor_assign, bit_xor_assign_naive, bitmul_accumulate_naive, bitmul_accumulate_u16,
    bitmul_accumulate_u8, clmul, clmul_u64, clmul_u8, clmul_u8_x86,
};
use faest::gf2::GF2Vector;
use faest::gf2psmall::{GF2p8, GF2p9};
use ff::Field;
use rand::{thread_rng, Rng};

pub fn bench_clmul(c: &mut Criterion) {
    let mut g = c.benchmark_group("clmul");

    g.bench_function("8/naive", |b| {
        let x = thread_rng().gen();
        let y = thread_rng().gen();
        b.iter(|| black_box(clmul_u8(x, y)));
    });
    g.bench_function("8/clever", |b| {
        let x = thread_rng().gen();
        let y = thread_rng().gen();
        b.iter(|| black_box(clmul_u8_x86(x, y)));
    });
    g.bench_function("64/naive", |b| {
        let x = thread_rng().gen();
        let y = thread_rng().gen();
        b.iter(|| black_box(clmul_u64(x, y)));
    });
    g.bench_function("64/clever", |b| {
        let x = thread_rng().gen();
        let y = thread_rng().gen();
        b.iter(|| black_box(clmul(x, y)));
    });
}

pub fn bench_bitmul_accumulate(c: &mut Criterion) {
    let n = 2000;
    let x_8 = GF2p8::random(thread_rng());
    let mut ys_8: Vec<_> = (0..n).map(|_| GF2p8::random(thread_rng())).collect();
    let x_16 = GF2p9::random(thread_rng());
    let mut ys_16: Vec<_> = (0..n).map(|_| GF2p9::random(thread_rng())).collect();
    let mut bs = GF2Vector::with_capacity(n);
    bs.bits.resize(n, false);
    thread_rng().fill(bs.as_raw_mut_slice());

    let mut g = c.benchmark_group("bitmul_accumulate");

    g.bench_function("clever-u8", |b| {
        b.iter(|| unsafe { bitmul_accumulate_u8(&mut ys_8, x_8, bs.as_raw_slice()) });
    });
    g.bench_function("clever-u16", |b| {
        b.iter(|| unsafe { bitmul_accumulate_u16(&mut ys_16, x_16, bs.as_raw_slice()) });
    });
    g.bench_function("naive-u8", |b| {
        b.iter(|| bitmul_accumulate_naive(&mut ys_8, x_8, bs.as_raw_slice()));
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

criterion_group!(
    benches,
    bench_clmul,
    bench_bitmul_accumulate,
    bench_bit_xor_assign
);
criterion_main!(benches);
