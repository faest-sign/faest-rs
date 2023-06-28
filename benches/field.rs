use criterion::{black_box, criterion_group, criterion_main, Criterion};
use faest::field::InnerProduct;
use faest::gf2p128::GF2p128Fast;
// use faest::gf2p128::GF2p128Naive;
use faest::gf2psmall::{GF2p10, GF2p11, GF2p7, GF2p8, GF2p9};
use ff::Field;
use rand::thread_rng;

fn bench_field<F>(name: &str, c: &mut Criterion)
where
    F: Field + for<'a> InnerProduct<&'a F>,
{
    let mut g = c.benchmark_group(name);
    let n = 1000;

    g.bench_function("add", |b| {
        let x = F::random(thread_rng());
        let y = F::random(thread_rng());
        b.iter(|| black_box(x + y));
    });
    g.bench_function("mul", |b| {
        let x = F::random(thread_rng());
        let y = F::random(thread_rng());
        b.iter(|| black_box(x * y));
    });
    g.bench_function("invert", |b| {
        let x = F::random(thread_rng());
        b.iter(|| black_box(x.invert()));
    });
    g.bench_function("sum", |b| {
        let x: Vec<_> = (0..n).map(|_| F::random(thread_rng())).collect();
        b.iter(|| {
            let z: F = x.iter().sum();
            black_box(z)
        });
    });
    g.bench_function("product", |b| {
        let x: Vec<_> = (0..n).map(|_| F::random(thread_rng())).collect();
        b.iter(|| {
            let z: F = x.iter().product();
            black_box(z)
        });
    });
    g.bench_function("inner_product", |b| {
        let x: Vec<_> = (0..n).map(|_| F::random(thread_rng())).collect();
        let y: Vec<_> = (0..n).map(|_| F::random(thread_rng())).collect();
        b.iter(|| {
            let z: F = InnerProduct::<&F>::inner_product(x.iter(), y.iter());
            black_box(z)
        });
    });
    g.bench_function("elementwise_add", |b| {
        let x: Vec<_> = (0..n).map(|_| F::random(thread_rng())).collect();
        let y: Vec<_> = (0..n).map(|_| F::random(thread_rng())).collect();
        b.iter(|| {
            let z: Vec<_> = x.iter().zip(y.iter()).map(|(&a, &b)| a + b).collect();
            black_box(z)
        });
    });
    g.bench_function("elementwise_mul", |b| {
        let x: Vec<_> = (0..n).map(|_| F::random(thread_rng())).collect();
        let y: Vec<_> = (0..n).map(|_| F::random(thread_rng())).collect();
        b.iter(|| {
            let z: Vec<_> = x.iter().zip(y.iter()).map(|(&a, &b)| a * b).collect();
            black_box(z)
        });
    });
}

fn bench_gf2p7(c: &mut Criterion) {
    bench_field::<GF2p7>("gf2p7", c)
}

fn bench_gf2p8(c: &mut Criterion) {
    bench_field::<GF2p8>("gf2p8", c)
}

fn bench_gf2p9(c: &mut Criterion) {
    bench_field::<GF2p9>("gf2p9", c)
}

fn bench_gf2p10(c: &mut Criterion) {
    bench_field::<GF2p10>("gf2p10", c)
}

fn bench_gf2p11(c: &mut Criterion) {
    bench_field::<GF2p11>("gf2p11", c)
}

// fn bench_gf2p128(c: &mut Criterion) {
//     bench_field::<GF2p128Naive>("gf2p128-naive", c)
// }

fn bench_gf2p128_fast(c: &mut Criterion) {
    bench_field::<GF2p128Fast>("gf2p128-fast", c)
}

criterion_group!(
    benches,
    bench_gf2p7,
    bench_gf2p8,
    bench_gf2p9,
    bench_gf2p10,
    bench_gf2p11,
    // bench_gf2p128,
    bench_gf2p128_fast
);
criterion_main!(benches);
