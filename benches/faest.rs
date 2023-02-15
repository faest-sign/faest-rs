use criterion::{
    black_box, criterion_group, criterion_main, measurement::WallTime, BenchmarkGroup, Criterion,
};
use homcomzk::faest::{Prover, Verifier};
use homcomzk::fiat_shamir::{SignatureVerifier, Signer};
use homcomzk::gf2psmall::{GF2p10, GF2p11, GF2p7, GF2p8, GF2p9, SmallGF};
use homcomzk::{keygen, FaestProver, FaestSignatureVerifier, FaestSigner, FaestVerifier};

pub fn bench_keygen(c: &mut Criterion) {
    c.bench_function("keygen", |b| {
        b.iter(|| black_box(keygen()));
    });
}

pub fn bench_faest_interactive<F: SmallGF>(g: &mut BenchmarkGroup<WallTime>) {
    g.bench_function("faest_interactive", |b| {
        let (secret_key, public_key) = keygen();
        b.iter(|| {
            let mut prover = FaestProver::<F>::new(secret_key, public_key);
            let mut verifier = FaestVerifier::<F>::new(public_key);

            let commitment = prover.commit();
            let challenge = verifier.commit_send_challenge(commitment);
            let response = prover.commit_prove_consistency(challenge);
            verifier.commit_receive_response(response);

            let challenge = verifier.send_challenge();
            let proof = prover.prove(challenge);
            let choice = verifier.choose(proof);
            let decommitment = prover.transfer(choice);
            let result = verifier.verify(decommitment);
            black_box(result);
        });
    });
}

pub fn bench_faest_sign<F: SmallGF>(g: &mut BenchmarkGroup<WallTime>) {
    g.bench_function("faest_sign", |b| {
        let message = "I see a ship in the harbor";
        let (secret_key, public_key) = keygen();
        b.iter(|| {
            let signer = FaestSigner::<F>::new(secret_key, public_key);
            let signature = signer.sign(message.as_bytes());
            black_box(signature);
        });
    });
}

pub fn bench_faest_verify_signature<F: SmallGF>(g: &mut BenchmarkGroup<WallTime>) {
    g.bench_function("faest_verify_signature", |b| {
        let message = "I see a ship in the harbor";
        let (secret_key, public_key) = keygen();
        let signer = FaestSigner::<F>::new(secret_key, public_key);
        let signature = signer.sign(message.as_bytes());
        b.iter(|| {
            let verifier = FaestSignatureVerifier::<F>::new(public_key);
            let result = verifier.verify(&signature, message.as_bytes());
            black_box(result);
        });
    });
}

pub fn bench_faest_instance<F: SmallGF>(c: &mut Criterion) {
    let mut g = c.benchmark_group(format!("faest_q={}", F::ORDER));
    bench_faest_interactive::<F>(&mut g);
    bench_faest_sign::<F>(&mut g);
    bench_faest_verify_signature::<F>(&mut g);
    g.finish()
}

pub fn bench_faest(c: &mut Criterion) {
    bench_faest_instance::<GF2p7>(c);
    bench_faest_instance::<GF2p8>(c);
    bench_faest_instance::<GF2p9>(c);
    bench_faest_instance::<GF2p10>(c);
    bench_faest_instance::<GF2p11>(c);
}

criterion_group!(benches, bench_keygen, bench_faest,);
criterion_main!(benches);
