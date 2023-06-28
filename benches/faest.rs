use criterion::{
    black_box, criterion_group, criterion_main, measurement::WallTime, BenchmarkGroup, Criterion,
};
use faest::aes::Aes;
use faest::faest::{Prover, Verifier};
use faest::fiat_shamir::{SignatureVerifier, Signer};
use faest::gf2psmall::{GF2p10, GF2p11, GF2p7, GF2p8, GF2p9, SmallGF};
use faest::{keygen, FaestProver, FaestSignatureVerifier, FaestSigner, FaestVerifier};

pub fn bench_keygen(c: &mut Criterion) {
    let mut g = c.benchmark_group("keygen");
    for aes in [Aes::Aes128, Aes::Aes192, Aes::Aes256] {
        g.bench_function(format!("{aes:?}"), |b| {
            b.iter(|| black_box(keygen(aes)));
        });
    }
    g.finish();
}

pub fn bench_faest_interactive<F: SmallGF>(g: &mut BenchmarkGroup<WallTime>, aes: Aes) {
    g.bench_function("faest_interactive", |b| {
        let (secret_key, public_key) = keygen(aes);
        b.iter(|| {
            let mut prover = FaestProver::<F>::new(secret_key, public_key);
            let mut verifier = FaestVerifier::<F>::new(aes, public_key);

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

pub fn bench_faest_sign<F: SmallGF>(g: &mut BenchmarkGroup<WallTime>, aes: Aes) {
    g.bench_function("faest_sign", |b| {
        let message = "I see a ship in the harbor";
        let (secret_key, public_key) = keygen(aes);
        b.iter(|| {
            let signer = FaestSigner::<F>::new(secret_key, public_key);
            let signature = signer.sign(message.as_bytes());
            black_box(signature);
        });
    });
}

pub fn bench_faest_verify_signature<F: SmallGF>(g: &mut BenchmarkGroup<WallTime>, aes: Aes) {
    g.bench_function("faest_verify_signature", |b| {
        let message = "I see a ship in the harbor";
        let (secret_key, public_key) = keygen(aes);
        let signer = FaestSigner::<F>::new(secret_key, public_key);
        let signature = signer.sign(message.as_bytes());
        b.iter(|| {
            let verifier = FaestSignatureVerifier::<F>::new(aes, public_key);
            let result = verifier.verify(&signature, message.as_bytes());
            black_box(result);
        });
    });
}

pub fn bench_faest_instance<F: SmallGF>(c: &mut Criterion, aes: Aes) {
    let mut g = c.benchmark_group(format!("faest-{aes:?}-q={}", F::ORDER));
    bench_faest_interactive::<F>(&mut g, aes);
    bench_faest_sign::<F>(&mut g, aes);
    bench_faest_verify_signature::<F>(&mut g, aes);
    g.finish()
}

pub fn bench_faest(c: &mut Criterion) {
    for aes in [Aes::Aes128, Aes::Aes192, Aes::Aes256] {
        bench_faest_instance::<GF2p7>(c, aes);
        bench_faest_instance::<GF2p8>(c, aes);
        bench_faest_instance::<GF2p9>(c, aes);
        bench_faest_instance::<GF2p10>(c, aes);
        bench_faest_instance::<GF2p11>(c, aes);
    }
}

criterion_group!(benches, bench_keygen, bench_faest,);
criterion_main!(benches);
