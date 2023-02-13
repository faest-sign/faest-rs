use criterion::{black_box, criterion_group, criterion_main, Criterion};
use homcomzk::faest::{Prover, Verifier};
use homcomzk::fiat_shamir::{SignatureVerifier, Signer};
use homcomzk::{keygen, FaestProver, FaestSignatureVerifier, FaestSigner, FaestVerifier};

pub fn bench_keygen(c: &mut Criterion) {
    c.bench_function("keygen", |b| {
        b.iter(|| black_box(keygen()));
    });
}

pub fn bench_faest_interactive(c: &mut Criterion) {
    c.bench_function("faest_interactive", |b| {
        let (secret_key, public_key) = keygen();
        b.iter(|| {
            let mut prover = FaestProver::new(secret_key, public_key);
            let mut verifier = FaestVerifier::new(public_key);

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

pub fn bench_faest_sign(c: &mut Criterion) {
    c.bench_function("faest_sign", |b| {
        let message = "I see a ship in the harbor";
        let (secret_key, public_key) = keygen();
        b.iter(|| {
            let signer = FaestSigner::new(secret_key, public_key);
            let signature = signer.sign(message.as_bytes());
            black_box(signature);
        });
    });
}

pub fn bench_faest_verify_signature(c: &mut Criterion) {
    c.bench_function("faest_verify_signature", |b| {
        let message = "I see a ship in the harbor";
        let (secret_key, public_key) = keygen();
        let signer = FaestSigner::new(secret_key, public_key);
        let signature = signer.sign(message.as_bytes());
        b.iter(|| {
            let verifier = FaestSignatureVerifier::new(public_key);
            let result = verifier.verify(&signature, message.as_bytes());
            black_box(result);
        });
    });
}

criterion_group!(
    benches,
    bench_keygen,
    bench_faest_interactive,
    bench_faest_sign,
    bench_faest_verify_signature
);
criterion_main!(benches);
