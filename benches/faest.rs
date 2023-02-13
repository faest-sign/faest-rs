use criterion::{black_box, criterion_group, criterion_main, Criterion};
use homcomzk::faest::{Prover, Verifier};
use homcomzk::{keygen, FaestProver, FaestVerifier};

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

criterion_group!(benches, bench_keygen, bench_faest_interactive);
criterion_main!(benches);
