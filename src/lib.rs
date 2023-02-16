pub mod aes;
pub mod arithmetic;
mod common;
pub mod faest;
pub mod fiat_shamir;
pub mod field;
#[macro_use]
mod field_helpers;
pub mod gf2;
pub mod gf2p128;
pub mod gf2psmall;
pub mod homcom;
mod primitives;
mod traits;
pub mod veccom;
pub mod voleith;

pub use faest::{keygen, PublicKey, SecretKey};

type VC = veccom::GgmVecCom<
    common::Block128,
    primitives::Aes128CtrLdPrg,
    blake3::Hasher,
    primitives::Blake3LE<common::Block128>,
>;
pub type FaestProver<F> = faest::FaestProverFromHC<
    homcom::HomCom128SenderFromVitH<voleith::VoleInTheHeadSenderFromVC<F, VC>>,
>;
pub type FaestVerifier<F> = faest::FaestVerifierFromHC<
    homcom::HomCom128ReceiverFromVitH<voleith::VoleInTheHeadReceiverFromVC<F, VC>>,
>;

pub type FaestSignature<F> = fiat_shamir::FsSignature<FaestProver<F>>;
pub type FaestSigner<F> = fiat_shamir::FsSigner<FaestProver<F>, FaestVerifier<F>>;
pub type FaestSignatureVerifier<F> = fiat_shamir::FsVerifier<FaestProver<F>, FaestVerifier<F>>;
