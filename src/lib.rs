mod aes;
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
pub type FaestProver = faest::FaestProverFromHC<
    homcom::HomCom128SenderFromVitH<voleith::VoleInTheHeadSenderFromVC<gf2psmall::GF2p8, VC>>,
>;
pub type FaestVerifier = faest::FaestVerifierFromHC<
    homcom::HomCom128ReceiverFromVitH<voleith::VoleInTheHeadReceiverFromVC<gf2psmall::GF2p8, VC>>,
>;

pub type FaestSignature = fiat_shamir::FsSignature<FaestProver>;
pub type FaestSigner = fiat_shamir::FsSigner<FaestProver, FaestVerifier>;
pub type FaestSignatureVerifier = fiat_shamir::FsVerifier<FaestProver, FaestVerifier>;
