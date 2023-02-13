mod aes;
mod arithmetic;
mod common;
pub mod faest;
pub mod field;
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
    homcom::HomCom128SenderFromVitH<voleith::VoleInTheHeadSenderFromVC<VC>>,
>;
pub type FaestVerifier = faest::FaestVerifierFromHC<
    homcom::HomCom128ReceiverFromVitH<voleith::VoleInTheHeadReceiverFromVC<VC>>,
>;
