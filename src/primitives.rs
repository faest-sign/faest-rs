use crate::common::Block128;
use crate::traits::{LeafExpander, LengthDoubling};
use aes::cipher::crypto_common::Block;
use aes::cipher::{BlockEncrypt, KeyInit};
use aes::Aes128;
use core::marker::PhantomData;

pub struct Aes128CtrLdPrg {}

impl LengthDoubling for Aes128CtrLdPrg {
    type Block = Block128;

    fn expand(seed: Self::Block) -> (Self::Block, Self::Block) {
        let aes = Aes128::new_from_slice(seed.as_ref())
            .expect("does not fail since seed has the right size");
        let mut b0 = Block128([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]);
        let mut b1 = Block128([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]);
        aes.encrypt_block(Block::<Aes128>::from_mut_slice(b0.as_mut()));
        aes.encrypt_block(Block::<Aes128>::from_mut_slice(b1.as_mut()));
        (b0, b1)
    }
}

pub struct Blake3LE<Com> {
    _phantom_com: PhantomData<Com>,
}

impl<Com> LeafExpander for Blake3LE<Com>
where
    Com: Default + AsMut<[u8]>,
{
    type Com = Com;
    type XofR = blake3::OutputReader;
    fn expand_leaf(leaf: impl AsRef<[u8]>) -> (Self::Com, Self::XofR) {
        let mut hasher = blake3::Hasher::new();
        hasher.update(leaf.as_ref());
        let mut xof = hasher.finalize_xof();
        let mut com: Self::Com = Default::default();
        xof.fill(com.as_mut());
        (com, xof)
    }
}
