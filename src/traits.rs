use core::default::Default;
use core::fmt::Debug;
use digest::XofReader;
use rand::Rng;

pub trait Block: Debug + Copy + Default + AsRef<[u8]> + AsMut<[u8]> {
    const BIT_LENGTH: usize;
    fn random(rng: &mut impl Rng) -> Self;
}

pub trait LengthDoubling {
    type Block: Block;
    fn expand(seed: Self::Block) -> (Self::Block, Self::Block);
}

pub trait LeafExpander {
    type Com;
    type XofR: XofReader;
    fn expand_leaf(leaf: impl AsRef<[u8]>) -> (Self::Com, Self::XofR);
}
