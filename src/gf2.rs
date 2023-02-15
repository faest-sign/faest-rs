use bitvec::{slice::BitSlice, vec::BitVec};

pub type GF2View = BitSlice<u8>;
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct GF2Vector {
    pub bits: BitVec<u8>,
}

impl GF2Vector {
    pub fn new() -> Self {
        Self {
            bits: BitVec::new(),
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            bits: BitVec::with_capacity(capacity),
        }
    }

    pub fn from_vec(vec: Vec<u8>) -> Self {
        Self {
            bits: BitVec::from_vec(vec),
        }
    }

    pub fn len(&self) -> usize {
        self.bits.len()
    }

    pub fn is_empty(&self) -> bool {
        self.bits.is_empty()
    }

    pub fn resize(&mut self, size: usize, value: bool) {
        self.bits.resize(size, value)
    }

    pub fn as_raw_slice(&self) -> &[u8] {
        self.bits.as_raw_slice()
    }

    pub fn as_raw_mut_slice(&mut self) -> &mut [u8] {
        self.bits.as_raw_mut_slice()
    }
}

impl From<BitVec<u8>> for GF2Vector {
    fn from(bits: BitVec<u8>) -> Self {
        Self { bits }
    }
}

impl From<GF2Vector> for BitVec<u8> {
    fn from(vec: GF2Vector) -> Self {
        vec.bits
    }
}

impl AsRef<BitSlice<u8>> for GF2Vector {
    fn as_ref(&self) -> &BitSlice<u8> {
        &self.bits
    }
}

impl AsMut<BitSlice<u8>> for GF2Vector {
    fn as_mut(&mut self) -> &mut BitSlice<u8> {
        &mut self.bits
    }
}

impl bincode::Encode for GF2Vector {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> core::result::Result<(), bincode::error::EncodeError> {
        debug_assert_eq!(
            {
                let mut cpy = self.bits.clone();
                cpy.force_align();
                cpy.into_vec()
            },
            {
                let cpy = self.bits.clone();
                cpy.into_vec()
            }
        );
        bincode::Encode::encode(&self.bits.len(), encoder)?;
        bincode::Encode::encode(&self.bits.as_raw_slice(), encoder)?;
        Ok(())
    }
}
