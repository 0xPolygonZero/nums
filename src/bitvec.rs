use crate::util::ceil_div_usize;
use alloc::vec;
use alloc::vec::Vec;
use core::convert::identity;
use core::fmt::{Debug, Formatter};
use core::ops::BitXorAssign;

#[derive(Clone, Eq, PartialEq, Hash)]
pub struct BitVec {
    bytes: Vec<u8>,
    len: usize,
}

impl BitVec {
    #[must_use]
    pub fn new(len: usize) -> Self {
        let bytes = vec![0; ceil_div_usize(len, 8)];
        Self { bytes, len }
    }

    #[inline]
    pub fn set(&mut self, i: usize, v: bool) {
        assert!(i < self.len);
        self.bytes[i / 8] |= (v as u8) << (i % 8);
    }

    #[must_use]
    #[inline]
    pub fn get(&self, i: usize) -> bool {
        assert!(i < self.len);
        let byte = self.bytes[i / 8];
        byte & (1 << (i % 8)) != 0
    }

    #[must_use]
    #[inline]
    pub fn count_zeros(&self) -> usize {
        self.len() - self.count_ones()
    }

    #[must_use]
    #[inline]
    pub fn count_ones(&self) -> usize {
        self.bytes.iter().map(|byte| byte.count_ones()).sum::<u32>() as usize
    }

    #[must_use]
    #[inline]
    pub fn leading_zeros(&self) -> usize {
        self.iter().position(identity).unwrap_or(self.len)
    }

    #[must_use]
    #[inline]
    pub fn len(&self) -> usize {
        self.len
    }

    #[must_use]
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    #[inline]
    fn resize(&mut self, len: usize) {
        // Shrinking might leave dirty bits, which we'd have to clear, or account for in Hash etc.
        assert!(len >= self.len, "Shrinking not supported for now");

        let num_bytes = ceil_div_usize(len, 8);
        self.bytes.resize(num_bytes, 0);
        self.len = len;
    }

    #[must_use]
    pub fn iter(&self) -> BitVecIter {
        BitVecIter {
            bit_vec: self,
            position: 0,
        }
    }

    /// Rotate `k` bits to the left, i.e. in the direction of smaller indices.
    pub fn rotl(&self, k: usize) -> Self {
        let l = self.len;
        Self::from_iter((0..l).map(|i| self.get((i + k) % l)))
    }

    /// Rotate `k` bits to the right, i.e. in the direction of greater indices.
    pub fn rotr(&self, k: usize) -> Self {
        let l = self.len;
        let k = k % l;
        Self::from_iter((0..l).map(|i| self.get((l + i - k) % l)))
    }
}

impl<const N: usize> From<[bool; N]> for BitVec {
    fn from(value: [bool; N]) -> Self {
        let mut result = BitVec::new(N);
        for (i, bit) in value.into_iter().enumerate() {
            result.set(i, bit);
        }
        result
    }
}

impl BitXorAssign<&Self> for BitVec {
    fn bitxor_assign(&mut self, rhs: &Self) {
        assert_eq!(self.len(), rhs.len());
        self.bytes
            .iter_mut()
            .zip(rhs.bytes.iter())
            .for_each(|(self_byte, rhs_byte)| *self_byte ^= rhs_byte);
    }
}

pub struct BitVecIter<'a> {
    bit_vec: &'a BitVec,
    position: usize,
}

impl<'a> Iterator for BitVecIter<'a> {
    type Item = bool;

    fn next(&mut self) -> Option<bool> {
        if self.position == self.bit_vec.len {
            None
        } else {
            let b = self.bit_vec.get(self.position);
            self.position += 1;
            Some(b)
        }
    }
}

impl FromIterator<bool> for BitVec {
    fn from_iter<T: IntoIterator<Item = bool>>(iter: T) -> Self {
        let mut res = Self::new(0);
        for b in iter {
            let i = res.len;
            res.resize(i + 1);
            res.set(i, b);
        }
        res
    }
}

impl Debug for BitVec {
    fn fmt(&self, f: &mut Formatter<'_>) -> core::fmt::Result {
        for i in 0..self.len {
            write!(f, "{}", self.get(i) as u8)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::bitvec::BitVec;

    #[test]
    fn test_count_zeros() {
        assert_eq!(BitVec::from([false, true, false, true]).count_zeros(), 2);
        assert_eq!(
            BitVec::from([false, true, false, true, false, true, false, true]).count_zeros(),
            4
        );
    }

    #[test]
    fn test_count_ones() {
        assert_eq!(BitVec::from([false, true, false, true]).count_ones(), 2);
        assert_eq!(
            BitVec::from([false, true, false, true, false, true, false, true]).count_ones(),
            4
        );
    }

    #[test]
    fn test_rotl() {
        assert_eq!(BitVec::from([true, false, false, false, true]).rotl(1), BitVec::from([false, false, false, true, true]));
    }

    #[test]
    fn test_rotr() {
        assert_eq!(BitVec::from([true, false, false, false, true]).rotr(1), BitVec::from([true, true, false, false, false]));
    }
}
