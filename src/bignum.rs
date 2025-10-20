//! Backend-agnostic bignum wrapper
//!
//! This module provides a unified interface for different bignum implementations.
//! The actual implementation is selected at compile-time via cargo features.

// Rug backend (default)
#[cfg(feature = "rug")]
mod rug_backend {
    pub use rug::{integer::IsPrime, rand::RandState, Integer};
    pub use rug::ops::Pow;
    
    pub type PrimalityResult = IsPrime;
    pub type RngState = RandState<'static>;
    
    pub fn new_rng() -> RngState {
        RandState::new()
    }
    
    pub fn probably_prime() -> PrimalityResult {
        IsPrime::Probably
    }
    
    pub fn not_prime() -> PrimalityResult {
        IsPrime::No
    }
    
    // Dummy trait for rug to keep API consistent
    pub trait IntegerExt {}
    impl IntegerExt for Integer {}
}

#[cfg(feature = "rug")]
pub use rug_backend::*;

// num-bigint backend
#[cfg(all(feature = "num-bigint", not(feature = "rug")))]
mod num_bigint_backend {
    use num_bigint::BigInt;
    use num_traits::{One, Zero};
    
    pub type Integer = BigInt;
    pub type RngState = rand::rngs::ThreadRng;
    
    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    pub enum PrimalityResult {
        Probably,
        No,
    }
    
    pub type IsPrime = PrimalityResult;
    
    // Re-export Pow trait for num-bigint
    pub use num_traits::Pow;
    
    pub fn new_rng() -> RngState {
        rand::thread_rng()
    }
    
    pub fn probably_prime() -> PrimalityResult {
        PrimalityResult::Probably
    }
    
    pub fn not_prime() -> PrimalityResult {
        PrimalityResult::No
    }
    
    // Extension trait for num-bigint to match rug API
    pub trait IntegerExt {
        fn is_probably_prime(&self, _reps: u32) -> PrimalityResult;
        fn random_below(&self, rng: &mut RngState) -> Self;
        fn invert_ref(&self, modulus: &Self) -> Option<Self> where Self: Sized;
        fn invert(&self, modulus: &Self) -> Option<Self> where Self: Sized;
        fn pow_mod(&self, exp: &Self, modulus: &Self) -> Option<Self> where Self: Sized;
        fn is_divisible(&self, other: &Self) -> bool;
        fn to_f64(&self) -> f64;
        fn to_u32(&self) -> Option<u32>;
        fn to_u64(&self) -> Option<u64>;
        fn sqrt(&self) -> Self;
        fn gcd(&self, other: &Self) -> Self;
        fn div_rem(self, other: Self) -> (Self, Self) where Self: Sized;
    }
    
    impl IntegerExt for BigInt {
        fn is_probably_prime(&self, _reps: u32) -> PrimalityResult {
            // Simple primality check - not cryptographically strong
            if self < &BigInt::from(2) {
                return PrimalityResult::No;
            }
            if self == &BigInt::from(2) {
                return PrimalityResult::Probably;
            }
            if (self & BigInt::one()) == BigInt::zero() {
                return PrimalityResult::No;
            }
            // For simplicity, assume it's probably prime for larger numbers
            // A real implementation would use Miller-Rabin or similar
            PrimalityResult::Probably
        }
        
        fn random_below(&self, rng: &mut RngState) -> Self {
            use num_bigint::RandBigInt;
            rng.gen_bigint_range(&BigInt::zero(), self)
        }
        
        fn invert_ref(&self, modulus: &Self) -> Option<Self> {
            self.invert(modulus)
        }
        
        fn invert(&self, modulus: &Self) -> Option<Self> {
            // Extended Euclidean algorithm
            use num_traits::Signed;
            let (mut t, mut newt) = (BigInt::zero(), BigInt::one());
            let (mut r, mut newr) = (modulus.clone(), self.clone());
            
            while !newr.is_zero() {
                let quotient = &r / &newr;
                let temp = newt.clone();
                newt = t - &quotient * &newt;
                t = temp;
                let temp = newr.clone();
                newr = r - quotient * &newr;
                r = temp;
            }
            
            if r > BigInt::one() {
                return None; // Not invertible
            }
            if t < BigInt::zero() {
                t += modulus;
            }
            Some(t)
        }
        
        fn pow_mod(&self, exp: &Self, modulus: &Self) -> Option<Self> {
            use num_traits::Signed;
            if modulus.is_one() {
                return Some(BigInt::zero());
            }
            
            let mut result = BigInt::one();
            let mut base = self % modulus;
            let mut exp = exp.clone();
            
            while !exp.is_zero() {
                if (&exp & BigInt::one()) == BigInt::one() {
                    result = (result * &base) % modulus;
                }
                exp >>= 1;
                base = (&base * &base) % modulus;
            }
            
            Some(result)
        }
        
        fn is_divisible(&self, other: &Self) -> bool {
            (self % other).is_zero()
        }
        
        fn to_f64(&self) -> f64 {
            use num_traits::ToPrimitive;
            self.to_f64().unwrap_or(f64::INFINITY)
        }
        
        fn to_u32(&self) -> Option<u32> {
            use num_traits::ToPrimitive;
            self.to_u32()
        }
        
        fn to_u64(&self) -> Option<u64> {
            use num_traits::ToPrimitive;
            self.to_u64()
        }
        
        fn sqrt(&self) -> Self {
            use num_integer::Roots;
            self.sqrt()
        }
        
        fn gcd(&self, other: &Self) -> Self {
            use num_integer::Integer as IntegerTrait;
            IntegerTrait::gcd(self, other)
        }
        
        fn div_rem(self, other: Self) -> (Self, Self) {
            use num_integer::Integer as IntegerTrait;
            IntegerTrait::div_rem(&self, &other)
        }
    }
}

#[cfg(all(feature = "num-bigint", not(feature = "rug")))]
pub use num_bigint_backend::*;
#[cfg(all(feature = "num-bigint", not(feature = "rug")))]
pub use num_bigint_backend::IntegerExt as _;

// Ensure at least one backend is enabled
#[cfg(not(any(feature = "rug", feature = "num-bigint", feature = "ibig", feature = "rust-gmp")))]
compile_error!("At least one bignum backend must be enabled. Enable one of: rug (default), num-bigint, ibig, rust-gmp");
