//! Deterministic inputs shared by the benchmarks.

#![allow(dead_code)]

use rug::{integer::IsPrime, ops::Pow, rand::RandState, Integer};
use std::{collections::HashMap, str::FromStr};

/// Seed used for every generated input, so all runs measure exactly the same work.
pub const SEED: u64 = 0xD15C;

/// Number of targets solved per instance: the randomized algorithms (Pollard's rho, index
/// calculus) always draw the same random numbers, so the luck is averaged over several targets
/// instead of several seeds.
pub const TARGETS: usize = 3;

/// A discrete logarithm problem `b^x = a (mod n)` for each `a` of `targets`, where `b` has
/// order `order` modulo `n`.
#[derive(Clone, Debug)]
pub struct Instance {
    pub n: Integer,
    pub b: Integer,
    pub order: Integer,
    pub targets: Vec<Integer>,
}

impl Instance {
    /// Instance with the given targets and the order of `b` computed.
    pub fn new(n: &str, b: &str, targets: &[&str]) -> Self {
        let n = int(n);
        let b = int(b);
        let order = discrete_logarithm::n_order(&b, &n).unwrap();
        let targets = targets.iter().map(|a| int(a)).collect();
        Self {
            n,
            b,
            order,
            targets,
        }
    }

    /// Instance with a known order and a single target.
    pub fn with_order(n: &str, a: &str, b: &str, order: &str) -> Self {
        Self {
            n: int(n),
            b: int(b),
            order: int(order),
            targets: vec![int(a)],
        }
    }

    /// `n = 2q + 1` safe prime of `bits` bits and `b = 4`, of prime order `q`: every algorithm
    /// for prime orders applies. The targets are random powers of `b`.
    pub fn safe_prime(bits: u32, seed: u64) -> Self {
        let mut q = prime_bits(bits - 1, seed);
        let n = loop {
            let n = Integer::from(&q * 2u32) + 1u32;
            if n.is_probably_prime(30) != IsPrime::No {
                break n;
            }
            q.next_prime_mut();
        };
        assert_eq!(n.significant_bits(), bits);
        let b = Integer::from(4);

        let mut rand = RandState::new();
        rand.seed(&Integer::from(seed));
        let targets = (0..TARGETS)
            .map(|_| {
                let x = Integer::from(q.random_below_ref(&mut rand));
                b.clone().pow_mod(&x, &n).unwrap()
            })
            .collect();
        Self {
            n,
            b,
            order: q,
            targets,
        }
    }

    /// Checks that each `x` of `logs` is a discrete logarithm of its target.
    pub fn check(&self, logs: &[Integer]) {
        assert_eq!(logs.len(), self.targets.len());
        for (a, x) in self.targets.iter().zip(logs) {
            assert_eq!(
                &self.b.clone().pow_mod(x, &self.n).unwrap(),
                a,
                "wrong logarithm {x} of {a} in base {} modulo {}",
                self.b,
                self.n
            );
        }
    }
}

pub fn int(s: &str) -> Integer {
    Integer::from_str(s).unwrap()
}

/// Random prime of exactly `bits` bits, derived from `seed`.
pub fn prime_bits(bits: u32, seed: u64) -> Integer {
    let mut rand = RandState::new();
    rand.seed(&Integer::from(seed));
    loop {
        let mut x = Integer::from(Integer::random_bits(bits, &mut rand));
        x.set_bit(bits - 1, true);
        let p = x.next_prime();
        if p.significant_bits() == bits {
            return p;
        }
    }
}

/// Checks that `factors` is the prime factorization of `n`.
pub fn check_factors(n: &Integer, factors: &HashMap<Integer, usize>) {
    let product = factors.iter().fold(Integer::from(1), |acc, (p, e)| {
        acc * p.clone().pow(*e as u32)
    });
    assert_eq!(&product, n, "wrong factorization of {n}: {factors:?}");
    for p in factors.keys() {
        assert_ne!(
            p.is_probably_prime(30),
            IsPrime::No,
            "composite factor {p} of {n}"
        );
    }
}

/// The 108-digit instance of the tests (smooth group order: Pohlig-Hellman).
pub const DIGITS_108: (&str, &str, &str) = (
    "22708823198678103974314518195029102158525052496759285596453269189798311427475159776411276642277139650833937",
    "17463946429475485293747680247507700244427944625055089103624311227422110546803452417458985046168310373075327",
    "123456",
);
