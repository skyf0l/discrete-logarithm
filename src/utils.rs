use std::collections::HashMap;

use primal::Primes;
use rug::Integer;

pub fn fast_factor(n: &Integer) -> HashMap<Integer, usize> {
    let mut factors: HashMap<Integer, usize> = HashMap::new();
    let mut n: Integer = n.clone();
    for prime in Primes::all().take(1_000_000) {
        let prime = Integer::from(prime);
        if n.clone().div_rem(prime.clone()).1 == 0 {
            // factors.insert(prime.clone(), 1);
            while n.clone().div_rem(prime.clone()).1 == 0 {
                n /= &prime;
                *factors.entry(prime.clone()).or_insert(0) += 1;
            }
        }
    }

    if n != 1 {
        *factors.entry(n).or_insert(0) += 1;
    }

    factors
}

pub fn crt(residues: &[Integer], modulli: &[Integer]) -> Option<Integer> {
    let prod = modulli.iter().product::<Integer>();
    let mut sum = Integer::ZERO;

    for (residue, modulus) in residues.iter().zip(modulli) {
        let p = prod.clone() / modulus;
        sum += residue * Integer::from(p.invert_ref(modulus)?) * p
    }

    Some(sum % prod)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chinese_remainder_theorem() {
        assert_eq!(
            crt(
                &[3.into(), 5.into(), 7.into()],
                &[2.into(), 3.into(), 1.into()]
            ),
            Some(Integer::from(5))
        );
        assert_eq!(
            crt(
                &[1.into(), 4.into(), 6.into()],
                &[3.into(), 5.into(), 7.into()]
            ),
            Some(Integer::from(34))
        );
        assert_eq!(
            crt(
                &[1.into(), 4.into(), 6.into()],
                &[1.into(), 2.into(), 0.into()]
            ),
            None
        );
        assert_eq!(
            crt(
                &[2.into(), 5.into(), 7.into()],
                &[6.into(), 9.into(), 15.into()]
            ),
            None
        );
    }
}
