use rug::Integer;

/// Chinese remainder theorem: the smallest non negative `x` with `x = residues[i] (mod
/// modulli[i])` for every `i`, `None` if the moduli are not coprime.
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
