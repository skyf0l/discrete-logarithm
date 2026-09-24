use primal::Primes;
use rug::{rand::RandState, Integer};

use crate::Error;

/// Check if a number can be factored using the given factor base.
/// Returns the exponents vector if smooth, None otherwise.
pub fn is_smooth(mut n: Integer, factorbase: &[u32]) -> Option<Vec<u32>> {
    let mut factors = vec![0u32; factorbase.len()];

    for (exponent, &p) in factors.iter_mut().zip(factorbase) {
        // Nothing is left to divide.
        if n == 1 {
            break;
        }
        while n.is_divisible_u(p) {
            *exponent += 1;
            n.div_exact_u_mut(p);
        }
    }

    // The number does not factor completely over the factor base.
    (n == 1).then_some(factors)
}

/// Index Calculus algorithm for computing the discrete logarithm of `a` in base `b` modulo `n`.
///
/// The group order must be given and prime. It is not suitable for small orders
/// and the algorithm might fail to find a solution in such situations.
///
/// This algorithm is particularly efficient for large prime orders when
/// exp(2*sqrt(log(n)*log(log(n)))) < sqrt(order).
///
/// # Examples
///
/// ```
/// use discrete_logarithm::discrete_log_index_calculus;
/// use rug::Integer;
///
/// let n = Integer::from(24570203447_u64);
/// let a = Integer::from(23859756228_u64);
/// let b = Integer::from(2);
/// let order = Integer::from(12285101723_u64);
///
/// let x = discrete_log_index_calculus(&n, &a, &b, Some(&order)).unwrap();
/// assert_eq!(x, Integer::from(4519867240_u64));
/// ```
///
/// If the order of the group is known, it must be passed as `order`.
pub fn discrete_log_index_calculus(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
) -> Result<Integer, Error> {
    solve(n, a, b, order, &mut RandState::new())
}

/// Index Calculus algorithm for computing the discrete logarithm of `a` in base `b` modulo `n`.
///
/// Same as [`discrete_log_index_calculus`], with the random generator seeded with `seed`: the
/// same seed always looks for the relations in the same order, and a retry with another seed
/// looks for others.
pub fn discrete_log_index_calculus_with_seed(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
    seed: u64,
) -> Result<Integer, Error> {
    let mut rand_state = RandState::new();
    rand_state.seed(&Integer::from(seed));
    solve(n, a, b, order, &mut rand_state)
}

fn solve(
    n: &Integer,
    a: &Integer,
    b: &Integer,
    order: Option<&Integer>,
    rand_state: &mut RandState<'_>,
) -> Result<Integer, Error> {
    let a = a.clone().modulo(n);
    let b = b.clone().modulo(n);

    let order = match order {
        Some(order) => order.clone(),
        None => return Err(Error::LogDoesNotExist),
    };

    // Compute the bound B for the factorbase using the heuristic from the sympy implementation
    // B = exp(0.5 * sqrt(log(n) * log(log(n))) * (1 + 1/log(log(n))))
    let n_f64 = n.to_f64();
    let log_n = n_f64.ln();
    let log_log_n = log_n.ln();
    let b_bound = (0.5 * (log_n * log_log_n).sqrt() * (1.0 + 1.0 / log_log_n)).exp();
    // Such a factorbase cannot be sieved, let alone held in memory.
    if b_bound > u32::MAX as f64 {
        return Err(Error::LogDoesNotExist);
    }
    let b_bound = b_bound as u32;

    // Compute the factorbase - all primes up to B (exclusive, matching sympy's primerange(B))
    let factorbase: Vec<u32> = Primes::all()
        .take_while(|&p| p < b_bound as usize)
        .map(|p| p as u32)
        .collect();
    let lf = factorbase.len();

    if lf == 0 {
        return Err(Error::LogDoesNotExist);
    }

    // Maximum number of tries to find a relation
    let max_tries = 5 * u128::from(b_bound) * u128::from(b_bound);

    // First, find a relation for a
    let mut relationa: Option<(Vec<Integer>, Integer)> = None;
    let mut abx = a.clone();

    for x in 0..order.to_u64().unwrap_or(u64::MAX) {
        if abx == 1 {
            return Ok((order.clone() - x) % &order);
        }

        if let Some(factors) = is_smooth(abx.clone(), &factorbase) {
            // Convert to Integer and compute modulo order
            let factors_int: Vec<Integer> =
                factors.iter().map(|&f| Integer::from(f) % &order).collect();
            relationa = Some((factors_int, Integer::from(x)));
            break;
        }

        abx = abx * &b % n;
    }

    let (mut relationa, relationa_x) = match relationa {
        Some(r) => r,
        None => return Err(Error::LogDoesNotExist),
    };
    relationa.push(relationa_x);

    // Now find relations for the factorbase elements
    let mut relations: Vec<Option<Vec<Integer>>> = vec![None; lf];
    let mut k = 1; // number of relations found
    let mut kk = 0; // number of consecutive failures

    let order_minus_1: Integer = order.clone() - 1u32;

    while k < 3 * lf && kk < max_tries {
        // Generate random exponent x in [1, order-1]
        let x = order_minus_1.clone().random_below(rand_state) + 1u32;

        // Compute b^x mod n
        let bx = b.clone().pow_mod(&x, n).unwrap();

        // Try to factor it over the factorbase
        let relation = match is_smooth(bx, &factorbase) {
            Some(factors) => {
                let mut rel: Vec<Integer> =
                    factors.iter().map(|&f| Integer::from(f) % &order).collect();
                rel.push(x);
                rel
            }
            None => {
                kk += 1;
                continue;
            }
        };

        k += 1;
        kk = 0;

        // Gaussian elimination step
        let mut relation = relation;
        let mut index = lf; // index of first nonzero entry

        for i in 0..lf {
            let ri = relation[i].clone() % &order;

            if ri > 0 && relations[i].is_some() {
                // Make this entry zero using existing relation
                let existing = relations[i].as_ref().unwrap();
                for j in 0..=lf {
                    let diff = relation[j].clone() - &ri * &existing[j];
                    relation[j] = (diff % &order + &order) % &order;
                }
            } else {
                relation[i] = ri.clone();
            }

            if relation[i] > 0 && index == lf {
                index = i;
            }
        }

        if index == lf || relations[index].is_some() {
            // No new information
            continue;
        }

        // Normalize the relation
        let rinv = relation[index].clone().invert(&order).unwrap();
        for item in relation.iter_mut().skip(index) {
            *item = (rinv.clone() * &*item) % &order;
        }

        relations[index] = Some(relation.clone());

        // Reduce relationa with the new relation
        for i in 0..lf {
            if relationa[i] > 0 && relations[i].is_some() {
                let rbi = relationa[i].clone();
                let existing = relations[i].as_ref().unwrap();
                for j in 0..=lf {
                    let diff = relationa[j].clone() - &rbi * &existing[j];
                    relationa[j] = (diff % &order + &order) % &order;
                }
            }
            if relationa[i] > 0 {
                break; // We have a nonzero entry, don't need to continue reducing
            }
        }

        // Check if all unknowns are eliminated
        let all_zero = (0..lf).all(|i| relationa[i] == 0);
        if all_zero {
            let x = (order.clone() - &relationa[lf]) % &order;

            // Verify the result
            if b.clone().pow_mod(&x, n).unwrap() == a {
                return Ok(x);
            }

            return Err(Error::LogDoesNotExist);
        }
    }

    Err(Error::LogDoesNotExist)
}

#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use rug::ops::Pow;

    use super::*;

    fn int(s: &str) -> Integer {
        Integer::from_str(s).unwrap()
    }

    #[test]
    fn index_calculus() {
        // Test case from sympy documentation
        assert_eq!(
            discrete_log_index_calculus(
                &int("24570203447"),
                &int("23859756228"),
                &2.into(),
                Some(&int("12285101723"))
            )
            .unwrap(),
            int("4519867240")
        );
    }

    #[test]
    fn index_calculus_small() {
        // Small test cases
        assert_eq!(
            discrete_log_index_calculus(
                &587.into(),
                &(Integer::from(2).pow(9)),
                &2.into(),
                Some(&293.into())
            )
            .unwrap(),
            9
        );
    }

    #[test]
    fn sympy_cases() {
        for (n, a, b, order, x) in [
            ("983", "948", "2", "491", "183"),
            ("633383", "21794", "2", "316691", "68048"),
            ("941762639", "68822582", "2", "470881319", "338029275"),
            (
                "999231337607",
                "888188918786",
                "2",
                "499615668803",
                "142811376514",
            ),
            // A base much larger than the primes of the factorbase.
            (
                "47747730623",
                "19410045286",
                "43425105668",
                "645239603",
                "590504662",
            ),
        ] {
            assert_eq!(
                discrete_log_index_calculus(&int(n), &int(a), &int(b), Some(&int(order))).unwrap(),
                int(x),
                "n = {n}"
            );
        }
    }

    #[test]
    fn seeded() {
        // The same seed always looks for the same relations, and the result does not depend on it.
        for seed in [0, 1, 42] {
            assert_eq!(
                discrete_log_index_calculus_with_seed(
                    &int("941762639"),
                    &int("68822582"),
                    &2.into(),
                    Some(&int("470881319")),
                    seed
                )
                .unwrap(),
                int("338029275")
            );
        }
    }

    #[test]
    fn no_order() {
        // The order is required.
        assert_eq!(
            discrete_log_index_calculus(&int("24570203447"), &int("23859756228"), &2.into(), None),
            Err(Error::LogDoesNotExist)
        );
    }

    #[test]
    fn smoothness() {
        let factorbase = [2, 3, 5, 7];
        assert_eq!(
            is_smooth(Integer::from(2 * 2 * 3 * 7 * 7), &factorbase),
            Some(vec![2, 1, 0, 2])
        );
        assert_eq!(is_smooth(Integer::from(1), &factorbase), Some(vec![0; 4]));
        assert_eq!(is_smooth(Integer::from(11), &factorbase), None);
        assert_eq!(is_smooth(Integer::from(2 * 11), &factorbase), None);
    }
}
