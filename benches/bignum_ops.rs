use criterion::{criterion_group, criterion_main, Criterion};

// Test values for benchmarking - large primes
const P: &str = "13407807929942597099574024998205846127479365820592393377723561443721764030073546976801874298166903427690031858186486050853753882811946569946433649006084171";
const G: &str = "11717829880366207009516117596335367088558084999998952205599979459063929499736583746670572176471460312928594829675428279466566527115212748467589894601965568";
const H: &str = "3239475104050450443565264378728065788649097520952449527834792452971981976143292558073856937958553180532878928001494706097394108577585732452307673444020333";

fn mulmod(c: &mut Criterion) {
    let mut group = c.benchmark_group("mulmod");
    
    // rug (default)
    group.bench_function("rug", |b| {
        use rug::Integer;
        let p = Integer::from_str_radix(P, 10).unwrap();
        let g = Integer::from_str_radix(G, 10).unwrap();
        let h = Integer::from_str_radix(H, 10).unwrap();
        b.iter(|| Integer::from(&g * &h) % &p);
    });

    #[cfg(feature = "bench-num-bigint")]
    group.bench_function("num-bigint", |b| {
        use num_bigint::BigInt;
        use num_traits::Num;
        let p = BigInt::from_str_radix(P, 10).unwrap();
        let g = BigInt::from_str_radix(G, 10).unwrap();
        let h = BigInt::from_str_radix(H, 10).unwrap();
        b.iter(|| (&g * &h) % &p);
    });

    #[cfg(feature = "bench-ramp")]
    group.bench_function("ramp", |b| {
        use ramp::Int;
        let p = Int::from_str_radix(P, 10).unwrap();
        let g = Int::from_str_radix(G, 10).unwrap();
        let h = Int::from_str_radix(H, 10).unwrap();
        b.iter(|| (&g * &h) % &p);
    });

    #[cfg(feature = "bench-ibig")]
    group.bench_function("ibig", |b| {
        use ibig::IBig;
        let p = IBig::from_str_radix(P, 10).unwrap();
        let g = IBig::from_str_radix(G, 10).unwrap();
        let h = IBig::from_str_radix(H, 10).unwrap();
        b.iter(|| (&g * &h) % &p);
    });

    #[cfg(feature = "bench-rust-gmp")]
    group.bench_function("rust-gmp", |b| {
        use gmp::mpz::Mpz;
        let p = Mpz::from_str_radix(P, 10).unwrap();
        let g = Mpz::from_str_radix(G, 10).unwrap();
        let h = Mpz::from_str_radix(H, 10).unwrap();
        b.iter(|| (&g * &h) % &p);
    });

    group.finish();
}

fn powmod(c: &mut Criterion) {
    let mut group = c.benchmark_group("powmod");
    
    // rug (default)
    group.bench_function("rug", |b| {
        use rug::Integer;
        let p = Integer::from_str_radix(P, 10).unwrap();
        let g = Integer::from_str_radix(G, 10).unwrap();
        let exp = Integer::from(65537);
        b.iter(|| g.clone().pow_mod(&exp, &p).unwrap());
    });

    #[cfg(feature = "bench-num-bigint")]
    group.bench_function("num-bigint", |b| {
        use num_bigint::BigInt;
        use num_traits::Num;
        let p = BigInt::from_str_radix(P, 10).unwrap();
        let g = BigInt::from_str_radix(G, 10).unwrap();
        let exp = BigInt::from(65537);
        b.iter(|| g.modpow(&exp, &p));
    });

    #[cfg(feature = "bench-ramp")]
    group.bench_function("ramp", |b| {
        use ramp::Int;
        let p = Int::from_str_radix(P, 10).unwrap();
        let g = Int::from_str_radix(G, 10).unwrap();
        let exp = Int::from(65537);
        b.iter(|| g.pow_mod(&exp, &p));
    });

    // ibig doesn't have a straightforward powmod API, skipping for now
    #[cfg(feature = "bench-ibig")]
    group.bench_function("ibig", |b| {
        use ibig::UBig;
        // Use a simple repeated squaring for power mod
        let p = UBig::from_str_radix(P, 10).unwrap();
        let g = UBig::from_str_radix(G, 10).unwrap();
        let exp = 65537u32;
        b.iter(|| {
            let mut result = UBig::from(1u32);
            let mut base = g.clone();
            let mut e = exp;
            while e > 0 {
                if e % 2 == 1 {
                    result = (&result * &base) % &p;
                }
                base = (&base * &base) % &p;
                e /= 2;
            }
            result
        });
    });

    #[cfg(feature = "bench-rust-gmp")]
    group.bench_function("rust-gmp", |b| {
        use gmp::mpz::Mpz;
        let p = Mpz::from_str_radix(P, 10).unwrap();
        let g = Mpz::from_str_radix(G, 10).unwrap();
        let exp = Mpz::from(65537);
        b.iter(|| g.powm(&exp, &p));
    });

    group.finish();
}

criterion_group!(benches, mulmod, powmod);
criterion_main!(benches);
