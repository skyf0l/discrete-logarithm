use criterion::{criterion_group, criterion_main, Criterion};
use discrete_logarithm::*;
use rug::{ops::Pow, Integer};

fn small_discrete_log(c: &mut Criterion) {
    let mut group = c.benchmark_group("discrete_log_small");

    group.bench_function("trial_mul", |b| {
        let n = Integer::from(587);
        let a = Integer::from(2).pow(9);
        let b_val = Integer::from(2);
        b.iter(|| discrete_log_trial_mul(&n, &a, &b_val, None).unwrap());
    });

    group.bench_function("shanks_steps", |b| {
        let n = Integer::from(442879);
        let a = Integer::from(7).pow(2);
        let b_val = Integer::from(7);
        b.iter(|| discrete_log_shanks_steps(&n, &a, &b_val, None).unwrap());
    });

    group.finish();
}

fn medium_discrete_log(c: &mut Criterion) {
    let mut group = c.benchmark_group("discrete_log_medium");
    group.sample_size(10);

    group.bench_function("pollard_rho", |b| {
        let n = Integer::from(24567899_u64);
        let a = Integer::from(3).pow(333);
        let b_val = Integer::from(3);
        b.iter(|| discrete_log_pollard_rho(&n, &a, &b_val, None).unwrap());
    });

    group.bench_function("pohlig_hellman", |b| {
        let n = Integer::from(98376431_u64);
        let a = Integer::from(11).pow(9);
        let b_val = Integer::from(11);
        b.iter(|| discrete_log_pohlig_hellman(&n, &a, &b_val, None).unwrap());
    });

    group.finish();
}

fn auto_discrete_log(c: &mut Criterion) {
    let mut group = c.benchmark_group("discrete_log_auto");
    group.sample_size(10);

    group.bench_function("auto_small", |b| {
        let n = Integer::from(587);
        let a = Integer::from(2).pow(9);
        let b_val = Integer::from(2);
        b.iter(|| discrete_log(&n, &a, &b_val).unwrap());
    });

    group.bench_function("auto_medium", |b| {
        let n = Integer::from(5779);
        let a = Integer::from(3528);
        let b_val = Integer::from(6215);
        b.iter(|| discrete_log(&n, &a, &b_val).unwrap());
    });

    group.finish();
}

criterion_group!(benches, small_discrete_log, medium_discrete_log, auto_discrete_log);
criterion_main!(benches);
