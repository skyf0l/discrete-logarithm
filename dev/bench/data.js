window.BENCHMARK_DATA = {
  "lastUpdate": 1790315007416,
  "repoUrl": "https://github.com/skyf0l/discrete-logarithm",
  "entries": {
    "Instruction counts": [
      {
        "commit": {
          "author": {
            "email": "59019720+skyf0l@users.noreply.github.com",
            "name": "Skyf0l",
            "username": "skyf0l"
          },
          "committer": {
            "email": "59019720+skyf0l@users.noreply.github.com",
            "name": "Skyf0l",
            "username": "skyf0l"
          },
          "distinct": true,
          "id": "5cf3fbd1a89066452f6d520b7f0710b3ca7f4c37",
          "message": "feat(ci): one PR comment with the worst instruction count changes first\n\n  - gungraun-comment.sh: markdown table sorted by change, unchanged benchmarks collapsed\n  - report the results even when the job fails on a regression\n  - build the benchmarks without debug info in CI\n  - bench-lto profile for the wall-clock benchmarks",
          "timestamp": "2026-09-24T19:27:58+04:00",
          "tree_id": "6456314aa0596471e983f55461328ed1bcf22721",
          "url": "https://github.com/skyf0l/discrete-logarithm/commit/5cf3fbd1a89066452f6d520b7f0710b3ca7f4c37"
        },
        "date": 1790264047731,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "e2e::e2e::solve::digits_108",
            "value": 3238003530,
            "unit": "instructions"
          },
          {
            "name": "e2e::e2e::solve::n_2456747",
            "value": 2004053350,
            "unit": "instructions"
          },
          {
            "name": "e2e::e2e::solve::n_32942478",
            "value": 4996647045,
            "unit": "instructions"
          },
          {
            "name": "e2e::e2e::solve::n_5779",
            "value": 2992455152,
            "unit": "instructions"
          },
          {
            "name": "e2e::e2e::solve::n_587",
            "value": 1994911641,
            "unit": "instructions"
          },
          {
            "name": "ops::composite_order::pohlig_hellman::digits_108",
            "value": 1058570300,
            "unit": "instructions"
          },
          {
            "name": "ops::composite_order::pohlig_hellman::n_32942478",
            "value": 997664174,
            "unit": "instructions"
          },
          {
            "name": "ops::factorization::factor::large_prime_cofactor",
            "value": 1021496419,
            "unit": "instructions"
          },
          {
            "name": "ops::factorization::factor::small",
            "value": 997444330,
            "unit": "instructions"
          },
          {
            "name": "ops::factorization::factor::smooth",
            "value": 997495100,
            "unit": "instructions"
          },
          {
            "name": "ops::factorization::index_calculus_smoothness::not_smooth",
            "value": 31270,
            "unit": "instructions"
          },
          {
            "name": "ops::factorization::index_calculus_smoothness::smooth",
            "value": 34739,
            "unit": "instructions"
          },
          {
            "name": "ops::n_order_group::order::composite",
            "value": 3998983425,
            "unit": "instructions"
          },
          {
            "name": "ops::n_order_group::order::prime",
            "value": 2001508344,
            "unit": "instructions"
          },
          {
            "name": "ops::n_order_group::order::prime_power",
            "value": 1994969482,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::index_calculus::n_47747730623",
            "value": 459568627,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::index_calculus::n_633383",
            "value": 4406355,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::index_calculus::n_941762639",
            "value": 118723518,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::index_calculus::n_983",
            "value": 1719807,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::index_calculus::n_999231337607",
            "value": 1705613612,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::prime_order_index_calculus::bits_28",
            "value": 189535747,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::prime_order_index_calculus::bits_34",
            "value": 1398320055,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::prime_order_pollard_rho::bits_28",
            "value": 212784192,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::prime_order_pollard_rho::bits_34",
            "value": 1623899704,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::prime_order_shanks_steps::bits_28",
            "value": 88066730,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::prime_order_shanks_steps::bits_34",
            "value": 493353182,
            "unit": "instructions"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "29139614+renovate[bot]@users.noreply.github.com",
            "name": "renovate[bot]",
            "username": "renovate[bot]"
          },
          "committer": {
            "email": "29139614+renovate[bot]@users.noreply.github.com",
            "name": "renovate[bot]",
            "username": "renovate[bot]"
          },
          "distinct": true,
          "id": "5be443433972e188926ce20290839d4cdcd81538",
          "message": "chore(deps): update rust crate thiserror to v2.0.21",
          "timestamp": "2026-09-25T05:37:13Z",
          "tree_id": "0688e79cce586fc064976a816d7c49693a63b88b",
          "url": "https://github.com/skyf0l/discrete-logarithm/commit/5be443433972e188926ce20290839d4cdcd81538"
        },
        "date": 1790315006915,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "e2e::e2e::solve::digits_108",
            "value": 3238011496,
            "unit": "instructions"
          },
          {
            "name": "e2e::e2e::solve::n_2456747",
            "value": 2004054463,
            "unit": "instructions"
          },
          {
            "name": "e2e::e2e::solve::n_32942478",
            "value": 4996646849,
            "unit": "instructions"
          },
          {
            "name": "e2e::e2e::solve::n_5779",
            "value": 2992455141,
            "unit": "instructions"
          },
          {
            "name": "e2e::e2e::solve::n_587",
            "value": 1994911641,
            "unit": "instructions"
          },
          {
            "name": "ops::composite_order::pohlig_hellman::digits_108",
            "value": 1058522762,
            "unit": "instructions"
          },
          {
            "name": "ops::composite_order::pohlig_hellman::n_32942478",
            "value": 997663158,
            "unit": "instructions"
          },
          {
            "name": "ops::factorization::factor::large_prime_cofactor",
            "value": 1021496420,
            "unit": "instructions"
          },
          {
            "name": "ops::factorization::factor::small",
            "value": 997444331,
            "unit": "instructions"
          },
          {
            "name": "ops::factorization::factor::smooth",
            "value": 997495156,
            "unit": "instructions"
          },
          {
            "name": "ops::factorization::index_calculus_smoothness::not_smooth",
            "value": 31271,
            "unit": "instructions"
          },
          {
            "name": "ops::factorization::index_calculus_smoothness::smooth",
            "value": 34740,
            "unit": "instructions"
          },
          {
            "name": "ops::n_order_group::order::composite",
            "value": 3998983487,
            "unit": "instructions"
          },
          {
            "name": "ops::n_order_group::order::prime",
            "value": 2001508345,
            "unit": "instructions"
          },
          {
            "name": "ops::n_order_group::order::prime_power",
            "value": 1994969483,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::index_calculus::n_47747730623",
            "value": 459568624,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::index_calculus::n_633383",
            "value": 4406352,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::index_calculus::n_941762639",
            "value": 118723515,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::index_calculus::n_983",
            "value": 1719804,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::index_calculus::n_999231337607",
            "value": 1705613609,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::prime_order_index_calculus::bits_28",
            "value": 189535747,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::prime_order_index_calculus::bits_34",
            "value": 1398320055,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::prime_order_pollard_rho::bits_28",
            "value": 212784192,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::prime_order_pollard_rho::bits_34",
            "value": 1623899704,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::prime_order_shanks_steps::bits_28",
            "value": 88069244,
            "unit": "instructions"
          },
          {
            "name": "ops::prime_order::prime_order_shanks_steps::bits_34",
            "value": 493298409,
            "unit": "instructions"
          }
        ]
      }
    ]
  }
}