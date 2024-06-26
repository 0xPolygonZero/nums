# nums

This library contains some number theoretic functions, such as primality testing and factorization, for `BigUint`s.


## Status

Primality tests:
- [x] Trial division
- [x] Miller-Rabin

Factorization:
- [x] Trial division
- [x] Pollard's rho
- [x] Quadratic sieve
    - [ ] SIQS (or other methods to mitigate growth)
    - [ ] Preprocessing to shrink exponent matrix, looking for primes that occur 0, 1 or 2 times
    - [ ] Large prime optimization
    - [ ] Replace Guassian elimination with block-Lanczos or block-Wiedemann nullspace algorithm
- [ ] General number sieve


## License

Licensed under either of

* Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
* MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.


### Licensing

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in the work by you, as defined in the
Apache-2.0 license, shall be dual licensed as above, without any
additional terms or conditions.
