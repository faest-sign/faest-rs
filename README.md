# The FAEST Signature Scheme (Rust, Crypto 2023 Version)

This repository contains the intial Rust implementation of FAEST developed for the benchmarks in our paper:

"Publicly Verifiable Zero-Knowledge and Post-Quantum Signatures From VOLE-in-the-Head". By *Carsten Baum, Lennart Braun, Cyprien Delpech de Saint Guilhem, Michael Klooß, Emmanuela Orsini, Lawrence Roy, and Peter Scholl*. [Crypto 2023](https://crypto.iacr.org/2023/). [Full version on ePrint](https://eprint.iacr.org/2023/996).

We have tested the code with Rust v1.67.0. It requires an x86 processor with the AESNI and AVX2 instruction set extensions. To compile and run the benchmarks, execute `cargo bench faest`.

Note that this implements an early version of FAEST. For a full specification and implementations of our improved scheme, see https://faest.info, however, it is no longer compatible with the implementation in this repository.
