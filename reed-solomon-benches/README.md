# reed-solomon-novelpoly

An implementation of [reed solomon erasure encoding in a novel polynomial base (README)](../reed-solomon-novelpoly/README.md) for speed.


# Validation and Verification
## erasure coding

This repo only exists, to check various erasure coding implementation algorithms in order to determin the most time and space efficient ones.

### test

All benches are also tests with smaller data samples to verify integrity.

```sh
cargo test
```

must always pass.


### criterion

For benchmarking the implementation against itself and the naive implementation,
`cargo criterion` is used.

### bench

```sh
cargo bench
```

will use `valgrind` to run the bench binaries, which will show various metrics, and their changes relative to the previous run.

### flamegraph

```sh
cargo run
```

runs a test case with 10 MB of randomly sampled data which is the recommended way to retrieve a `flamegraph` via `cargo flamegraph` (`cargo install flamegraph` to install).


### fuzzing

Currently `honggfuzz` is used.

To build that a `clang` based toolchain is required.

Install `cargo install honggfuzz` and run with

Run the fuzzer with `cargo hfuzz run fuzzit`.
