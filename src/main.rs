use rs_ec_perf::*;

use rand::distributions::Uniform;
use rand::distributions::Distribution;

fn main() {
	let mut rng = rand::thread_rng();
    let dice = Uniform::<u8>::new_inclusive(0, 255);
	let data = dice.sample_iter(&mut rng)
		.take(10_000_000)
		.collect::<Vec<_>>();
    roundtrip(status_quo_encode, status_quo_reconstruct, &data);
}
