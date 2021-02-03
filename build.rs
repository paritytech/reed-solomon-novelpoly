use fs_err::OpenOptions;
use rand::{self, distributions::Uniform, prelude::Distribution};
use std::env;
use std::io::Write;

fn gen_10mb_rand_data() -> Result<(), std::io::Error> {
	let mut rng = rand::thread_rng();
	let dice = Uniform::<u8>::new_inclusive(0, 255);
	let data = dice.sample_iter(&mut rng).take(10_000_000).collect::<Vec<_>>();

	let dest =
		std::path::PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR is set by cargo after process launch. qed"))
			.join("rand_data.bin");

	let mut f = OpenOptions::new().truncate(true).write(true).create(true).open(&dest)?;

	f.write_all(&data)?;

	f.flush()?;
	Ok(())
}

fn main() -> Result<(), std::io::Error> {
	gen_10mb_rand_data()
}
