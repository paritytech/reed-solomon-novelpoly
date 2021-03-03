use fs_err::OpenOptions;
use rand::{self, distributions::Uniform, prelude::Distribution};
use std::env;
use std::io::Write;
use std::path::PathBuf;

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

fn gen_ffi_novel_poly_basis_lib() {
	cc::Build::new().file("cxx/RSErasureCode.c").file("cxx/sha-256.c").include("cxx").compile("novelpolycxxffi");
}

fn gen_ffi_novel_poly_basis_bindgen() {
	println!("cargo:rustc-link-lib=novelpolycxxffi");

	println!("cargo:rerun-if-changed=wrapper.h");

	let bindings = bindgen::Builder::default()
		.header("wrapper.h")
		.parse_callbacks(Box::new(bindgen::CargoCallbacks))
		.generate()
		.expect("Unable to generate bindings");

	// Write the bindings to the $OUT_DIR/bindings.rs file.
	let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
	bindings.write_to_file(out_path.join("bindings.rs")).expect("Couldn't write bindings!");
}

fn main() -> Result<(), std::io::Error> {
	gen_ffi_novel_poly_basis_lib();
	gen_ffi_novel_poly_basis_bindgen();
	gen_10mb_rand_data()
}
