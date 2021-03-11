
use std::env;
use std::io::{self,Write,Result};
use std::path::{Path,PathBuf};
use std::fmt;

use fs_err::{self as fs, OpenOptions};
use rand::{self, distributions::Uniform, prelude::Distribution};

include!("src/f2e16.rs");


/// Write Rust `const` declaration
pub fn write_const<W,T>(mut w: W, name: &str, value: &T, type_name: &str) -> Result<()>
where
    W: std::io::Write,
    T: std::fmt::Debug
{
    write!(w, "pub(crate) static {}: {} = {:#?};\n\n", name, type_name, value)
}


/// Compute tables determined solely by the field, which never depend
/// upon the FFT domain or erasure coding paramaters.
///
/// We compute `LOG_TABLE` and `EXP_TABLE` here of course.  We compute
/// the Walsh transform table `LOG_WALSH` here too because we never figured
/// out how to shrink `LOG_WALSH` below the size of the full field (TODO).
/// We thus assume it depends only upon the field for now.
#[allow(dead_code)]
fn write_field_tables<W: std::io::Write>(mut w: W) -> std::io::Result<()> {
    let mut log_table: [GFSymbol; FIELD_SIZE] = [0_u16; FIELD_SIZE];
    let mut exp_table: [GFSymbol; FIELD_SIZE] = [0_u16; FIELD_SIZE];

	let mas: Elt = (1 << FIELD_BITS - 1) - 1;
	let mut state: usize = 1;
	for i in 0_usize..(ONEMASK as usize) {
		exp_table[state] = i as Elt;
		if (state >> FIELD_BITS - 1) != 0 {
			state &= mas as usize;
			state = state << 1_usize ^ GENERATOR as usize;
		} else {
			state <<= 1;
		}
	}
	exp_table[0] = ONEMASK;

	log_table[0] = 0;
	for i in 0..FIELD_BITS {
		for j in 0..(1 << i) {
			log_table[j + (1 << i)] = log_table[j] ^ BASE[i];
		}
	}
	for i in 0..FIELD_SIZE {
		log_table[i] = exp_table[log_table[i] as usize];
	}

	for i in 0..FIELD_SIZE {
		exp_table[log_table[i] as usize] = i as GFSymbol;
	}
	exp_table[ONEMASK as usize] = exp_table[0];

    write_const(&mut w,"LOG_TABLE",&log_table, "[u16; FIELD_SIZE]") ?;
    write_const(&mut w,"EXP_TABLE",&exp_table, "[u16; FIELD_SIZE]") ?;

	// mem_cpy(&mut log_walsh[..], &log_table[..]);
    let log_walsh = log_table.clone();
    let mut log_walsh = unsafe { core::mem::transmute::<_,[Multiplier; FIELD_SIZE]>(log_walsh) };
	log_walsh[0] = Multiplier(0);
	walsh(&mut log_walsh[..], FIELD_SIZE);

    write_const(w, "LOG_WALSH", &log_walsh, "[Multiplier; FIELD_SIZE]")?;
	Ok(())
}



/// Create tables file
///
/// We'll eventually need a seperate tables.rs build target because cargo
/// dislikes build artifacts appearing outside env!("OUT_DIR") and we
/// require tables to build other tables.
/// ref.  https://doc.rust-lang.org/cargo/reference/build-scripts.html#outputs-of-the-build-script
pub fn gen_field_tables() -> std::io::Result<()> {
	// to avoid a circular loop, we need to import a dummy
	// table, such that we do not depend on the thing we are
	// about to spawn
	println!("cargo:rustc-cfg=table_bootstrap_complete");


	let out = env::var("OUT_DIR").expect("OUT_DIR is set by cargo after process launch. qed");

	let path = PathBuf::from(out).join(format!("table_{}.rs", "f2e16"));
	let mut f = OpenOptions::new().create(true).truncate(true).write(true).open(path)?;
	write_field_tables(f)?;

	Ok(())
}

fn gen_10mb_rand_data() -> Result<()> {
	let mut rng = rand::thread_rng();
	let dice = Uniform::<u8>::new_inclusive(0, 255);
	let data = dice.sample_iter(&mut rng).take(10_000_000).collect::<Vec<_>>();

	let out = env::var("OUT_DIR").expect("OUT_DIR is set by cargo after process launch. qed");
	let dest = PathBuf::from(out).join("rand_data.bin");

	let mut f = OpenOptions::new().truncate(true).write(true).create(true).open(&dest)?;

	f.write_all(&data)?;

	f.flush()?;

	Ok(())
}

#[cfg(feature = "cmp-with-cxx")]
fn gen_ffi_novel_poly_basis_lib() {
	cc::Build::new().file("cxx/RSErasureCode.c").file("cxx/sha-256.c").include("cxx").compile("novelpolycxxffi");
}

#[cfg(feature = "cmp-with-cxx")]
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

fn main() -> Result<()> {
    gen_field_tables()?;

	#[cfg(feature = "cmp-with-cxx")]
	{
		gen_ffi_novel_poly_basis_lib();
		gen_ffi_novel_poly_basis_bindgen();
	}
	gen_10mb_rand_data()
}
