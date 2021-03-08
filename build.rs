
extern crate derive_more;

use std::env;
use std::io::{Write,Result};
use std::path::{Path,PathBuf};

use fs_err::OpenOptions;
use rand::{self, distributions::Uniform, prelude::Distribution};


/// Create tables file
///
/// We'll eventually need a seperate tables.rs build target because cargo
/// dislikes build artifacts appearing outside env!("OUT_DIR") and we
/// require tables to build other tables.
/// ref.  https://doc.rust-lang.org/cargo/reference/build-scripts.html#outputs-of-the-build-script
fn create_tables(field: impl AsRef<Path>) -> Result<std::fs::File> {
	let mut name = PathBuf::from("src").join(field);
    name.set_extension("tables.rs");
    std::fs::File::create(name)
}

/// Write Rust `const` declaration
pub fn write_const<W,T>(mut w: W, name: &str, value: &T, ty: Option<&str>) -> Result<()> 
where
    W: std::io::Write,
    T: std::fmt::Debug
{
    let type_name = ty.unwrap_or(std::any::type_name::<T>());
    write!(w,"pub(crate) static {}: {} = {:#?};\n\n", name, type_name, value)
}

mod f2e8 { 
}

mod f2e16 {
    use super::write_const;

    include!("src/f2e16.rs");

    pub fn gen_tables() -> std::io::Result<()> {
        let f = super::create_tables("f2e16") ?;
        write_tables(f)
    }    
}



fn gen_10mb_rand_data() -> Result<()> {
	let mut rng = rand::thread_rng();
	let dice = Uniform::<u8>::new_inclusive(0, 255);
	let data = dice.sample_iter(&mut rng).take(10_000_000).collect::<Vec<_>>();

	let dest = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR is set by cargo after process launch. qed"))
		.join("rand_data.bin");

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
    // f2e8::gen_tables() ?;
    f2e16::gen_tables() ?;
    println!("cargo:rustc-cfg=table_build_done");

	#[cfg(feature = "cmp-with-cxx")]
	{
		gen_ffi_novel_poly_basis_lib();
		gen_ffi_novel_poly_basis_bindgen();
	}
	gen_10mb_rand_data()
}
