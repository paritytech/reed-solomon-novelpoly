// Encoding/erasure decoding for Reed-Solomon codes over binary extension fields
//
// Derived impl of `RSAErasureCode.c`.
//
// Lin, Han and Chung, "Novel Polynomial Basis and Its Application to Reed-Solomon Erasure Codes," FOCS14.
// (http://arxiv.org/abs/1404.3458)

#![allow(dead_code)]

use super::*;

use std::{
	io::{BufRead, Read},
	slice::from_raw_parts,
};

pub type GFSymbol = u16;

pub const FIELD_BITS: usize = 16;

pub const GENERATOR: GFSymbol = 0x2D; //x^16 + x^5 + x^3 + x^2 + 1

// Cantor basis
pub const BASE: [GFSymbol; FIELD_BITS] =
	[1_u16, 44234, 15374, 5694, 50562, 60718, 37196, 16402, 27800, 4312, 27250, 47360, 64952, 64308, 65336, 39198];

pub const FIELD_SIZE: usize = 1_usize << FIELD_BITS;

pub const MODULO: GFSymbol = (FIELD_SIZE - 1) as GFSymbol;

static mut LOG_TABLE: [GFSymbol; FIELD_SIZE] = [0_u16; FIELD_SIZE];
static mut EXP_TABLE: [GFSymbol; FIELD_SIZE] = [0_u16; FIELD_SIZE];

//-----Used in decoding procedure-------
//twisted factors used in FFT
static mut SKEW_FACTOR: [GFSymbol; MODULO as usize] = [0_u16; MODULO as usize];

//factors used in formal derivative
static mut B: [GFSymbol; FIELD_SIZE >> 1] = [0_u16; FIELD_SIZE >> 1];

//factors used in the evaluation of the error locator polynomial
static mut LOG_WALSH: [GFSymbol; FIELD_SIZE] = [0_u16; FIELD_SIZE];

//return a*EXP_TABLE[b] over GF(2^r)
pub fn mul_table(a: GFSymbol, b: GFSymbol) -> GFSymbol {
	if a != 0_u16 {
		unsafe {
			let offset = (LOG_TABLE[a as usize] as u32 + b as u32 & MODULO as u32)
				+ (LOG_TABLE[a as usize] as u32 + b as u32 >> FIELD_BITS);
			EXP_TABLE[offset as usize]
		}
	} else {
		0_u16
	}
}

pub const fn log2(mut x: usize) -> usize {
	let mut o: usize = 0;
	while x > 1 {
		x >>= 1;
		o += 1;
	}
	o
}

pub const fn is_power_of_2(x: usize) -> bool {
	x > 0_usize && x & (x - 1) == 0
}

//fast Walsh–Hadamard transform over modulo mod
pub fn walsh(data: &mut [GFSymbol], size: usize) {
	let mut depart_no = 1_usize;
	while depart_no < size {
		let mut j = 0;
		let depart_no_next = depart_no << 1;
		while j < size {
			for i in j..(depart_no + j) {
				let tmp2: u32 = data[i] as u32 + MODULO as u32 - data[i + depart_no] as u32;
				data[i] = ((data[i] as u32 + data[i + depart_no] as u32 & MODULO as u32)
					+ (data[i] as u32 + data[i + depart_no] as u32 >> FIELD_BITS)) as GFSymbol;
				data[i + depart_no] = ((tmp2 & MODULO as u32) + (tmp2 >> FIELD_BITS)) as GFSymbol;
			}
			j += depart_no_next;
		}
		depart_no = depart_no_next;
	}
}

//formal derivative of polynomial in the new basis
pub fn formal_derivative(cos: &mut [GFSymbol], size: usize) {
	for i in 1..size {
		let length = ((i ^ i - 1) + 1) >> 1;
		for j in (i - length)..i {
			cos[j] ^= cos.get(j + length).copied().unwrap_or_default();
		}
	}
	let mut i = size;
	while i < FIELD_SIZE && i < cos.len() {
		for j in 0..size {
			cos[j] ^= cos.get(j + i).copied().unwrap_or_default();
		}
		i <<= 1;
	}
}

// We want the low rate scheme given in
// https://www.citi.sinica.edu.tw/papers/whc/5524-F.pdf
// and https://github.com/catid/leopard/blob/master/docs/LowRateDecoder.pdf
// but this code resembles https://github.com/catid/leopard which
// implements the high rate decoder in
// https://github.com/catid/leopard/blob/master/docs/HighRateDecoder.pdf
// We're hunting for the differences and trying to undersrtand the algorithm.

//IFFT in the proposed basis
pub fn inverse_fft_in_novel_poly_basis(data: &mut [GFSymbol], size: usize, index: usize) {
	// All line references to Algorithm 2 page 6288 of
	// https://www.citi.sinica.edu.tw/papers/whc/5524-F.pdf

	// Depth of the recursion on line 7 and 8 is given by depart_no
	// aka 1 << ((k of Algorithm 2) - (i of Algorithm 2)) where
	// k of Algorithm 1 is read as FIELD_BITS here.
	// Recusion base layer implicitly imports d_r aka ala line 1.
	// After this, we start at depth (i of Algorithm 2) = (k of Algorithm 2) - 1
	// and progress through FIELD_BITS-1 steps, obtaining \Psi_\beta(0,0).
	let mut depart_no = 1_usize;
	while depart_no < size {
		// Agrees with for loop (j of Algorithm 2) in (0..2^{k-i-1}) from line 3,
		// except we've j in (depart_no..size).step_by(2*depart_no), meaning
		// the doubled step compensated for the halve size exponent, and
		// somehow this j captures the subscript on \omega_{j 2^{i+1}}.	 (TODO)
		let mut j = depart_no;
		while j < size {
			// At this point loops over i in (j - depart_no)..j give a bredth
			// first loop across the recursion branches from lines 7 and 8,
			// so the i loop corresponds to r in Algorithm 2.  In fact,
			// data[i] and data[i + depart_no] together cover everything,
			// thanks to the outer j loop.

			// Loop on line 3, so i corresponds to j in Algorithm 2
			for i in (j - depart_no)..j {
				// Line 4, justified by (34) page 6288, but
				// adding depart_no acts like the r+2^i superscript.
				data[i + depart_no] ^= data[i];
			}

			// Algorithm 2 indexs the skew factor in line 5 page 6288
			// by i and \omega_{j 2^{i+1}}, but not by r explicitly.
			// We further explore this confusion below. (TODO)
			let skew = unsafe { SKEW_FACTOR[j + index - 1] };
			// It's reasonale to skip the loop if skew is zero, but doing so with
			// all bits set requires justification.	 (TODO)
			if skew != MODULO {
				// Again loop on line 3, except skew should depend upon i aka j in Algorithm 2 (TODO)
				for i in (j - depart_no)..j {
					// Line 5, justified by (35) page 6288, but
					// adding depart_no acts like the r+2^i superscript.
					data[i] ^= mul_table(data[i + depart_no], skew);
				}
			}

			// Increment by double depart_no in agreement with
			// our updating 2*depart_no elements at this depth.
			j += depart_no << 1;
		}
		depart_no <<= 1;
	}
}

//FFT in the proposed basis
pub fn fft_in_novel_poly_basis(data: &mut [GFSymbol], size: usize, index: usize) {
	// All line references to Algorithm 1 page 6287 of
	// https://www.citi.sinica.edu.tw/papers/whc/5524-F.pdf

	// Depth of the recursion on line 3 and 4 is given by depart_no
	// aka 1 << ((k of Algorithm 1) - (i of Algorithm 1)) where
	// k of Algorithm 1 is read as FIELD_BITS here.
	// Recusion base layer implicitly imports d_r aka ala line 1.
	// After this, we start at depth (i of Algorithm 1) = (k of Algorithm 1) - 1
	// and progress through FIELD_BITS-1 steps, obtaining \Psi_\beta(0,0).
	let mut depart_no = size >> 1_usize;
	while depart_no > 0 {
		// Agrees with for loop (j of Algorithm 1) in (0..2^{k-i-1}) from line 5,
		// except we've j in (depart_no..size).step_by(2*depart_no), meaning
		// the doubled step compensated for the halve size exponent, and
		// somehow this j captures the subscript on \omega_{j 2^{i+1}}.	 (TODO)
		let mut j = depart_no;
		while j < size {
			// At this point loops over i in (j - depart_no)..j give a bredth
			// first loop across the recursion branches from lines 3 and 4,
			// so the i loop corresponds to r in Algorithm 1.  In fact,
			// data[i] and data[i + depart_no] together cover everything,
			// thanks to the outer j loop.

			// Algorithm 1 indexs the skew factor in line 6 aka (28) page 6287
			// by i and \omega_{j 2^{i+1}}, but not by r explicitly.
			// We doubt the lack of explicit dependence upon r justifies
			// extracting the skew factor outside the loop here.
			// As indexing by \omega_{j 2^{i+1}} appears absolute elsewhere,
			// we think r actually appears but the skew factor repeats itself
			// like in (19) in the proof of Lemma 4.  (TODO)
			// We should understand the rest of this basis story, like (8) too.	 (TODO)
			let skew = unsafe { SKEW_FACTOR[j + index - 1] };
			// It's reasonale to skip the loop if skew is zero, but doing so with
			// all bits set requires justification.	 (TODO)
			if skew != MODULO {
				// Loop on line 5, except skew should depend upon i aka j in Algorithm 1 (TODO)
				for i in (j - depart_no)..j {
					// Line 6, explained by (28) page 6287, but
					// adding depart_no acts like the r+2^i superscript.
					data[i] ^= mul_table(data[i + depart_no], skew);
				}
			}

			// Again loop on line 5, so i corresponds to j in Algorithm 1
			for i in (j - depart_no)..j {
				// Line 7, explained by (31) page 6287, but
				// adding depart_no acts like the r+2^i superscript.
				data[i + depart_no] ^= data[i];
			}

			// Increment by double depart_no in agreement with
			// our updating 2*depart_no elements at this depth.
			j += depart_no << 1;
		}
		depart_no >>= 1;
	}
}

//initialize LOG_TABLE[], EXP_TABLE[]
unsafe fn init() {
	let mas: GFSymbol = (1 << FIELD_BITS - 1) - 1;
	let mut state: usize = 1;
	for i in 0_usize..(MODULO as usize) {
		EXP_TABLE[state] = i as GFSymbol;
		if (state >> FIELD_BITS - 1) != 0 {
			state &= mas as usize;
			state = state << 1_usize ^ GENERATOR as usize;
		} else {
			state <<= 1;
		}
	}
	EXP_TABLE[0] = MODULO;

	LOG_TABLE[0] = 0;
	for i in 0..FIELD_BITS {
		for j in 0..(1 << i) {
			LOG_TABLE[j + (1 << i)] = LOG_TABLE[j] ^ BASE[i];
		}
	}
	for i in 0..FIELD_SIZE {
		LOG_TABLE[i] = EXP_TABLE[LOG_TABLE[i] as usize];
	}

	for i in 0..FIELD_SIZE {
		EXP_TABLE[LOG_TABLE[i] as usize] = i as GFSymbol;
	}
	EXP_TABLE[MODULO as usize] = EXP_TABLE[0];
}

//initialize SKEW_FACTOR[], B[], LOG_WALSH[]
unsafe fn init_dec() {
	let mut base: [GFSymbol; FIELD_BITS - 1] = Default::default();

	for i in 1..FIELD_BITS {
		base[i - 1] = 1 << i;
	}

	// We construct SKW_FACTOR to be \bar{s}_j(omega) from page 6285
	// for all omega in the field.
	for m in 0..(FIELD_BITS - 1) {
		let step = 1 << (m + 1);
		SKEW_FACTOR[(1 << m) - 1] = 0;
		for i in m..(FIELD_BITS - 1) {
			let s = 1 << (i + 1);

			let mut j = (1 << m) - 1;
			while j < s {
				// Justified by (5) page 6285, except..
				// we expect SKEW_FACTOR[j ^ field_base[i]] or similar
				SKEW_FACTOR[j + s] = SKEW_FACTOR[j] ^ base[i];
				j += step;
			}
		}

		let idx = mul_table(base[m], LOG_TABLE[(base[m] ^ 1_u16) as usize]);
		base[m] = MODULO - LOG_TABLE[idx as usize];

		for i in (m + 1)..(FIELD_BITS - 1) {
			let b = LOG_TABLE[(base[i] as u16 ^ 1_u16) as usize] as u32 + base[m] as u32;
			let b = b % MODULO as u32;
			base[i] = mul_table(base[i], b as u16);
		}
	}
	for i in 0..(MODULO as usize) {
		SKEW_FACTOR[i] = LOG_TABLE[SKEW_FACTOR[i] as usize];
	}

	base[0] = MODULO - base[0];
	for i in 1..(FIELD_BITS - 1) {
		base[i] = ((MODULO as u32 - base[i] as u32 + base[i - 1] as u32) % MODULO as u32) as GFSymbol;
	}

	B[0] = 0;
	for i in 0..(FIELD_BITS - 1) {
		let depart = 1 << i;
		for j in 0..depart {
			B[j + depart] = ((B[j] as u32 + base[i] as u32) % MODULO as u32) as GFSymbol;
		}
	}

	mem_cpy(&mut LOG_WALSH[..], &LOG_TABLE[..]);
	LOG_WALSH[0] = 0;
	walsh(&mut LOG_WALSH[..], FIELD_SIZE);
}

/// Setup both decoder and encoder.
pub fn setup() {
	use std::sync::Once;

	static SETUP: Once = Once::new();

	SETUP.call_once(|| unsafe {
		init();
		init_dec();
	});
}

// Encoding alg for k/n < 0.5: message is a power of two
pub fn encode_low(data: &[GFSymbol], k: usize, codeword: &mut [GFSymbol], n: usize) {
	assert!(k + k <= n);
	assert_eq!(codeword.len(), n);
	assert_eq!(data.len(), n);

	assert!(is_power_of_2(n));
	assert!(is_power_of_2(k));

	// k | n is guaranteed
	assert_eq!((n / k) * k, n);

	// move the data to the codeword
	mem_cpy(&mut codeword[0..], &data[0..]);

	// split after the first k
	let (codeword_first_k, codeword_skip_first_k) = codeword.split_at_mut(k);

	inverse_fft_in_novel_poly_basis(codeword_first_k, k, 0);

	// the first codeword is now the basis for the remaining transforms
	// denoted `M_topdash`

	for shift in (k..n).into_iter().step_by(k) {
		let codeword_at_shift = &mut codeword_skip_first_k[(shift - k)..shift];
		// copy `M_topdash` to the position we are currently at, the n transform
		mem_cpy(codeword_at_shift, codeword_first_k);
		fft_in_novel_poly_basis(codeword_at_shift, k, shift);
	}

	// restore `M` from the derived ones
	mem_cpy(&mut codeword[0..k], &data[0..k]);
}

fn mem_zero(zerome: &mut [GFSymbol]) {
	for i in 0..zerome.len() {
		zerome[i] = 0_u16;
	}
}

fn mem_cpy(dest: &mut [GFSymbol], src: &[GFSymbol]) {
	let sl = src.len();
	debug_assert_eq!(dest.len(), sl);
	for i in 0..sl {
		dest[i] = src[i];
	}
}

//data: message array. parity: parity array. mem: buffer(size>= n-k)
//Encoding alg for k/n>0.5: parity is a power of two.
pub fn encode_high(data: &[GFSymbol], k: usize, parity: &mut [GFSymbol], mem: &mut [GFSymbol], n: usize) {
	let t: usize = n - k;

	mem_zero(&mut parity[0..t]);

	let mut i = t;
	while i < n {
		mem_cpy(&mut mem[..t], &data[(i - t)..t]);

		inverse_fft_in_novel_poly_basis(mem, t, i);
		for j in 0..t {
			parity[j] ^= mem[j];
		}
		i += t;
	}
	fft_in_novel_poly_basis(parity, t, 0);
}

// Compute the evaluations of the error locator polynomial
// `fn decode_init`
// since this has only to be called once per reconstruction
pub fn eval_error_polynomial(erasure: &[bool], log_walsh2: &mut [GFSymbol], n: usize) {
	let z = std::cmp::min(n, erasure.len());
	for i in 0..z {
		log_walsh2[i] = erasure[i] as GFSymbol;
	}
	for i in z..n {
		log_walsh2[i] = 0 as GFSymbol;
	}
	walsh(log_walsh2, FIELD_SIZE);
	for i in 0..n {
		let tmp = log_walsh2[i] as u32 * unsafe { LOG_WALSH[i] } as u32;
		log_walsh2[i] = (tmp % MODULO as u32) as GFSymbol;
	}
	walsh(log_walsh2, FIELD_SIZE);
	for i in 0..z {
		if erasure[i] {
			log_walsh2[i] = MODULO - log_walsh2[i];
		}
	}
}

/// recover determines how many shards to recover (starting from 0)
// technically we only need to recover
// the first `k` instead of all `n` which
// would include parity chunks.
fn decode_main(codeword: &mut [GFSymbol], recover_up_to: usize, erasure: &[bool], log_walsh2: &[GFSymbol], n: usize) {
	assert_eq!(codeword.len(), n);
	assert!(n >= recover_up_to);
	assert_eq!(erasure.len(), n);

	for i in 0..n {
		codeword[i] = if erasure[i] { 0_u16 } else { mul_table(codeword[i], log_walsh2[i]) };
	}

	inverse_fft_in_novel_poly_basis(codeword, n, 0);

	// formal derivative
	for i in (0..n).into_iter().step_by(2) {
		let b = MODULO - unsafe { B[i >> 1] };
		codeword[i] = mul_table(codeword[i], b);
		codeword[i + 1] = mul_table(codeword[i + 1], b);
	}

	formal_derivative(codeword, n);

	for i in (0..n).into_iter().step_by(2) {
		let b = unsafe { B[i >> 1] };
		codeword[i] = mul_table(codeword[i], b);
		codeword[i + 1] = mul_table(codeword[i + 1], b);
	}

	fft_in_novel_poly_basis(codeword, n, 0);

	for i in 0..recover_up_to {
		codeword[i] = if erasure[i] { mul_table(codeword[i], log_walsh2[i]) } else { 0_u16 };
	}
}
use itertools::Itertools;

pub const fn next_higher_power_of_2(k: usize) -> usize {
	if !is_power_of_2(k) {
		1 << (log2(k) + 1)
	} else {
		k
	}
}

pub const fn next_lower_power_of_2(k: usize) -> usize {
	if !is_power_of_2(k) {
		1 << log2(k)
	} else {
		k
	}
}

#[derive(Debug, Clone, Copy)]
struct ReedSolomon {
	/// total number of message symbols to send
	n: usize,
	/// number of information containing chunks
	k: usize,
}

impl ReedSolomon {
	/// Create a new reed solomon erasure encoding wrapper
	fn new(validator_count: usize) -> ReedSolomon {
		// we need to be able to reconstruct from 1/3 - eps
		let k = validator_count / 3;
		let k = next_lower_power_of_2(k);
		let n = next_higher_power_of_2(validator_count);
		Self { n, k }
	}
}

pub fn encode(bytes: &[u8]) -> Vec<WrappedShard> {
	setup();
	let validator_count = N_VALIDATORS;
	dbg!((bytes.len(), N_VALIDATORS));
	let rs = ReedSolomon::new(validator_count);

	// setup the shards, n is likely _larger_, so use the truely required number of shards

	// shard length in GF(2^16) symbols
	let shard_len = ((bytes.len() + 1) / 2 + rs.k - 1) / rs.k;
	// collect all sub encoding runs

	let k2 = rs.k * 2;
	// prepare one wrapped shard per validator
	let mut shards = vec![
		WrappedShard::new({
			let mut v = Vec::<u8>::with_capacity(shard_len * 2);
			unsafe { v.set_len(shard_len * 2) }
			v
		});
		validator_count
	];

	for (chunk_idx, i) in (0..bytes.len()).into_iter().step_by(k2).enumerate() {
		let end = std::cmp::min(i + k2, bytes.len());
		dbg!((i, end));
		assert_ne!(i, end);
		assert_ne!(i + 1, end);
		let mut encoding_run = encode_sub(&bytes[i..end], rs.n, rs.k);
		for val_idx in 0..validator_count {
			AsMut::<[[u8; 2]]>::as_mut(&mut shards[val_idx])[chunk_idx] = encoding_run[val_idx].to_be_bytes();
		}
	}

	shards
}

/// each shard contains one symbol of one run of erasure coding
pub fn reconstruct(received_shards: Vec<Option<WrappedShard>>) -> Option<Vec<u8>> {
	setup();
	let validator_count = N_VALIDATORS;
	let rs = ReedSolomon::new(validator_count);

	// obtain a sample of a shard length and assume that is the truth
	// XXX make sure all shards have equal length
	let shard_len = received_shards
		.iter()
		.find_map(|x| {
			x.as_ref().map(|x| {
				let x = AsRef::<[[u8; 2]]>::as_ref(x);
				x.len()
			})
		})
		.unwrap();

	let mut acc = Vec::<u8>::with_capacity(shard_len * 2 * rs.k);
	for i in 0..shard_len {
		// take the i-th element of all shards and try to recover
		let mut decoding_run = received_shards
			.iter()
			.map(|x| {
				x.as_ref().map(|x| {
					let z = AsRef::<[[u8; 2]]>::as_ref(&x)[i];
					u16::from_be_bytes(z)
				})
			})
			.chain(
				// reconstruct_sub expects a `rs.n` length symbol slice
				std::iter::repeat(None).take((rs.n).saturating_sub(received_shards.len())),
			)
			.collect::<Vec<Option<GFSymbol>>>();

		assert_eq!(decoding_run.len(), rs.n);

		// reconstruct from one set of symbols which was spread over all erasure chunks
		let piece = reconstruct_sub(&decoding_run[..], rs.n, rs.k).unwrap();
		acc.extend_from_slice(&piece[..]);
	}
	Some(acc)
}

/// Bytes shall only contain payload data
pub fn encode_sub(bytes: &[u8], n: usize, k: usize) -> Vec<GFSymbol> {
	assert!(is_power_of_2(n), "Algorithm only works for 2^i sizes for N");
	assert!(is_power_of_2(k), "Algorithm only works for 2^i sizes for K");
	assert!(bytes.len() <= k << 1);
	assert!(k <= n / 2);

	// must be power of 2
	let dl = bytes.len();
	let l = if is_power_of_2(dl) {
		dl
	} else {
		let l = log2(dl);
		let l = 1 << l;
		let l = if l >= dl { l } else { l << 1 };
		l
	};
	assert!(is_power_of_2(l));
	assert!(l >= dl);

	// pad the incoming bytes with trailing 0s
	// so we get a buffer of size `N` in `GF` symbols
	let zero_bytes_to_add = n * 2 - dl;
	let data: Vec<GFSymbol> = bytes
		.into_iter()
		.copied()
		.chain(std::iter::repeat(0u8).take(zero_bytes_to_add))
		.tuple_windows()
		.step_by(2)
		.map(|(a, b)| u16::from_le_bytes([a, b]))
		.collect::<Vec<GFSymbol>>();

	println!("Data:");
	for data in data.iter() {
		print!("{:04x} ", data);
	}
	println!("");

	// update new data bytes with zero padded bytes
	// `l` is now `GF(2^16)` symbols
	let l = data.len();
	assert_eq!(l, n, "Zer0 padding to full buffer size always works. qed");

	let mut codeword = data.clone();
	assert_eq!(codeword.len(), n);

	encode_low(&data[..], k, &mut codeword[..], n);

	println!("Codeword:");
	for codeword in codeword.iter() {
		print!("{:04x} ", codeword);
	}
	println!("");

	codeword
}

pub fn reconstruct_sub(codewords: &[Option<GFSymbol>], n: usize, k: usize) -> Option<Vec<u8>> {
	assert!(is_power_of_2(n), "Algorithm only works for 2^i sizes for N");
	assert!(is_power_of_2(k), "Algorithm only works for 2^i sizes for K");
	assert_eq!(codewords.len(), n);
	assert!(k <= n / 2);

	// collect all `None` values
	let mut existential_count = 0;
	let erasures = codewords
		.iter()
		.map(|x| x.is_some())
		.inspect(|v| {
			if !*v {
				existential_count += 1;
			}
		})
		.map(|x| !x)
		.collect::<Vec<bool>>();

	assert!(existential_count <= n);

	if existential_count < k {
		println!("Nothing exists");
		return None;
	}

	// the first k suffice for the original k message codewords
	let recover_up_to = n; // k;

	// The recovered _data_ chunks AND parity chunks
	let mut recovered = vec![0 as GFSymbol; recover_up_to];

	// get rid of all `None`s
	let mut codeword = codewords
		.into_iter()
		.enumerate()
		.map(|(idx, sym)| {
			// fill the gaps with `0_u16` codewords
			if let Some(sym) = sym {
				(idx, *sym)
			} else {
				(idx, 0_u16)
			}
		})
		.map(|(idx, codeword)| {
			if idx < recovered.len() {
				recovered[idx] = codeword;
			}
			codeword
		})
		.collect::<Vec<u16>>();

	// filled up the remaining spots with 0s
	assert_eq!(codeword.len(), n);

	let recover_up_to = k; // the first k would suffice for the original k message codewords

	//---------Erasure decoding----------------
	let mut log_walsh2: [GFSymbol; FIELD_SIZE] = [0_u16; FIELD_SIZE];

	// Evaluate error locator polynomial
	eval_error_polynomial(&erasures[..], &mut log_walsh2[..], FIELD_SIZE);

	//---------main processing----------
	decode_main(&mut codeword[..], recover_up_to, &erasures[..], &log_walsh2[..], n);

	println!("Decoded result:");
	for idx in 0..recover_up_to {
		print!("{:04x} ", codeword[idx]);
		if erasures[idx] {
			recovered[idx] = codeword[idx];
		};
	}

	let recovered = unsafe {
		// TODO assure this does not leak memory
		let x = from_raw_parts(recovered.as_ptr() as *const u8, recover_up_to * 2);
		std::mem::forget(recovered);
		x
	};
	let k2 = k * 2;
	Some(recovered[0..k2].to_vec())
}

#[cfg(test)]
mod test {
	use rand::seq::index::IndexVec;

	use super::*;

	fn print_sha256(txt: &'static str, data: &[GFSymbol]) {
		use sha2::Digest;
		let data = unsafe { ::std::slice::from_raw_parts(data.as_ptr() as *const u8, data.len() * 2) };

		let mut digest = sha2::Sha256::new();
		digest.update(data);
		println!("sha256(rs|{}):", txt);
		for byte in digest.finalize().into_iter() {
			print!("{:02x}", byte);
		}
		println!("")
	}

	/// Generate a random index
	fn rand_gf_element() -> GFSymbol {
		use rand::distributions::{Distribution, Uniform};
		use rand::thread_rng;

		let mut rng = thread_rng();
		let uni = Uniform::<GFSymbol>::new_inclusive(0, MODULO);
		uni.sample(&mut rng)
	}

	#[test]
	fn base_2_powers_of_2() {
		assert!(!is_power_of_2(0));
		for i in 0..20 {
			assert!(is_power_of_2(1 << i));
		}
		for i in 0..20 {
			assert!(!is_power_of_2(7 << i));
		}
		let mut f = 3;
		for i in 0..20 {
			f *= 7;
			assert!(!is_power_of_2(f));
		}
		assert_eq!(is_power_of_2(3), false);
	}

	#[test]
	fn base_2_upper_bound() {
		for i in 1_usize..=1024 {
			let upper = next_higher_power_of_2(i);
			if is_power_of_2(i) {
				assert_eq!(upper, i);
			} else {
				assert!(upper > i);
			}
		}
	}

	#[test]
	fn k_n_construction() {
		for validator_count in 3_usize..=4096 {
			let ReedSolomon { n, k } = ReedSolomon::new(validator_count);

			assert!(validator_count <= n);
			assert!(validator_count / 3 >= k);
			assert!(validator_count >= k * 3);
		}
	}

	#[test]
	fn flt_back_and_forth() {
		const N: usize = 128;

		let mut data = (0..N).into_iter().map(|_x| rand_gf_element()).collect::<Vec<GFSymbol>>();
		let expected = data.clone();

		fft_in_novel_poly_basis(&mut data, N, N / 4);

		// make sure something is done
		assert!(data.iter().zip(expected.iter()).filter(|(a, b)| { a != b }).count() > 0);

		inverse_fft_in_novel_poly_basis(&mut data, N, N / 4);

		itertools::assert_equal(data, expected);
	}

	#[test]
	fn sub_encode_decode() {
		setup();
		let mut rng = rand::thread_rng();

		const N: usize = 32;
		const K: usize = 4;

		let mut data = [0u8; K * 2];
		rng.fill_bytes(&mut data[..]);

		let codewords = encode_sub(&data, N, K);
		let mut codewords = codewords.into_iter().map(|x| Some(x)).collect::<Vec<_>>();
		assert_eq!(codewords.len(), N);
		codewords[0] = None;
		codewords[1] = None;
		codewords[2] = None;
		codewords[N - 3] = None;
		codewords[N - 2] = None;
		codewords[N - 1] = None;

		let reconstructed = reconstruct_sub(&codewords, N, K).unwrap();
		itertools::assert_equal(data.iter(), reconstructed.iter().take(K * 2));
	}

	#[test]
	fn flt_rountrip_small() {
		const N: usize = 16;
		const EXPECTED: [GFSymbol; N] = [1, 2, 3, 5, 8, 13, 21, 44, 65, 0, 0xFFFF, 2, 3, 5, 7, 11];

		let mut data = EXPECTED.clone();

		fft_in_novel_poly_basis(&mut data, N, N / 4);

		println!("novel basis(rust):");
		data.iter().for_each(|sym| {
			print!(" {:04X}", sym);
		});
		println!("");

		inverse_fft_in_novel_poly_basis(&mut data, N, N / 4);
		itertools::assert_equal(data.iter(), EXPECTED.iter());
	}

	#[test]
	fn ported_c_test() {
		const N: usize = 256;
		const K: usize = 8;

		setup();

		//-----------Generating message----------
		//message array
		let mut data: [GFSymbol; N] = [0; N];

		for i in 0..K {
			//filled with random numbers
			data[i] = (i * i % MODULO as usize) as u16;
			// data[i] = rand_gf_element();
		}

		assert_eq!(data.len(), N);

		println!("Message(Last n-k are zeros): ");
		for i in 0..K {
			print!("{:04x} ", data[i]);
		}
		println!("");
		print_sha256("data", &data[..]);

		//---------encoding----------
		let mut codeword = [0_u16; N];

		if K + K > N && false {
			let (data_till_t, data_skip_t) = data.split_at_mut(N - K);
			encode_high(data_skip_t, K, data_till_t, &mut codeword[..], N);
		} else {
			encode_low(&data[..], K, &mut codeword[..], N);
		}

		// println!("Codeword:");
		// for i in K..(K+100) {
		// print!("{:04x} ", codeword[i]);
		// }
		// println!("");

		print_sha256("encoded", &codeword);

		//--------erasure simulation---------

		// Array indicating erasures
		let mut erasure = [false; N];

		let erasures_iv = if false {
			// erase random `(N-K)` codewords
			let mut rng = rand::thread_rng();
			let erasures_iv: IndexVec = rand::seq::index::sample(&mut rng, N, N - K);

			erasures_iv
		} else {
			IndexVec::from((0..(N - K)).into_iter().collect::<Vec<usize>>())
		};
		assert_eq!(erasures_iv.len(), N - K);

		for i in erasures_iv {
			//erasure codeword symbols
			erasure[i] = true;
			codeword[i] = 0 as GFSymbol;
		}

		print_sha256("erased", &codeword);

		//---------Erasure decoding----------------
		let mut log_walsh2: [GFSymbol; FIELD_SIZE] = [0_u16; FIELD_SIZE];

		eval_error_polynomial(&erasure[..], &mut log_walsh2[..], FIELD_SIZE);

		print_sha256("log_walsh2", &log_walsh2);

		decode_main(&mut codeword[..], K, &erasure[..], &log_walsh2[..], N);

		print_sha256("decoded", &codeword[0..K]);

		println!("Decoded result:");
		for i in 0..N {
			// the data word plus a few more
			print!("{:04x} ", codeword[i]);
		}
		println!("");

		for i in 0..K {
			//Check the correctness of the result
			if data[i] != codeword[i] {
				println!("🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍");
				panic!("Decoding ERROR! value at [{}] should={:04x} vs is={:04x}", i, data[i], codeword[i]);
			}
		}
		println!(
			r#">>>>>>>>> 🎉🎉🎉🎉
>>>>>>>>> > Decoding is **SUCCESS** ful! 🎈
>>>>>>>>>"#
		);
	}
}
