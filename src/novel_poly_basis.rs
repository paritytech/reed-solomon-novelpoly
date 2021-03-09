// Encoding/erasure decoding for Reed-Solomon codes over binary extension fields
//
// Derived impl of `RSAErasureCode.c`.
//
// Lin, Han and Chung, "Novel Polynomial Basis and Its Application to Reed-Solomon Erasure Codes," FOCS14.
// (http://arxiv.org/abs/1404.3458)

#![allow(dead_code)]

use super::*;

use super::f2e16::*;


use std::slice::from_raw_parts;



//-----Used in decoding procedure-------
//twisted factors used in FFT
static mut SKEW_FACTOR: [GFSymbol; ONEMASK as usize] = [0_u16; ONEMASK as usize];

//factors used in formal derivative
static mut B: [GFSymbol; FIELD_SIZE >> 1] = [0_u16; FIELD_SIZE >> 1];




//return a*EXP_TABLE[b] over GF(2^r)
pub fn mul_table(a: GFSymbol, b: GFSymbol) -> GFSymbol {
    if a != 0_u16 {
        unsafe {
            let ab_log = (LOG_TABLE[a as usize] as u32) + b as u32;
            let offset = (ab_log & ONEMASK as u32) + (ab_log >> FIELD_BITS);
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
pub fn inverse_afft_in_novel_poly_basis(data: &mut [GFSymbol], size: usize, index: usize) {
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
			let skew = Multiplier(unsafe { SKEW_FACTOR[j + index - 1] });
			// It's reasonale to skip the loop if skew is zero, but doing so with
			// all bits set requires justification.	 (TODO)
			if skew.0 != ONEMASK {
				// Again loop on line 3, except skew should depend upon i aka j in Algorithm 2 (TODO)
				for i in (j - depart_no)..j {
					// Line 5, justified by (35) page 6288, but
					// adding depart_no acts like the r+2^i superscript.
					data[i] ^= Additive(data[i + depart_no]).mul(skew).0;
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
pub fn afft_in_novel_poly_basis(data: &mut [GFSymbol], size: usize, index: usize) {
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
			let skew = Multiplier(unsafe { SKEW_FACTOR[j + index - 1] });
			// It's reasonale to skip the loop if skew is zero, but doing so with
			// all bits set requires justification.	 (TODO)
			if skew.0 != ONEMASK {
				// Loop on line 5, except skew should depend upon i aka j in Algorithm 1 (TODO)
				for i in (j - depart_no)..j {
					// Line 6, explained by (28) page 6287, but
					// adding depart_no acts like the r+2^i superscript.
					data[i] ^= Additive(data[i + depart_no]).mul(skew).0;
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


//initialize SKEW_FACTOR and B
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
		// let idx = base[m].mul( (base[m] ^ 1_u16).to_multiplier() );
		base[m] = ONEMASK - LOG_TABLE[idx as usize];

		for i in (m + 1)..(FIELD_BITS - 1) {
			let b = LOG_TABLE[(base[i] as u16 ^ 1_u16) as usize] as u32 + base[m] as u32;
			// let b = (base[i] as u16 ^ 1_u16).to_multiplier().to_wide() as u32 + base[m] as u32;
			let b = b % ONEMASK as u32;
			base[i] = mul_table(base[i], b as u16);
		}
	}
	for i in 0..(ONEMASK as usize) {
		SKEW_FACTOR[i] = LOG_TABLE[SKEW_FACTOR[i] as usize];
	}

	base[0] = ONEMASK - base[0];
	for i in 1..(FIELD_BITS - 1) {
		base[i] = ((ONEMASK as u32 - base[i] as u32 + base[i - 1] as u32) % ONEMASK as u32) as GFSymbol;
	}

	B[0] = 0;
	for i in 0..(FIELD_BITS - 1) {
		let depart = 1 << i;
		for j in 0..depart {
			B[j + depart] = ((B[j] as u32 + base[i] as u32) % ONEMASK as u32) as GFSymbol;
		}
	}
}

/// Setup both decoder and encoder.
pub fn setup() {
	use std::sync::Once;

	static SETUP: Once = Once::new();

	SETUP.call_once(|| unsafe {
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

	inverse_afft_in_novel_poly_basis(codeword_first_k, k, 0);

	// the first codeword is now the basis for the remaining transforms
	// denoted `M_topdash`

	for shift in (k..n).into_iter().step_by(k) {
		let codeword_at_shift = &mut codeword_skip_first_k[(shift - k)..shift];
		// copy `M_topdash` to the position we are currently at, the n transform
		mem_cpy(codeword_at_shift, codeword_first_k);
		afft_in_novel_poly_basis(codeword_at_shift, k, shift);
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

		inverse_afft_in_novel_poly_basis(mem, t, i);
		for j in 0..t {
			parity[j] ^= mem[j];
		}
		i += t;
	}
	afft_in_novel_poly_basis(parity, t, 0);
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
		log_walsh2[i] = (tmp % ONEMASK as u32) as GFSymbol;
	}
	walsh(log_walsh2, FIELD_SIZE);
	for i in 0..z {
		if erasure[i] {
			log_walsh2[i] = ONEMASK - log_walsh2[i];
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

	inverse_afft_in_novel_poly_basis(codeword, n, 0);

	// formal derivative

	#[cfg(feature = "b_not_one")]
	for i in (0..n).into_iter().step_by(2) {
		let b = ONEMASK - unsafe { B[i >> 1] };
    	// #[cfg(not(feature = "b_not_one"))]
		// assert_eq!(b, ONEMASK);
		codeword[i] = mul_table(codeword[i], b);
		codeword[i + 1] = mul_table(codeword[i + 1], b);
	}

	formal_derivative(codeword, n);

	#[cfg(feature = "b_not_one")]
	for i in (0..n).into_iter().step_by(2) {
		let b = unsafe { B[i >> 1] };
		codeword[i] = mul_table(codeword[i], b);
		codeword[i + 1] = mul_table(codeword[i + 1], b);
	}

	afft_in_novel_poly_basis(codeword, n, 0);

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

/// Params for the encoder / decoder
/// derived from a target validator count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CodeParams {
	/// total number of message symbols to send
	/// Invariant is a power of base 2
	n: usize,
	/// number of information containing chunks
	/// Invariant is a power of base 2, `k < n`
	k: usize,
	/// Avoid copying unnecessary chunks.
	validator_count: usize,
}

impl CodeParams {
	/// Create a new reed solomon erasure encoding wrapper
	pub fn derive_from_validator_count(validator_count: usize) -> Result<Self> {
		if validator_count < 2 {
			return Err(Error::ValidatorCountTooLow(validator_count))
		}
		// we need to be able to reconstruct from 1/3 - eps
		let k = std::cmp::max(validator_count / 3, 1); // for the odd case of 2 validators
		let k = next_lower_power_of_2(k);
		let n = next_higher_power_of_2(validator_count);
		if n > FIELD_SIZE as usize {
			return Err(Error::ValidatorCountTooLow(validator_count))
		}
		Ok(Self { n, k, validator_count })
	}

	// make a reed-solomon instance.
	pub fn make_encoder(&self) -> ReedSolomon {
		ReedSolomon::new(self.n, self.k, self.validator_count)
			.expect("this struct is not created with invalid shard number; qed")
	}
}

pub fn encode(bytes: &[u8], validator_count: usize) -> Result<Vec<WrappedShard>> {
	let params = CodeParams::derive_from_validator_count(validator_count)?;

	let rs = params.make_encoder();
	rs.encode(bytes)
}

/// each shard contains one symbol of one run of erasure coding
pub fn reconstruct(received_shards: Vec<Option<WrappedShard>>, validator_count: usize) -> Result<Vec<u8>> {
	let params = CodeParams::derive_from_validator_count(validator_count)?;

	let rs = params.make_encoder();
	rs.reconstruct(received_shards)
}

pub struct ReedSolomon {
	n: usize,
	k: usize,
	validator_count: usize,
}

impl ReedSolomon {
	/// Returns the size per shard in bytes
	pub fn shard_len(&self, payload_size: usize) -> usize {
		let payload_symbols = (payload_size + 1) / 2;
		let shard_symbols_ceil = (payload_symbols + self.k - 1) / self.k;
		let shard_bytes = shard_symbols_ceil * 2;
		shard_bytes
	}

	pub fn new(n: usize, k: usize, validator_count: usize) -> Result<Self> {
		setup();
		if !is_power_of_2(n) && !is_power_of_2(k) {
			Err(Error::ParamterMustBePowerOf2 { n, k })
		} else {
			Ok(Self { validator_count, n, k })
		}
	}

	pub fn encode(&self, bytes: &[u8]) -> Result<Vec<WrappedShard>> {
		if bytes.is_empty() {
			return Err(Error::PayloadSizeIsZero);
		}

		// setup the shards, n is likely _larger_, so use the truely required number of shards

		// required shard length in bytes, rounded to full symbols
		let shard_len = self.shard_len(bytes.len());
		assert!(shard_len > 0);
		// collect all sub encoding runs

		let validator_count = self.validator_count;
		let k2 = self.k * 2;
		// prepare one wrapped shard per validator
		let mut shards = vec![
			WrappedShard::new({
				let mut v = Vec::<u8>::with_capacity(shard_len);
				unsafe { v.set_len(shard_len) }
				v
			});
			validator_count
		];

		for (chunk_idx, i) in (0..bytes.len()).into_iter().step_by(k2).enumerate() {
			let end = std::cmp::min(i + k2, bytes.len());
			assert_ne!(i, end);
			let data_piece = &bytes[i..end];
			assert!(!data_piece.is_empty());
			assert!(data_piece.len() <= k2);
			let encoding_run = encode_sub(data_piece, self.n, self.k)?;
			for val_idx in 0..validator_count {
				AsMut::<[[u8; 2]]>::as_mut(&mut shards[val_idx])[chunk_idx] = encoding_run[val_idx].to_be_bytes();
			}
		}

		Ok(shards)
	}

	/// each shard contains one symbol of one run of erasure coding
	pub fn reconstruct(&self, received_shards: Vec<Option<WrappedShard>>) -> Result<Vec<u8>> {
		let gap = self.n.saturating_sub(received_shards.len());
		let received_shards =
			received_shards.into_iter().take(self.n).chain(std::iter::repeat(None).take(gap)).collect::<Vec<_>>();

		assert_eq!(received_shards.len(), self.n);

		// obtain a sample of a shard length and assume that is the truth
		// XXX make sure all shards have equal length
		let shard_len_in_syms = received_shards
			.iter()
			.find_map(|x| {
				x.as_ref().map(|x| {
					let x = AsRef::<[[u8; 2]]>::as_ref(x);
					x.len()
				})
			})
			.unwrap();

		// TODO check shard length is what we'd expect

		let mut existential_count = 0_usize;
		let erasures = received_shards
			.iter()
			.map(|x| x.is_none())
			.inspect(|erased| existential_count += !*erased as usize)
			.collect::<Vec<bool>>();

		if existential_count < self.k {
			return Err(Error::NeedMoreShards { have: existential_count, min: self.k, all: self.n });
		}

		// Evaluate error locator polynomial only once
		let mut error_poly_in_log = [0_u16 as GFSymbol; FIELD_SIZE];
		eval_error_polynomial(&erasures[..], &mut error_poly_in_log[..], FIELD_SIZE);

		let mut acc = Vec::<u8>::with_capacity(shard_len_in_syms * 2 * self.k);
		for i in 0..shard_len_in_syms {
			// take the i-th element of all shards and try to recover
			let decoding_run = received_shards
				.iter()
				.map(|x| {
					x.as_ref().map(|x| {
						let z = AsRef::<[[u8; 2]]>::as_ref(&x)[i];
						u16::from_be_bytes(z)
					})
				})
				.collect::<Vec<Option<GFSymbol>>>();

			assert_eq!(decoding_run.len(), self.n);

			// reconstruct from one set of symbols which was spread over all erasure chunks
			let piece = reconstruct_sub(&decoding_run[..], &erasures, self.n, self.k, &error_poly_in_log).unwrap();
			acc.extend_from_slice(&piece[..]);
		}

		Ok(acc)
	}
}

/// Bytes shall only contain payload data
pub fn encode_sub(bytes: &[u8], n: usize, k: usize) -> Result<Vec<GFSymbol>> {
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

	// update new data bytes with zero padded bytes
	// `l` is now `GF(2^16)` symbols
	let l = data.len();
	assert_eq!(l, n);

	let mut codeword = data.clone();
	assert_eq!(codeword.len(), n);

	encode_low(&data[..], k, &mut codeword[..], n);

	Ok(codeword)
}

pub fn reconstruct_sub(
	codewords: &[Option<GFSymbol>],
	erasures: &[bool],
	n: usize,
	k: usize,
	error_poly: &[GFSymbol; FIELD_SIZE],
) -> Result<Vec<u8>> {
	assert!(is_power_of_2(n), "Algorithm only works for 2^i sizes for N");
	assert!(is_power_of_2(k), "Algorithm only works for 2^i sizes for K");
	assert_eq!(codewords.len(), n);
	assert!(k <= n / 2);

	// the first k suffice for the original k message codewords
	let recover_up_to = k; // n;

	// The recovered _payload_ chunks AND parity chunks
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

	// the first k would suffice for the original k message codewords
	let recover_up_to = k;

	//---------Erasure decoding----------------

	decode_main(&mut codeword[..], recover_up_to, &erasures[..], &error_poly[..], n);

	for idx in 0..recover_up_to {
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
	Ok(recovered[0..k2].to_vec())
}

#[cfg(test)]
mod test {
	use rand::distributions::Uniform;
	use rand::seq::index::IndexVec;
	use assert_matches::assert_matches;

	use super::*;

	// If this passes then you do not require the b_not_one feature
	fn b_is_one() {
		let n = FIELD_SIZE >> 1;
		// Everything
		// for i in (0..n) {
		// Just like in decode_main
		for i in (0..n).into_iter().step_by(2) {
			let b = ONEMASK - unsafe { B[i >> 1] };
			assert_eq!(b, ONEMASK);
		}
	}

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
		let mut rng = thread_rng();
		let uni = Uniform::<GFSymbol>::new_inclusive(0, ONEMASK);
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
		for _i in 0..20 {
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
		// skip the two, it's a special case
		for validator_count in 3_usize..=8200 {
			let CodeParams { n, k, .. } = CodeParams::derive_from_validator_count(validator_count).unwrap();
			assert!(validator_count <= n, "vc={} <= n={} violated", validator_count, n);
			assert!(validator_count / 3 >= k, "vc={} / 3 >= k={} violated", validator_count, k);
			assert!(validator_count >= k * 3, "vc={} <= k={} *3  violated", validator_count, k);
		}
	}

	#[test]
	fn flt_back_and_forth() {
		const N: usize = 128;

		let mut data = (0..N).into_iter().map(|_x| rand_gf_element()).collect::<Vec<GFSymbol>>();
		let expected = data.clone();

		afft_in_novel_poly_basis(&mut data, N, N / 4);

		// make sure something is done
		assert!(data.iter().zip(expected.iter()).filter(|(a, b)| { a != b }).count() > 0);

		inverse_afft_in_novel_poly_basis(&mut data, N, N / 4);

		itertools::assert_equal(data, expected);
	}

	#[test]
	fn sub_encode_decode() -> Result<()> {
		setup();
		let mut rng = rand::thread_rng();

		const N: usize = 32;
		const K: usize = 4;

		const K2: usize = K * 2;
		let mut data = [0u8; K2];
		rng.fill_bytes(&mut data[..]);

		let codewords = encode_sub(&data, N, K)?;
		let mut codewords = codewords.into_iter().map(|x| Some(x)).collect::<Vec<_>>();
		assert_eq!(codewords.len(), N);
		codewords[0] = None;
		codewords[1] = None;
		codewords[2] = None;
		codewords[N - 3] = None;
		codewords[N - 2] = None;
		codewords[N - 1] = None;

		let erasures = codewords.iter().map(|x| x.is_none()).collect::<Vec<bool>>();

		// Evaluate error locator polynomial only once
		let mut error_poly_in_log = [0_u16 as GFSymbol; FIELD_SIZE];
		eval_error_polynomial(&erasures[..], &mut error_poly_in_log[..], FIELD_SIZE);

		let reconstructed = reconstruct_sub(&codewords[..], &erasures[..], N, K, &error_poly_in_log)?;
		itertools::assert_equal(data.iter(), reconstructed.iter().take(K2));
		Ok(())
	}

	fn deterministic_drop_shards<T: Sized, G: rand::SeedableRng + rand::Rng>(
		codewords: &mut [Option<T>],
		n: usize,
		k: usize,
		_rng: &mut G,
	) -> IndexVec {
		let l = codewords.len();
		let mut v = Vec::with_capacity(n - k);
		// k is a power of 2
		let half = (n - k) >> 1;
		for i in 0..half {
			codewords[i] = None;
			v.push(i);
		}
		// if the codewords is shorter than n
		// the remaining ones were
		// already dropped implicitly
		for i in n - half..n {
			if i < l {
				codewords[i] = None;
				v.push(i);
			}
		}
		IndexVec::from(v)
	}

	fn deterministic_drop_shards_clone<T: Sized + Clone>(
		codewords: &[T],
		n: usize,
		k: usize,
	) -> (Vec<Option<T>>, IndexVec) {
		let mut rng = SmallRng::from_seed(crate::SMALL_RNG_SEED);
		let mut codewords = codewords.into_iter().map(|x| Some(x.clone())).collect::<Vec<Option<T>>>();
		let idx = deterministic_drop_shards::<T, SmallRng>(&mut codewords, n, k, &mut rng);
		assert!(idx.len() <= n - k);
		(codewords, idx)
	}

	// for shards of length 1
	fn wrapped_shard_len1_as_gf_sym(w: &WrappedShard) -> GFSymbol {
		let val = AsRef::<[[u8; 2]]>::as_ref(w)[0];
		u16::from_be_bytes(val)
	}

	#[test]
	fn sub_eq_big_for_small_messages() {
		const N_VALIDATORS: usize = 128;
		const N: usize = N_VALIDATORS;
		const K: usize = 32;

		setup();

		const K2: usize = K * 2;

		// assure the derived sizes match
		let rs = CodeParams::derive_from_validator_count(N_VALIDATORS).unwrap();
		assert_eq!(rs.n, N);
		assert_eq!(rs.k, K);

		// create random predictable bytes
		// and create a message that results in 1 GF element symbols
		// per validator
		let data = {
			let mut rng = SmallRng::from_seed(crate::SMALL_RNG_SEED);
			let mut data = [0u8; K2];
			rng.fill_bytes(&mut data[..]);
			data
		};

		let mut codewords = encode(&data, rs.n).unwrap();
		let mut codewords_sub = encode_sub(&data, N, K).unwrap();

		itertools::assert_equal(codewords.iter().map(wrapped_shard_len1_as_gf_sym), codewords_sub.iter().copied());

		let (codewords, _) = deterministic_drop_shards_clone(&mut codewords, N, K);
		let (codewords_sub, _) = deterministic_drop_shards_clone(&mut codewords_sub, N, K);

		itertools::assert_equal(
			codewords.iter().map(|w| w.as_ref().map(wrapped_shard_len1_as_gf_sym)),
			codewords_sub.iter().copied(),
		);

		let erasures = codewords.iter().map(|x| x.is_none()).collect::<Vec<bool>>();

		// Evaluate error locator polynomial only once
		let mut error_poly_in_log = [0_u16 as GFSymbol; FIELD_SIZE];
		eval_error_polynomial(&erasures[..], &mut error_poly_in_log[..], FIELD_SIZE);

		let reconstructed_sub = reconstruct_sub(&codewords_sub[..], &erasures[..], N, K, &error_poly_in_log).unwrap();
		let reconstructed = reconstruct(codewords, rs.n).unwrap();
		itertools::assert_equal(reconstructed.iter().take(K2), reconstructed_sub.iter().take(K2));
		itertools::assert_equal(reconstructed.iter().take(K2), data.iter());
		itertools::assert_equal(reconstructed_sub.iter().take(K2), data.iter());
	}

	#[test]
	fn roundtrip_for_large_messages() -> Result<()> {
		const N_VALIDATORS: usize = 2000;
		const N: usize = 2048;
		const K: usize = 512;

		setup();

		const K2: usize = K * 2;

		// assure the derived sizes match
		let rs = CodeParams::derive_from_validator_count(N_VALIDATORS).unwrap();
		assert_eq!(rs.n, N);
		assert_eq!(rs.k, K);

		// make sure each shard is more than one byte to
		// test the shard size
		// in GF symbols
		let shard_length: usize = 23;

		let payload = &crate::BYTES[0..K2 * shard_length];
		// let payload = &crate::BYTES[..];

		let mut shards = encode(payload, N_VALIDATORS).unwrap();

		// for (idx, shard) in shards.iter().enumerate() {
		//	let sl = AsRef::<[[u8; 2]]>::as_ref(&shard).len();
		//	assert_eq!(shard_length, sl, "Shard #{} has an unxpected length {} (expected: {})", idx, sl, shard_length);
		// }
		let (received_shards, dropped_indices) = deterministic_drop_shards_clone(&mut shards, rs.n, rs.k);

		let reconstructed_payload = reconstruct(received_shards, N_VALIDATORS).unwrap();

		assert_recovery(payload, &reconstructed_payload, dropped_indices);

		// verify integrity with criterion tests
		roundtrip_w_drop_closure::<_, _, _, SmallRng>(
			encode,
			reconstruct,
			payload,
			N_VALIDATORS,
			deterministic_drop_shards::<WrappedShard, SmallRng>,
		)?;

		roundtrip_w_drop_closure::<_, _, _, SmallRng>(
			encode,
			reconstruct,
			payload,
			N_VALIDATORS,
			crate::drop_random_max,
		)?;

		Ok(())
	}

	macro_rules! simplicissimus {
		($name:ident: validators: $validator_count:literal, payload: $payload_size:literal; $matchmaker:pat) => {
			simplicissimus!($name: validators: $validator_count, payload: $payload_size; $matchmaker => {});
		};
		($name:ident: validators: $validator_count:literal, payload: $payload_size:literal) => {
			simplicissimus!($name: validators: $validator_count, payload: $payload_size; Ok(x) => { let _ = x; });
		};
		($name:ident: validators: $validator_count:literal, payload: $payload_size:literal; $matchmaker:pat => $assertive:expr) => {
			#[test]
			fn $name () {
				let res = roundtrip_w_drop_closure::<_,_,_,SmallRng>(
					encode,
					reconstruct,
					&BYTES[0..$payload_size], $validator_count,
					 deterministic_drop_shards::<WrappedShard, SmallRng>);
				assert_matches::assert_matches!(res, $matchmaker => {
					$assertive
				});
			}
		};
	}

	simplicissimus!(case_0: validators: 2003, payload: 0; Err(Error::PayloadSizeIsZero));

	// Roughly one GFSymbol per validator payload
	simplicissimus!(case_1: validators: 10, payload: 16);

	// Unit payload, but mayn validators
	simplicissimus!(case_2: validators: 100, payload: 1);

	// Common case of way ore payload than validators
	simplicissimus!(case_3: validators: 4, payload: 100);

	// Way more validators than payload bytes
	simplicissimus!(case_4: validators: 2003, payload: 17);

	#[test]
	fn flt_roundtrip_small() {
		const N: usize = 16;
		const EXPECTED: [GFSymbol; N] = [1, 2, 3, 5, 8, 13, 21, 44, 65, 0, 0xFFFF, 2, 3, 5, 7, 11];

		let mut data = EXPECTED.clone();

		afft_in_novel_poly_basis(&mut data, N, N / 4);

		println!("novel basis(rust):");
		data.iter().for_each(|sym| {
			print!(" {:04X}", sym);
		});
		println!("");

		inverse_afft_in_novel_poly_basis(&mut data, N, N / 4);
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
			data[i] = (i * i % ONEMASK as usize) as u16;
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


	#[test]
	fn test_code_params() {
		assert_matches!(CodeParams::derive_from_validator_count(0), Err(_));

		assert_matches!(CodeParams::derive_from_validator_count(1), Err(_));

		assert_eq!(CodeParams::derive_from_validator_count(2), Ok(CodeParams {
			n: 2,
			k: 1,
			validator_count: 2,
		}));

		assert_eq!(CodeParams::derive_from_validator_count(3), Ok(CodeParams {
			n: 4,
			k: 1,
			validator_count: 3,
		}));

		assert_eq!(CodeParams::derive_from_validator_count(4), Ok(CodeParams {
			n: 4,
			k: 1,
			validator_count: 4,
		}));

		assert_eq!(CodeParams::derive_from_validator_count(100), Ok(CodeParams {
			n: 128,
			k: 32,
			validator_count: 100,
		}));
	}

	#[test]
	fn shard_len_is_reasonable() {
		let rs = CodeParams {
			n: 16,
			k: 4,
			validator_count: 5,
		}.make_encoder();

		// since n must be a power of 2
		// the chunk sizes becomes slightly larger
		// than strictly necessary
		assert_eq!(rs.shard_len(100), 26);
		assert_eq!(rs.shard_len(99), 26);

		// see if it rounds up to 2.
		assert_eq!(rs.shard_len(95), 24);
		assert_eq!(rs.shard_len(94), 24);

		assert_eq!(rs.shard_len(90), 24);

		// needs 3 bytes to fit, rounded up to next even number.
		assert_eq!(rs.shard_len(19), 6);
	}
}
