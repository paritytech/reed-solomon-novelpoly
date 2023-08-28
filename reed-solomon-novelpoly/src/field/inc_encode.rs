#[inline(always)]
pub fn encode_low(data: &[Additive], k: usize, codeword: &mut [Additive], n: usize) {
	if k >= 16 && k % 8 == 0 && n % 8 == 0 && (n-k) % 8 == 0 {
		encode_low_faster8_adaptor(data, k, codeword, n);		
	} else {
		encode_low_plain(data, k, codeword, n);
	}
}


// Encoding alg for k/n < 0.5: message is a power of two
pub fn encode_low_plain(data: &[Additive], k: usize, codeword: &mut [Additive], n: usize) {
	assert!(k + k <= n);
	assert_eq!(codeword.len(), n);
	assert_eq!(data.len(), n);

	assert!(is_power_of_2(n));
	assert!(is_power_of_2(k));

	// k | n is guaranteed
	assert_eq!((n / k) * k, n);

	// move the data to the codeword
    codeword.copy_from_slice(data);

	// split after the first k
	let (codeword_first_k, codeword_skip_first_k) = codeword.split_at_mut(k);

    inverse_afft(codeword_first_k, k, 0);

	// dbg!(&codeword_first_k);
	// the first codeword is now the basis for the remaining transforms
	// denoted `M_topdash`

	for shift in (k..n).into_iter().step_by(k) {
		let codeword_at_shift = &mut codeword_skip_first_k[(shift - k)..shift];
		// copy `M_topdash` to the position we are currently at, the n transform
		codeword_at_shift.copy_from_slice(codeword_first_k);
		// dbg!(&codeword_at_shift);
		afft(codeword_at_shift, dbg!(k), dbg!(shift));
		// let post = &codeword_at_shift;
		// dbg!(post);
	}

	// restore `M` from the derived ones
	(&mut codeword[0..k]).copy_from_slice(&data[0..k]);
}

pub fn encode_low_faster8_adaptor(data: &[Additive], k: usize, codeword: &mut [Additive], n: usize) {
	assert_eq!(n % Additive8x::LANE, 0);
	let mut codeword8x = vec![Additive8x::zero(); n / Additive8x::LANE];
	convert_to_faster8(&data[..k], &mut codeword8x[..]);
	let data8x = codeword8x.clone();
	encode_low_faster8(&data8x[..], k, &mut codeword8x[..], n);
	convert_from_faster8(&codeword8x[..], &mut codeword[..]);	
}

pub fn encode_low_faster8(data8x: &[Additive8x], k: usize, codeword8x: &mut [Additive8x], n: usize) {
	assert!(k + k <= n);
	assert_eq!(codeword8x.len(), n/Additive8x::LANE);
	assert_eq!(data8x.len(), n/Additive8x::LANE); // FIXME, we never read data beyond 0..k

	assert!(is_power_of_2(n));
	assert!(is_power_of_2(k));

	assert_eq!(k % Additive8x::LANE, 0);
	assert_eq!((n-k) % Additive8x::LANE, 0);

	// k | n is guaranteed
	assert_eq!((n / k) * k, n);

	
	let k_8x = k / Additive8x::LANE;
	let n_8x = n / Additive8x::LANE;
	
	// move the data to the codeword
	// split after the first k
	let (codeword8x_first_k, codeword8x_skip_first_k) = codeword8x.split_at_mut(k_8x);

	assert!((k >> 1) >= Additive8x::LANE);
    inverse_afft_faster8(codeword8x_first_k, k, 0);

	// dbg!(&codeword8x_first_k);

	// the first codeword is now the basis for the remaining transforms
	// denoted `M_topdash`

	
	for shift_8x in (k_8x..n_8x).into_iter().step_by(k_8x) {
		let codeword8x_at_shift = &mut codeword8x_skip_first_k[(shift_8x - k_8x)..][..k_8x];
		// copy `M_topdash` to the position we are currently at, the n transform
		codeword8x_at_shift.copy_from_slice(codeword8x_first_k);
		// dbg!(&codeword8x_at_shift);
		afft_faster8(codeword8x_at_shift, dbg!(k), shift_8x * Additive8x::LANE);
		// let post = &codeword8x_at_shift;
		// dbg!(post);
	}

	// restore `M` from the derived ones
	(&mut codeword8x[0..k_8x]).copy_from_slice(&data8x[0..k_8x]);
}



//data: message array. parity: parity array. mem: buffer(size>= n-k)
//Encoding alg for k/n>0.5: parity is a power of two.
#[inline(always)]
pub fn encode_high(data: &[Additive], k: usize, parity: &mut [Additive], mem: &mut [Additive], n: usize) {
	if (n-k) % Additive8x::LANE == 0 && n % Additive8x::LANE == 0 && k % Additive8x::LANE == 0 {
		encode_high_faster8_adapter(data, k, parity, mem, n);
	} else {
		encode_high_plain(data, k, parity, mem, n);
	}
}

//data: message array. parity: parity array. mem: buffer(size>= n-k)
//Encoding alg for k/n>0.5: parity is a power of two.
pub fn encode_high_plain(data: &[Additive], k: usize, parity: &mut [Additive], mem: &mut [Additive], n: usize) {
	let t: usize = n - k;
	
	// mem_zero(&mut parity[0..t]);
	for i in 0..t {
		parity[i] = Additive(0);
	}

	let mut i = t;
	while i < n {
		(&mut mem[..t]).copy_from_slice(&data[(i - t)..t]);

		inverse_afft(mem, t, i);
		for j in 0..t {
			parity[j] ^= mem[j];
		}
		i += t;
	}
	afft(parity, t, 0);
}

/// Avoid using this, it incures quite some cost
pub fn encode_high_faster8_adapter(data: &[Additive], k: usize, parity: &mut [Additive], mem: &mut [Additive], n: usize) {
	let mut mem8 = vec![Additive8x::zero(); mem.len()/Additive8x::LANE];
	let data8 = Vec::from_iter(data.iter().step_by(Additive8x::LANE).enumerate().map(|(piece_idx, _offset)| {
		Additive8x::load(&data[(piece_idx*Additive8x::LANE)..][..Additive8x::LANE])
	}));
	let mut parity8 = Vec::from_iter(data.iter().step_by(Additive8x::LANE).enumerate().map(|(piece_idx, _offset)| {
		Additive8x::load(&parity[(piece_idx*Additive8x::LANE)..][..Additive8x::LANE])
	}));
	
	encode_high_faster8(&data8, k, &mut parity8, &mut mem8, n);
	
	for (i, mem8) in mem8.iter().enumerate() {
		mem8.copy_to_slice(&mut mem[(Additive8x::LANE * i)..][..8]);
	}
}

pub fn encode_high_faster8(data: &[Additive8x], k: usize, parity: &mut [Additive8x], mem: &mut [Additive8x], n: usize) {
	let t: usize = n - k;
	assert!(t >= 8);
	assert_eq!(t % 8, 0);

	let t8s = t >> 3;
	for i in 0..t8s {
		parity[i] = Additive8x::zero();
	}

	let mut i = t8s;
	while i < n {
		(&mut mem[..t8s]).copy_from_slice(&data[(i - t8s)..t]);

		inverse_afft_faster8(mem, t8s, i);
		for j in 0..t8s {
			parity[j] ^= mem[j];
		}
		i += t8s;
	}
	afft_faster8(parity, t8s, 0);
}

pub fn encode_sub(bytes: &[u8], n: usize, k: usize) -> Result<Vec<Additive>> {
	if (k % Additive8x::LANE) == 0 && (k >> 1) >= Additive8x::LANE {
		encode_sub_faster8(bytes, n, k)
	} else {
		encode_sub_plain(bytes, n, k)
	}
}

/// Bytes shall only contain payload data
pub fn encode_sub_plain(bytes: &[u8], n: usize, k: usize) -> Result<Vec<Additive>> {
	assert!(is_power_of_2(n), "Algorithm only works for 2^i sizes for N");
	assert!(is_power_of_2(k), "Algorithm only works for 2^i sizes for K");
	assert!(bytes.len() <= k << 1);
	assert!(k <= n / 2);

	// must be power of 2
	let bytes_len = bytes.len();
	let upper_len = if is_power_of_2(bytes_len) {
		bytes_len
	} else {
		let loglen = log2(bytes_len);
		let upper_len = 1 << loglen;
		let upper_len = if upper_len >= bytes_len { upper_len } else { upper_len << 1 };
		upper_len
	};
	assert!(is_power_of_2(upper_len));
	assert!(upper_len >= bytes_len);

    // tuple_windows are only used here
    use itertools::Itertools;

	// pad the incoming bytes with trailing 0s
	// so we get a buffer of size `N` in `GF` symbols
	let zero_bytes_to_add = n * 2 - bytes_len;
	let elm_data = Vec::<Additive>::from_iter(bytes
		.into_iter()
		.copied()
		.chain(std::iter::repeat(0u8).take(zero_bytes_to_add))
		.tuple_windows()
		.step_by(2)
		.map(|(a, b)| Additive(Elt::from_be_bytes([a, b]))));

	// update new data bytes with zero padded bytes
	// `l` is now `GF(2^16)` symbols
	let elm_len = elm_data.len();
	assert_eq!(elm_len, n);

	let mut codeword = elm_data.clone();
	assert_eq!(codeword.len(), n);

	encode_low_plain(&elm_data[..], k, &mut codeword[..], n);

	Ok(codeword)
}


/// Bytes shall only contain payload data
pub fn encode_sub_faster8(bytes: &[u8], n: usize, k: usize) -> Result<Vec<Additive>> {
	assert!(is_power_of_2(n), "Algorithm only works for 2^i sizes for N");
	assert!(is_power_of_2(k), "Algorithm only works for 2^i sizes for K");
	assert!(bytes.len() <= k << 1);
	assert!(k <= n / 2);
	assert_eq!(k % Additive8x::LANE, 0);
	assert!((k >> 1) >= Additive8x::LANE);

	// must be power of 2
	let bytes_len = bytes.len();
	let upper_len = if is_power_of_2(bytes_len) {
		bytes_len
	} else {
		let loglen = log2(std::cmp::max(Additive8x::LANE, bytes_len));
		let upper_len = 1 << loglen;
		let upper_len = if upper_len >= bytes_len { upper_len } else { upper_len << 1 };
		upper_len
	};
	assert!(is_power_of_2(upper_len));
	assert!(upper_len >= bytes_len);

    // tuple_windows are only used here
    use itertools::Itertools;

	// pad the incoming bytes with trailing 0s
	// so we get a buffer of size `N` in `GF` symbols
	let zero_bytes_to_add = n * 2 - bytes_len;
	let data = Vec::<Additive>::from_iter(
		bytes
		.into_iter()
		.copied()
		.chain(std::iter::repeat(0u8).take(zero_bytes_to_add))
		.tuple_windows()
		.step_by(2)
		.map(|(a,b)| Additive(Elt::from_be_bytes([a, b])))
	);
	let data8x = Vec::<Additive8x>::from_iter(data.chunks(Additive8x::LANE).map(|sliced| Additive8x::load(sliced)));

	// update new data bytes with zero padded bytes
	// `l` is now `GF(2^16)` symbols
	let elm_len = data8x.len();
	assert_eq!(elm_len * Additive8x::LANE, n);

	let mut codeword8x = data8x.clone();
	assert_eq!(codeword8x.len(), elm_len);

	encode_low_faster8(&data8x[..], k, &mut codeword8x[..], n);

	let mut codeword = Vec::<Additive>::with_capacity(Additive8x::LANE * codeword8x.len());
	unsafe {
		codeword.set_len(codeword.capacity());
	}

	codeword8x
		.into_iter()
		.enumerate()
		.for_each(
	|(idx, item8x)| {
		item8x.copy_to_slice(&mut codeword[(idx * Additive8x::LANE)..][..Additive8x::LANE]);
	});
	Ok(codeword)
}

#[cfg(test)]
mod tests_plain_vs_faster8 {
	use super::*;
	
	#[test]
	fn encode_low_output_plain_eq_faster8() {
		// k must be larger, since the afft is only accelerated by lower values
		let n: usize = 32;
		let k: usize = 16;
		let data1 = vec![Additive(0x1234); n];
		let data2 = data1.clone();
		
		let mut parity1 = vec![Additive::zero(); n];
		encode_low_plain(&data1[..], k, &mut parity1[..], n);

		let mut parity2 = vec![Additive::zero(); n];
		encode_low_faster8_adaptor(&data2[..], k, &mut parity2[..], n);

		assert_eq!(parity1, parity2);
	}
	
	
	#[test]
	fn encode_sub_output_plain_eq_faster8() {
		let n = 64;
		let k = 32; // smallest supported size
		let bytes = vec![0x2A_u8; 64];
		let bytes = bytes.as_slice();
		let fe_plain = encode_sub_plain(bytes, n, k).unwrap();
		let fe_faster8: Vec<Additive> = encode_sub_faster8(bytes, n, k).unwrap();
		
		assert_eq!(fe_plain, fe_faster8);
	}

}
