#[inline(always)]
pub fn encode_low(data: &[Additive], k: usize, codeword: &mut [Additive], n: usize) {
	#[cfg(target_feature = "avx")]
	if k >= 16 && k % 8 == 0 && n % 8 == 0 && (n-k) % 8 == 0 {
		encode_low_faster8(data, k, codeword, n);		
	} else {
		encode_low_plain(data, k, codeword, n);
	}
	
	#[cfg(not(target_feature = "avx"))]
	encode_low_plain(data, k, codeword, n);
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

	for shift in (k..n).step_by(k) {
		let codeword_at_shift = &mut codeword_skip_first_k[(shift - k)..shift];
		// copy `M_topdash` to the position we are currently at, the n transform
		codeword_at_shift.copy_from_slice(codeword_first_k);
		// dbg!(&codeword_at_shift);
		afft(codeword_at_shift, k, shift);
		// let post = &codeword_at_shift;
		// dbg!(post);
	}

	// restore `M` from the derived ones
	codeword[0..k].copy_from_slice(&data[0..k]);
}


#[cfg(target_feature = "avx")]
pub fn encode_low_faster8(data: &[Additive], k: usize, codeword: &mut [Additive], n: usize) {
	assert!(k + k <= n);
	assert_eq!(codeword.len(), n);
	assert_eq!(data.len(), n); // FIXME, we never read data beyond 0..k
	
	assert!(is_power_of_2(n));
	assert!(is_power_of_2(k));
	
	assert_eq!(k % Additive8x::LANE, 0);
	assert_eq!(n % Additive8x::LANE, 0);
	assert_eq!((n-k) % Additive8x::LANE, 0);

	// k | n is guaranteed
	assert_eq!((n / k) * k, n);

	
	// move the data to the codeword
	codeword.copy_from_slice(data);

	// split after the first k
	let (codeword_first_k, codeword_skip_first_k) = codeword.split_at_mut(k);

	assert!((k >> 1) >= Additive8x::LANE);
    inverse_afft_faster8(codeword_first_k, k, 0);


	// the first codeword is now the basis for the remaining transforms
	// denoted `M_topdash`

	
	for shift in (k..n).step_by(k) {
		let codeword_at_shift = &mut codeword_skip_first_k[(shift - k)..shift];
		// copy `M_topdash` to the position we are currently at, the n transform
		codeword_at_shift.copy_from_slice(codeword_first_k);

		afft_faster8(codeword_at_shift, k, shift);
		// let post = &codeword8x_at_shift;
	}

	// restore `M` from the derived ones
	
	codeword[0..k].copy_from_slice(&data[0..k]);
}



//data: message array. parity: parity array. mem: buffer(size>= n-k)
//Encoding alg for k/n>0.5: parity is a power of two.
#[inline(always)]
pub fn encode_high(data: &[Additive], k: usize, parity: &mut [Additive], mem: &mut [Additive], n: usize) {
	#[cfg(target_feature = "avx")]
	if (n-k) % Additive8x::LANE == 0 && n % Additive8x::LANE == 0 && k % Additive8x::LANE == 0 {
		encode_high_faster8(data, k, parity, mem, n);
	} else {
		encode_high_plain(data, k, parity, mem, n);
	}
	#[cfg(not(target_feature = "avx"))]
	encode_high_plain(data, k, parity, mem, n);
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
		mem[..t].copy_from_slice(&data[(i - t)..t]);

		inverse_afft(mem, t, i);
		for j in 0..t {
			parity[j] ^= mem[j];
		}
		i += t;
	}
	afft(parity, t, 0);
}

#[cfg(target_feature = "avx")]
pub fn encode_high_faster8_adapter(data: &[Additive], k: usize, parity: &mut [Additive], mem: &mut [Additive], n: usize) {
	encode_high_faster8(&data, k, parity, mem, n);
}

#[cfg(target_feature = "avx")]
pub fn encode_high_faster8(data: &[Additive], k: usize, parity: &mut [Additive], mem: &mut [Additive], n: usize) {
	let t: usize = n - k;
	assert!(t >= 8);
	assert_eq!(t % 8, 0);

	for i in 0..t {
		parity[i] = Additive::zero();
	}

	let mut i = t;
	while i < n {
		mem[..t].copy_from_slice(&data[(i - t)..t]);

		inverse_afft_faster8(mem, t, i);
		for j in 0..t {
			parity[j] ^= mem[j];
		}
		i += t;
	}
	afft_faster8(parity, t, 0);
}

pub fn encode_sub(bytes: &[u8], n: usize, k: usize) -> Result<Vec<Additive>> {
	#[cfg(target_feature = "avx")]
	if (k % Additive8x::LANE) == 0 && (k >> 1) >= Additive8x::LANE {
		encode_sub_faster8(bytes, n, k)
	} else {
		encode_sub_plain(bytes, n, k)
	}
	#[cfg(not(target_feature = "avx"))]
	encode_sub_plain(bytes, n, k)
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
		
		if upper_len >= bytes_len { upper_len } else { upper_len << 1 }
	};
	assert!(is_power_of_2(upper_len));
	assert!(upper_len >= bytes_len);

	// tuple are only used here
	use itertools::Itertools;

	// pad the incoming bytes with trailing 0s
	// so we get a buffer of size `N` in `GF` symbols
	let zero_bytes_to_add = n * 2 - bytes_len;
	let mut elm_data = Vec::with_capacity(n);
	let zeros = std::iter::repeat(&0u8).take(zero_bytes_to_add);
	for (first, second) in bytes.iter().chain(zeros).tuples() {
		elm_data.push(Additive(Elt::from_be_bytes([*first, *second])));
	}

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
#[cfg(target_feature = "avx")]
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
		
		if upper_len >= bytes_len { upper_len } else { upper_len << 1 }
	};
	assert!(is_power_of_2(upper_len));
	assert!(upper_len >= bytes_len);

    // tuples are only used here
    use itertools::Itertools;

	// pad the incoming bytes with trailing 0s
	// so we get a buffer of size `N` in `GF` symbols
	let zero_bytes_to_add = n * 2 - bytes_len;
	let mut elm_data = Vec::with_capacity(n);
	let zeros = std::iter::repeat(&0u8).take(zero_bytes_to_add);
	for (first, second) in bytes.iter().chain(zeros).tuples() {
		elm_data.push(Additive(Elt::from_be_bytes([*first, *second])));
	}

	// update new data bytes with zero padded bytes
	// `l` is now `GF(2^16)` symbols
	let elm_len = elm_data.len();
	assert_eq!(elm_len, n);

	let mut codeword = elm_data.clone();
	encode_low_faster8(&elm_data[..], k, &mut codeword[..], n);

	Ok(codeword)
}

#[cfg(target_feature = "avx")]
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
		encode_low_faster8(&data2[..], k, &mut parity2[..], n);

		assert_eq!(parity1, parity2);
	}
	
	
	#[cfg(target_feature = "avx")]
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
