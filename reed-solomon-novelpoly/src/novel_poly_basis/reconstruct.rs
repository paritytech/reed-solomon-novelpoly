use super::*;

pub fn reconstruct_sub(
	codewords: &[Option<Additive>],
	erasures: &[bool],
	n: usize,
	k: usize,
	error_poly: &[Multiplier; FIELD_SIZE],
) -> Result<Vec<u8>> {
	assert!(is_power_of_2(n), "Algorithm only works for 2^i sizes for N");
	assert!(is_power_of_2(k), "Algorithm only works for 2^i sizes for K");
	assert_eq!(codewords.len(), n);
	assert!(k <= n / 2);

	// the first k suffice for the original k message codewords
	let recover_up_to = k; // n;

	// The recovered _payload_ chunks AND parity chunks
	let mut recovered = vec![Additive(0); recover_up_to];

	// get rid of all `None`s
	let mut codeword = codewords
		.into_iter()
		.enumerate()
		.map(|(idx, sym)| {
			// fill the gaps with `0_u16` codewords
			if let Some(sym) = sym {
				(idx, *sym)
			} else {
				(idx, Additive(0))
			}
		})
		.map(|(idx, codeword)| {
			if idx < recovered.len() {
				recovered[idx] = codeword;
			}
			codeword
		})
		.collect::<Vec<Additive>>();

	// filled up the remaining spots with 0s
	assert_eq!(codeword.len(), n);

	//---------Erasure decoding----------------

	algorithms::decode_main(&mut codeword[..], recover_up_to, &erasures[..], &error_poly[..], n);

	for idx in 0..recover_up_to {
		if erasures[idx] {
			recovered[idx] = codeword[idx];
		};
	}

	let mut recovered_bytes = Vec::with_capacity(recover_up_to * 2);
	recovered.into_iter().take(k).for_each(|x| recovered_bytes.extend_from_slice(&x.0.to_be_bytes()[..]));
	Ok(recovered_bytes)
}


/// each shard contains one symbol of one run of erasure coding
pub fn reconstruct<'a, S: Shard>(received_shards: Vec<Option<S>>, validator_count: usize) -> Result<Vec<u8>> {
	let params = CodeParams::derive_parameters(validator_count, recoverablity_subset_size(validator_count))?;

	let rs = params.make_encoder();
	rs.reconstruct(received_shards)
}
