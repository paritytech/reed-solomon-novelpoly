use super::*;

use crate::field::f2e16;
use crate::WrappedShard;
use assert_matches::assert_matches;
use rand::distributions::Uniform;
use rand::prelude::*;
use rand::rngs::SmallRng;
use rand::seq::index::IndexVec;
use rand::thread_rng;
use reed_solomon_tester::*;
use quickcheck::{Arbitrary, Gen, QuickCheck};

/// Generate a random index
fn rand_gf_element() -> Additive {
	let mut rng = thread_rng();
	let uni = Uniform::<Elt>::new_inclusive(0, ONEMASK);
	Additive(uni.sample(&mut rng))
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
		assert_matches! {
			CodeParams::derive_parameters(validator_count, recoverablity_subset_size(validator_count)),
			Ok(CodeParams { n, k, wanted_n }) => {
				assert_eq!(wanted_n, validator_count);
				assert!(validator_count <= n, "vc={} <= n={} violated", validator_count, n);
				assert!(validator_count / 3 >= k - 1, "vc={} / 3 >= k={} violated", validator_count, k);
				assert!(validator_count >= (k-1) *3, "vc={} <= k={} *3  violated", validator_count, k);
			}
		}
	}
}

#[test]
fn flt_back_and_forth() {
	const N: usize = 128;

	let mut data = (0..N).map(|_x| rand_gf_element()).collect::<Vec<Additive>>();
	let expected = data.clone();

	afft(&mut data, N, N / 4);

	// make sure something is done
	assert!(data.iter().zip(expected.iter()).filter(|(a, b)| { a != b }).count() > 0);

	inverse_afft(&mut data, N, N / 4);

	itertools::assert_equal(data, expected);
}

#[test]
fn sub_encode_decode() -> Result<()> {
	let mut rng = rand::rngs::SmallRng::from_seed(SMALL_RNG_SEED);

	const N: usize = 32;
	const K: usize = 4;

	const K2: usize = K * 2;
	let mut data = [0u8; K2];
	rng.fill_bytes(&mut data[..]);

	let codewords = encode_sub_plain(&data, N, K)?;
	let mut codewords = codewords.into_iter().map(Some).collect::<Vec<_>>();
	assert_eq!(codewords.len(), N);
	codewords[0] = None;
	codewords[1] = None;
	codewords[2] = None;
	codewords[N - 3] = None;
	codewords[N - 2] = None;
	codewords[N - 1] = None;

	let erasures = codewords.iter().map(|x| x.is_none()).collect::<Vec<bool>>();

	// Evaluate error locator polynomial only once
	let mut error_poly_in_log = [Multiplier(0); FIELD_SIZE];
	eval_error_polynomial(&erasures[..], &mut error_poly_in_log[..], FIELD_SIZE);

	let reconstructed = reconstruct_sub(&codewords[..], &erasures[..], N, K, &error_poly_in_log)?;
	itertools::assert_equal(data.iter(), reconstructed.iter().take(K2));
	Ok(())
}

#[cfg(all(target_feature = "avx", feature = "avx"))]
#[test]
fn sub_encode_faster8_decode_plain() -> Result<()> {
	let mut rng = rand::rngs::SmallRng::from_seed(SMALL_RNG_SEED);

	const N: usize = 64;
	const K: usize = 16;

	const K2: usize = K * 2;
	let mut data = [0u8; K2];
	rng.fill_bytes(&mut data[..K2]);

	let codewords = encode_sub_faster8(&data, N, K)?;
	let codewords_expected = encode_sub_plain(&data, N, K)?;
	itertools::assert_equal(codewords.iter(), codewords_expected.iter());

	let mut codewords = codewords.into_iter().map(Some).collect::<Vec<_>>();
	assert_eq!(codewords.len(), N);
	codewords[0] = None;
	codewords[1] = None;
	codewords[2] = None;
	codewords[N - 3] = None;
	codewords[N - 2] = None;
	codewords[N - 1] = None;

	let erasures = codewords.iter().map(|x| x.is_none()).collect::<Vec<bool>>();

	// Evaluate error locator polynomial only once
	let mut error_poly_in_log = [Multiplier(0); FIELD_SIZE];
	eval_error_polynomial(&erasures[..], &mut error_poly_in_log[..], FIELD_SIZE);

	let reconstructed = reconstruct_sub(&codewords[..], &erasures[..], N, K, &error_poly_in_log)?;
	itertools::assert_equal(data.iter(), reconstructed.iter().take(K2));
	Ok(())
}

// for shards of length 1
fn wrapped_shard_len1_as_gf_sym(w: &WrappedShard) -> Additive {
	let val = AsRef::<[[u8; 2]]>::as_ref(w)[0];
	Additive(u16::from_be_bytes(val))
}

#[test]
fn sub_eq_big_for_small_messages() {
	const N_WANTED_SHARDS: usize = 128;
	const N: usize = N_WANTED_SHARDS;
	const K: usize = 32;

	const K2: usize = K * 2;

	// assure the derived sizes match
	let rs = CodeParams::derive_parameters(N_WANTED_SHARDS, N_WANTED_SHARDS / 3).unwrap();
	assert_eq!(rs.n, N);
	assert_eq!(rs.k, K);

	// create random predictable bytes
	// and create a message that results in 1 GF element symbols
	// per validator
	let data = {
		let mut rng = SmallRng::from_seed(SMALL_RNG_SEED);
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
	let mut error_poly_in_log = [Multiplier(0); FIELD_SIZE];
	eval_error_polynomial(&erasures[..], &mut error_poly_in_log[..], FIELD_SIZE);

	let reconstructed_sub = reconstruct_sub(&codewords_sub[..], &erasures[..], N, K, &error_poly_in_log).unwrap();
	let reconstructed = reconstruct(codewords, rs.n).unwrap();
	itertools::assert_equal(reconstructed.iter().take(K2), reconstructed_sub.iter().take(K2));
	itertools::assert_equal(reconstructed.iter().take(K2), data.iter());
	itertools::assert_equal(reconstructed_sub.iter().take(K2), data.iter());
}

#[test]
fn roundtrip_for_large_messages() -> Result<()> {
	const N: usize = 1024;
	const K: usize = 256;
	const N_WANTED_SHARDS: usize = K * 3 + 2; // construct a number that works for deriving parmeters

	const K2: usize = K * 2;

	// assure the derived sizes match
	let rs = CodeParams::derive_parameters(N_WANTED_SHARDS, recoverablity_subset_size(N_WANTED_SHARDS))
		.expect("Const test parameters are ok. qed");
	assert_eq!(rs.n, N);
	assert_eq!(rs.k, K);

	// make sure each shard is more than one byte to
	// test the shard size
	// in GF symbols
	let shard_length: usize = 10;

	let required_payload_size = K2 * shard_length;
	let payload = &BYTES[K2..][..required_payload_size];
	// let payload = &BYTES[..];

	let mut shards = encode::<WrappedShard>(payload, N_WANTED_SHARDS).expect("Const test parameters are ok. qed");

	for (idx, shard) in shards.iter().enumerate() {
		let raw_shard = AsRef::<[[u8; 2]]>::as_ref(&shard);
		assert_eq!(
			shard_length,
			raw_shard.len(),
			"Shard #{} has an unxpected length {} (expected: {})",
			idx,
			raw_shard.len(),
			shard_length
		);
	}

	let (received_shards, dropped_indices) = deterministic_drop_shards_clone(&mut shards, rs.n, rs.k);

	let reconstructed_payload = reconstruct::<WrappedShard>(received_shards, N_WANTED_SHARDS).unwrap();

	assert_recovery(payload, &reconstructed_payload, dropped_indices, rs.n, rs.k);

	// verify integrity with criterion tests
	roundtrip_w_drop_closure::<_, _, _, SmallRng, WrappedShard, _>(
		encode,
		reconstruct,
		payload,
		N_WANTED_SHARDS,
		deterministic_drop_shards,
	)?;

	roundtrip_w_drop_closure::<_, _, _, SmallRng, WrappedShard, _>(
		encode,
		reconstruct,
		payload,
		N_WANTED_SHARDS,
		drop_random_max,
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
            let res = roundtrip_w_drop_closure::<'_,_,_,_,SmallRng, WrappedShard, _>(
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

// Roughly one Elt per validator payload
simplicissimus!(case_1: validators: 10, payload: 16);

// Unit payload, but mayn validators
simplicissimus!(case_2: validators: 100, payload: 1);

// Common case of way ore payload than validators
simplicissimus!(case_3: validators: 3, payload: 10);
simplicissimus!(case_4: validators: 4, payload: 10);
simplicissimus!(case_5: validators: 4, payload: 2 /* bytes = 1 field element */);
simplicissimus!(case_6: validators: 4, payload: 100);
simplicissimus!(case_7: validators: 4, payload: 100);

// Way more validators than payload bytes
simplicissimus!(case_8: validators: 2003, payload: 17);

#[test]
fn flt_roundtrip_small() {
	const N: usize = 16;
	const EXPECTED: [Additive; N] =
		unsafe { std::mem::transmute([1_u16, 2, 3, 5, 8, 13, 21, 44, 65, 0, 0xFFFF, 2, 3, 5, 7, 11]) };

	let mut data = EXPECTED;

	f2e16::afft(&mut data, N, N / 4);

	println!("novel basis(rust):");
	data.iter().for_each(|sym| {
		print!(" {:04X}", sym.0);
	});
	println!();

	f2e16::inverse_afft(&mut data, N, N / 4);
	itertools::assert_equal(data.iter(), EXPECTED.iter());
}

#[test]
fn ported_c_test() {
	const N: usize = 256;
	const K: usize = 8;

	//-----------Generating message----------
	//message array
	let mut data = [Additive(0); N];

	for i in 0..K {
		//filled with random numbers
		data[i] = Additive((i * i % ONEMASK as usize) as u16);
		// data[i] = rand_gf_element();
	}

	assert_eq!(data.len(), N);

	println!("Message(Last n-k are zeros): ");
	for i in 0..K {
		print!("{:04x} ", data[i].0);
	}
	println!();

	//---------encoding----------
	let mut codeword = [Additive(0); N];

	if K + K > N && false {
		let (data_till_t, data_skip_t) = data.split_at_mut(N - K);
		f2e16::encode_high(data_skip_t, K, data_till_t, &mut codeword[..], N);
	} else {
		f2e16::encode_low(&data[..], K, &mut codeword[..], N);
	}

	// println!("Codeword:");
	// for i in K..(K+100) {
	// print!("{:04x} ", codeword[i]);
	// }
	// println!("");

	//--------erasure simulation---------

	// Array indicating erasures
	let mut erasure = [false; N];

	let erasures_iv = if false {
		// erase random `(N-K)` codewords
		let mut rng = rand::thread_rng();
		let erasures_iv: IndexVec = rand::seq::index::sample(&mut rng, N, N - K);

		erasures_iv
	} else {
		IndexVec::from((0..(N - K)).collect::<Vec<usize>>())
	};
	assert_eq!(erasures_iv.len(), N - K);

	for i in erasures_iv {
		//erasure codeword symbols
		erasure[i] = true;
		codeword[i] = Additive(0);
	}

	//---------Erasure decoding----------------
	let mut log_walsh2: [Multiplier; FIELD_SIZE] = [Multiplier(0); FIELD_SIZE];

	f2e16::eval_error_polynomial(&erasure[..], &mut log_walsh2[..], FIELD_SIZE);

	// TODO: Make print_sha256 polymorphic
	// print_sha256("log_walsh2", &log_walsh2);

	f2e16::decode_main(&mut codeword[..], K, &erasure[..], &log_walsh2[..], N);

	println!("Decoded result:");
	for i in 0..N {
		// the data word plus a few more
		print!("{:04x} ", codeword[i].0);
	}
	println!();

	for i in 0..K {
		//Check the correctness of the result
		if data[i] != codeword[i] {
			eprintln!("🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍🐍");
			panic!("Decoding ERROR! value at [{}] should={:04x} vs is={:04x}", i, data[i].0, codeword[i].0);
		}
	}
	println!(
		r#">>>>>>>>> 🎉🎉🎉🎉
>>>>>>>>> > Decoding is **SUCCESS**ful! 🎈
>>>>>>>>>"#
	);
}

#[test]
fn test_code_params() {
	assert_matches!(CodeParams::derive_parameters(0, recoverablity_subset_size(0)), Err(_));

	assert_matches!(CodeParams::derive_parameters(1, recoverablity_subset_size(1)), Err(_));

	assert_eq!(
		CodeParams::derive_parameters(2, recoverablity_subset_size(2)),
		Ok(CodeParams { n: 2, k: 1, wanted_n: 2 })
	);

	assert_eq!(
		CodeParams::derive_parameters(3, recoverablity_subset_size(3)),
		Ok(CodeParams { n: 4, k: 1, wanted_n: 3 })
	);

	assert_eq!(
		CodeParams::derive_parameters(4, recoverablity_subset_size(4)),
		Ok(CodeParams { n: 4, k: 2, wanted_n: 4 })
	);

	assert_eq!(
		CodeParams::derive_parameters(100, recoverablity_subset_size(100)),
		Ok(CodeParams { n: 128, k: 32, wanted_n: 100 })
	);
}

#[test]
fn shard_len_is_reasonable() {
	let rs = CodeParams { n: 16, k: 4, wanted_n: 5 }.make_encoder();

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

#[derive(Clone, Debug)]
struct ArbitraryData(Vec<u8>);

impl Arbitrary for ArbitraryData {
	fn arbitrary(g: &mut Gen) -> Self {
		// Limit the len to 1 mib, otherwise the test will take forever
		let len = (u32::arbitrary(g) % (1024 * 1024)).max(2);

		let data: Vec<u8> = (0..len).map(|_| u8::arbitrary(g)).collect();

		ArbitraryData(data)
	}
}

#[test]
fn round_trip_systematic_quickcheck() {
	fn property(available_data: ArbitraryData, n_validators: u16) {
		let n_validators = n_validators.max(2);
		let rs = CodeParams::derive_parameters(n_validators as usize, (n_validators as usize - 1) / 3 + 1)
			.unwrap()
			.make_encoder();
		let kpow2 = rs.k;
		let chunks = rs.encode::<WrappedShard>(&available_data.0).unwrap();
		let mut res = rs.reconstruct_from_systematic(chunks.into_iter().take(kpow2).collect()).unwrap();
		res.truncate(available_data.0.len());
		assert_eq!(res, available_data.0);
	}

	QuickCheck::new().quickcheck(property as fn(ArbitraryData, u16))
}

#[test]
fn round_trip_quickcheck() {
	fn property(available_data: ArbitraryData, n_validators: u16) {
		let n_validators = n_validators.max(2);
		let wanted_k = (n_validators as usize - 1) / 3 + 1;
		let rs = CodeParams::derive_parameters(n_validators as usize, wanted_k).unwrap().make_encoder();
		let chunks = rs.encode::<WrappedShard>(&available_data.0).unwrap();
		let mut res = rs.reconstruct(chunks.into_iter().take(wanted_k).map(|s| Some(s)).collect()).unwrap();
		res.truncate(available_data.0.len());
		assert_eq!(res, available_data.0);
	}

	QuickCheck::new().quickcheck(property as fn(ArbitraryData, u16))
}
