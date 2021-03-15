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

/// Does not care about power of 2 requirements.
///
/// Covers the 1/3 case.
pub const fn recoverablity_subset_size(n_wanted_shards: usize) -> usize {
	n_wanted_shards / 3
}


#[test]
fn one_third() {
	assert_eq!(recoverablity_subset_size(0), 0);
	assert_eq!(recoverablity_subset_size(1), 0);
	assert_eq!(recoverablity_subset_size(2), 0);
	assert_eq!(recoverablity_subset_size(3), 1);
	assert_eq!(recoverablity_subset_size(4), 1);
	assert_eq!(recoverablity_subset_size(5), 1);
	assert_eq!(recoverablity_subset_size(6), 2);
	assert_eq!(recoverablity_subset_size(8), 2);
	assert_eq!(recoverablity_subset_size(11), 3);

	assert_eq!(recoverablity_subset_size(173), 57);
	assert_eq!(recoverablity_subset_size(174), 58);
}
