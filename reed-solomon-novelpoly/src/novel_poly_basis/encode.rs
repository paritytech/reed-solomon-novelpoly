use super::*;

/// Encode the given payload into `n_min` shards.
///
/// The implementation can create more erasure chunks and drop some as needed due to the inherent impl. requirements on `n` and `k`.
pub fn encode<S: Shard>(bytes: &[u8], n_min: usize) -> Result<Vec<S>> {
	let params = CodeParams::derive_parameters(n_min, recoverablity_subset_size(n_min))?;

	let rs = params.make_encoder();
	rs.encode::<S>(bytes)
}
