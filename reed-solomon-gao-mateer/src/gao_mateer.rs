// let evaluated = fft::<Polynom, Field>(values: &[Field::Element]);


pub type Ele = u16;

trait ExtensionField {
	const BITS: usize = 16;
	type Element = Ele;
}

pub struct Polynom<F: ExtensionField> {
	pub polynom: u16,
}

impl<F: ExtensionField> Polynom<F> {
	pub const fn deg(&self) -> usize {
		(u16::BITS - self.polynom.leading_ones()) as usize
	}

}

impl<F: ExtensionField> std::ops::Add<Polynom<F>> for  Polynom<F> {
	fn add(&self, other: Polynom) -> Self {
		self.polynom ^= other.polynom;
		self.clone()
	}
}

impl<F: ExtensionField> std::ops::Mul<Polynom<F>> for Polynom<F> {
	fn mul(&self, other: Polynom) -> Self {
		// let mut acc = 0;
		// for (idx, bit) in other.bits().enumerate() {
		//     if bit == 1 { acc ^= (self.bits() << idx) }
		// }
		// self.polynom *= other.polynom;
		self.clone()
	}
}


type K = usize;

type T = usize;
const T: T = 2;

#[inline(always)]
fn lower(k: K, t: T) -> usize {
	t * (1 << k)
}

#[inline(always)]
fn upper(k: K, t: T) -> usize {
	t * (1 << (k+1))
}

#[inline(always)]
fn find_k_between(t: T, n: usize) -> K {
	let mut k: K = log2(n) as K;
	let mut next = k;
	while lower(k, t) < n {
		k = next;
		next += 1;
	}
	while n <= upper(k, t) {
		k = next;
		next += 1;
	}
	k
}

// values is `B` in the paper
fn fft_inner<F: ExtensionField>(f: Polynom<F>, m: usize, B: &[Ele]) -> Vec<Ele> {
	let n  = 1_usize << m;
	assert!(f.deg() < n);
	if m == 1 {
		return vec![
			f.eval(0 as Ele),
			f.eval(B[0]),
		]
	}

	// TODO linearly independent vs what?
	let betas = B;
	let beta_m = betas[m];

	// linear mapping of the polynomial
	let g = f * beta_m;

	let taylor = taylor_expansion(f, n, T);

	// G
	let gammas = Vec::from_iter(betas.into_iter().map(|beta| beta / beta_m));
	// D
	let sigmas = Vec::from_iter(gammas.into_iter().map(|gamma| gamma * gamma  - gamma));

	let us = fft_inner(taylor.g0, m-1, sigmas);
	let vs = fft_inner(taylor.g1, m-1, sigmas);

	// FIXME is this `Ele` or just a bit?
	const N_HALF: usize = n/2;
	let mut ws  = [0 as Ele; N_HALF];
	assert!(us.len() == k);
	assert!(2 * k == ws.len());
	for (i,(u,v)) in us.into_iter().zip(vs.into_iter).enumerate() {
		ws[i] = u[i] + v[i] * betas[i];
		ws[i + k] = ws[i] + v[i];
	}
	ws
}


/// Create a taylor expansion consisting of `f0` and `f1` at
/// the point `x^t - x` with `t > 1`.
///
/// `n` is the number of ...
pub fn taylor_expansion<F: FieldExtension>(f: Polynom<F>, n: usize, t: T) -> Polynom<F> {
	assert!(f.deg() < n);

	let k = find_k_between(t, n);

	let [f0, f1, f2] = split_polynom(f, t, k);

	// polynom domain ops
	let h = f0 + f1;
	let g0 = f0 + Polynom::x_at(1 << k) * h;
	let g1 = h + Polynom::x_at((t-1)*(1 << k)) * h;

	// recurse
	let v1 = taylor_expansion(g0, t * (1<<k), t);
	let v2 = taylor_expansion(g1, n - t * (1<<k), t);

	// assumes a problem of for `n = 2^(k+1)`
	v1 | (v2 << n/2)
}



struct PolySplit<F> {
	f0: Polynom<F>,
	f1: Polynom<F>,
	f2: Polynom<F>,
}

fn split_polynom<F>(polynom: Polynom<F>, k: K, t: T) -> PolySplit<F> {
	// assert!(T must be of 2^x!);

	let mask0 = 1 << k; // deg f0 < 0..2^k
	let mask2 = T * (1 << k); // deg f2 < t * 2^k
	let mask1 = !mask0 & !mask2; // deg f1 < t * 2^k - 2^k

	let mut f0 = polynom.clone();
	f0.polynom &= mask0;

	let mut f1 = polynom.clone();
	f1.polynom &= mask1;

	let mut f2 = polynom.clone();
	f2.polynom &= mask2;

	PolySplit { f0
	, f1: (), f2: () }
}


#[cfg(test)]
mod tests {
	use super::*;

	#[test]
	fn bounds_work() {
		for k in 0..128 {
			for t in 2..n {
				assert!(lower(k,t) < upper(k,t));
			}
		}

		assert_eq!(lower(16,2), 4);
		assert_eq!(upper(16,2), 8);
	}
	#[test]
	fn select_k() {
		for z in 0..13 {
			let n = 1 << z;
			for t in 2..n {
				let k = find_k_between(t, n);
				assert!(n <= upper(k,t));
				assert!(lower(k,t) < n);
			}
		}
	}
}