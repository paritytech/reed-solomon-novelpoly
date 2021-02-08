use std::cmp;
use std::fmt;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Clone, Copy, Eq, PartialOrd, Ord)]
#[repr(align(2))]
struct Element(u16);

impl fmt::Display for Element {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		write!(f, "{}", self.0)
	}
}

impl fmt::Debug for Element {
	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
		write!(f, "{}", self.0)
	}
}

// impl fmt::Debug for Vec<Element> {
//     fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
//         f.write_str("[")?;
//         let mut iter = self.iter();
//         if let Some(first) = iter.next() {
//             write!(f, "{}", self.0)?;
//             for element in iter {
//                 write!(f, ", {}", self.0)?;
//             }
//         }
//         f.write_str("]")?;
//         Ok(())
//     }
// }

impl Element {
	#[inline(always)]
	const fn zero() -> Self {
		Self(0u16)
	}

	#[inline(always)]
	const fn one() -> Self {
		Self(1u16)
	}

	#[inline(always)]
	fn log2(&self) -> Element {
		Self(log2(self.0))
	}

	#[inline(always)]
	fn is_power_of_2(&self) -> bool {
		is_power_of_2(self.0)
	}

	fn pow_mod(&self, mut exp: Element, modulo: u64) -> Self {
		let mut val = self.0 as u64;
		let mut res = 1_u64;
		while exp != 0 {
			if exp.0 & 0x1 != 0 {
				res *= val;
				res %= modulo;
			}
			val *= val;
			val %= modulo;

			exp >>= 1;
		}
		Self(res as u16)
	}
}

// #[target_feature = "step_trait"]
// impl Step for Element {
//     fn steps_between(start: &Self, end: &Self) -> Option<usize> {
//         let val = end.0.saturating_sub(start.0);
//         Some(val as usize + 1_usize)
//     }
// }

impl<R> std::ops::BitXor<R> for Element
where
	R: Into<Element>,
{
	type Output = Self;
	fn bitxor(self, rhs: R) -> Self::Output {
		Self(self.0 ^ rhs.into().0)
	}
}

impl<R> std::ops::BitXorAssign<R> for Element
where
	R: Into<Element>,
{
	fn bitxor_assign(&mut self, rhs: R) {
		self.0 ^= rhs.into().0
	}
}

impl<R> std::ops::Mul<R> for Element
where
	R: Into<Element>,
{
	type Output = Self;
	fn mul(self, rhs: R) -> Self::Output {
		let a = self.0;
		let b = rhs.into().0;
		Self(raw_mul(a, b))
	}
}

impl<R> std::ops::Rem<R> for Element
where
	R: Into<Element>,
{
	type Output = Self;
	fn rem(self, rhs: R) -> Self::Output {
		let a = self.0;
		let b = rhs.into().0;
		Self(raw_mod(a, b))
	}
}

impl<R> std::ops::Add<R> for Element
where
	R: Into<Element>,
{
	type Output = Self;
	fn add(self, rhs: R) -> Self::Output {
		let rhs = rhs.into();
		Self(self.0 + rhs.0)
	}
}

impl<R> std::ops::Sub<R> for Element
where
	R: Into<Element>,
{
	type Output = Self;
	fn sub(self, rhs: R) -> Self::Output {
		let rhs = rhs.into();
		Self(self.0 - rhs.0)
	}
}

impl<R> std::ops::Shl<R> for Element
where
	R: Into<Element>,
{
	type Output = Self;
	fn shl(self, rhs: R) -> Self::Output {
		let rhs = rhs.into();
		Self(self.0 << rhs.0)
	}
}
impl<R> std::ops::Shr<R> for Element
where
	R: Into<Element>,
{
	type Output = Self;
	fn shr(self, rhs: R) -> Self::Output {
		let rhs = rhs.into();
		Self(self.0 >> rhs.0)
	}
}
impl<R> std::ops::ShrAssign<R> for Element
where
	R: Into<Element>,
{
	fn shr_assign(&mut self, rhs: R) {
		let rhs = rhs.into();
		self.0 >>= rhs.0;
	}
}
impl<R> std::ops::ShlAssign<R> for Element
where
	R: Into<Element>,
{
	fn shl_assign(&mut self, rhs: R) {
		let rhs = rhs.into();
		self.0 <<= rhs.0;
	}
}

impl From<&u16> for Element {
	#[inline(always)]
	fn from(inner: &u16) -> Self {
		Self(*inner)
	}
}

impl From<u16> for Element {
	#[inline(always)]
	fn from(inner: u16) -> Self {
		Self(inner)
	}
}

impl From<&Element> for Element {
	#[inline(always)]
	fn from(inner: &Self) -> Self {
		*inner
	}
}

impl From<&usize> for Element {
	#[inline(always)]
	fn from(inner: &usize) -> Self {
		Self(*inner as u16)
	}
}

impl From<usize> for Element {
	#[inline(always)]
	fn from(inner: usize) -> Self {
		Self(inner as u16)
	}
}

impl From<i32> for Element {
	#[inline(always)]
	fn from(inner: i32) -> Self {
		Self(inner as u16)
	}
}

impl PartialEq<usize> for Element {
	fn eq(&self, other: &usize) -> bool {
		self.0.eq(&(*other as u16))
	}
}
impl PartialEq<i32> for Element {
	fn eq(&self, other: &i32) -> bool {
		self.0.eq(&(*other as u16))
	}
}
impl PartialEq<u16> for Element {
	fn eq(&self, other: &u16) -> bool {
		self.0.eq(other)
	}
}
impl PartialEq<Element> for Element {
	fn eq(&self, other: &Element) -> bool {
		self.0.eq(&other.0)
	}
}

const fn log2(mut x: u16) -> u16 {
	let mut o = 0;
	while x > 1 {
		x >>= 1;
		o += 1;
	}
	o
}

const fn is_power_of_2(x: u16) -> bool {
	return x > 0_u16 && x & (x - 1) == 0;
}

#[inline(always)]
fn raw_mul(a: u16, b: u16) -> u16 {
	if a.saturating_mul(b) == 0 {
		return 0;
	}
	let mut o = 0;
	for i in 0..(log2(b) + 1) {
		if (b & (1 << i)) != 0x0 {
			o ^= a << i
		}
	}
	o
}

#[inline(always)]
const fn raw_mod(mut a: u16, b: u16) -> u16 {
	let mut alog = log2(a);
	let blog = log2(b);
	while alog >= blog {
		if a & (1 << alog) != 0x0 {
			a ^= b << (alog - blog);
		}
		alog -= 1;
	}
	a
}

#[derive(Debug, thiserror::Error)]
enum Error {
	#[error("The provided modulues {0:?} is bad")]
	Badpd(Element),
}

//[derive(Debug)]
struct BinaryField {
	pd: Element,
	height: Element,
	order: Element,

	cache: Vec<Element>,
	invcache: Vec<Option<usize>>,
}

impl BinaryField {
	fn setup(&mut self) -> Result<()> {
		// XXX why 80?
		let pd = self.pd.0;
		let order = self.order.0 as usize;

		for base in 2..cmp::min(pd - 1, 80_u16) {
			let mut powers: Vec<Element> = vec![Element::one()];
			'p: while powers.len() < order + 2 {
				let previous = powers.last().unwrap();
				let val = raw_mod(raw_mul(previous.0, base), pd);
				powers.push(val.into());
				if val == 1 {
					break 'p;
				}
			}
			let _ = powers.pop();
			if powers.len() == order {
				self.cache = powers.clone();
				self.cache.extend(powers.iter().cloned());
				self.invcache = vec![None; order + 1];
				for (idx, p) in powers.into_iter().enumerate() {
					self.invcache[p.0 as usize] = Some(idx);
				}
				return Ok(());
			}
		}
		Err(Error::Badpd(self.pd))
	}

	pub fn new(pd: Element) -> Result<Self> {
		let height = pd.log2();
		let mut field = Self {
			pd,
			height,
			order: (Element::one() << height) - Element::one(),
			cache: Default::default(),
			invcache: Default::default(),
		};

		field.setup()?;

		Ok(field)
	}

	// binary field special
	fn add(&self, x: Element, y: Element) -> Element {
		x ^ y
	}

	fn sub(&self, x: Element, y: Element) -> Element {
		self.add(x, y)
	}

	fn mul(&self, x: Element, y: Element) -> Element {
		if x.0 as u32 * y.0 as u32 == 0 {
			Element::zero()
		} else {
			let idx = self.invcache[x.0 as usize].unwrap() + self.invcache[y.0 as usize].unwrap();
			self.cache[idx]
		}
	}

	fn sqr(&self, x: Element) -> Element {
		if x == Element::zero() {
			Element::zero()
		} else {
			let idx = (self.invcache[x.0 as usize].unwrap() << 1_usize) % self.order.0 as usize;
			self.cache[idx]
		}
	}

	fn div(&self, x: Element, y: Element) -> Element {
		if x == 0 {
			Element::zero()
		} else {
			let idx: Element = self.order + self.invcache[x.0 as usize].unwrap() - self.invcache[y.0 as usize].unwrap();
			self.cache[idx.0 as usize]
		}
	}

	fn inv(&self, x: Element) -> Element {
		assert_ne!(x, Element::zero());
		let idx = ((self.order - self.invcache[x.0 as usize].unwrap()) % self.order).0;
		self.cache[idx as usize]
	}

	fn exp(&self, x: Element, p: Element) -> Element {
		if p == Element::zero() {
			Element::one()
		} else if x == Element::zero() {
			Element::one()
		} else {
			let idx = (self.invcache[x.0 as usize].unwrap() * (p.0 as usize)) % self.order.0 as usize;
			self.cache[idx]
		}
	}

	fn multi_inv(&self, values: Vec<Option<Element>>) -> Vec<Element> {
		let mut partials = vec![Element::one()];
		for i in 0..values.len() {
			let last = partials.last().unwrap().clone();
			partials.push(self.mul(last, values[i].unwrap_or(Element::zero())))
		}
		let mut inv = self.inv(partials.last().unwrap().clone());
		let mut outputs = vec![Element::zero(); values.len()];
		for (i, value) in values.into_iter().enumerate().rev() {
			outputs[i] = if value.is_some() { self.mul(partials[i], inv) } else { Element::zero() };
			inv = self.mul(inv, value.unwrap_or(Element::one()))
		}
		outputs
	}

	// TODO XXX understand why there are two diff impls
	fn div_(&self, x: Element, y: Element) -> Element {
		self.mul(x, self.inv(y))
	}

	// Evaluate a polynomial at a point
	fn eval_poly_at(&self, p: &[Element], x: Element) -> Element {
		let mut y = Element::zero();
		let mut power_of_x = Element::one();
		for (i, &p_coeff) in p.into_iter().enumerate() {
			y ^= self.mul(power_of_x, p_coeff);
			power_of_x = self.mul(power_of_x, x);
		}
		y
	}

	// Arithmetic for polynomials
	fn add_polys(&self, a: Vec<Element>, b: Vec<Element>) -> Vec<Element> {
		let deg = cmp::max(a.len(), b.len());
		let mut res = a.clone();
		if deg < b.len() {
			res.extend(std::iter::repeat(Element::zero()).take(b.len() - deg))
		}
		for i in 0..deg {
			res[i] ^= b.get(i).cloned().unwrap_or(Element::zero());
		}
		res
	}

	fn sub_polys(&self, a: Vec<Element>, b: Vec<Element>) -> Vec<Element> {
		self.add_polys(a, b)
	}

	fn mul_by_const(&self, a: &[Element], c: Element) -> Vec<Element> {
		a.into_iter().map(move |x| self.mul(*x, c)).collect::<Vec<Element>>()
	}

	fn mul_polys(&self, a: Vec<Element>, b: Vec<Element>) -> Vec<Element> {
		let mut o = vec![Element::zero(); a.len() + b.len() - 1];
		for (i, &aval) in a.iter().enumerate() {
			for (j, &bval) in b.iter().enumerate() {
				o[i + j] ^= self.mul(aval, bval)
			}
		}
		o
	}

	fn div_polys(&self, mut a: Vec<Element>, b: Vec<Element>) -> Vec<Element> {
		assert!(a.len() >= b.len());
		let mut o = vec![];
		let mut apos = a.len() - 1_usize;
		let mut bpos = b.len() - 1_usize;
		let mut diff = apos as isize - bpos as isize;
		while diff >= 0_isize {
			let quot = self.div(a[apos], b[bpos]);
			o.insert(0, quot);
			for (i, b) in (0..bpos).into_iter().enumerate().rev() {
				a[(diff + i as isize) as usize] ^= self.mul(Element::from(b), quot);
			}
			apos -= 1_usize;
			diff -= 1_isize;
		}
		o
	}

	// Build a polynomial that returns 0 at all specified xs
	fn zpoly(&self, xs: &[Element]) -> Vec<Element> {
		let mut roots = vec![Element::one()];
		for &x in xs {
			roots.insert(0, Element::zero());
			let rl = roots.len() - 1;
			for j in 0..rl {
				let val = roots[j + 1];
				roots[j] ^= self.mul(val, x);
			}
		}
		roots
	}

	// Given p+1 y values && x values with no errors, recovers the original
	// p+1 degree polynomial.
	// Lagrange interpolation works roughly in the following way.
	// 1. Suppose you have a set of points, eg. x = [1, 2, 3], y = [2, 5, 10]
	// 2. For each x, generate a polynomial which equals its corresponding
	//    y coordinate at that point && 0 at all other points provided.
	// 3. Add these polynomials together.

	fn lagrange_interp(&self, xs: Vec<Element>, ys: Vec<Element>) -> Vec<Element> {
		// // Generate master numerator polynomial, eg. (x - x1) * (x - x2) * ... * (x - xn)
		let root = self.zpoly(&xs);
		assert_eq!(root.len(), ys.len() + 1);
		// // print(root)
		// // Generate per-value numerator polynomials, eg. for x=x2,
		// // (x - x1) * (x - x3) * ... * (x - xn), by dividing the master
		// // polynomial back by each x coordinate
		let mut nums =
			xs.iter().map(|&x| self.div_polys(root.clone(), vec![x, Element::one()])).collect::<Vec<Vec<Element>>>();
		// Generate denominators by evaluating numerator polys at each x
		let denoms = xs
			.iter()
			.zip(nums.iter())
			.map(|(&x, num)| Some(self.eval_poly_at(num, x)))
			.collect::<Vec<Option<Element>>>();
		let invdenoms = self.multi_inv(denoms);
		// Generate output polynomial, which is the sum of the per-value numerator
		// polynomials rescaled to have the right y values
		let mut b = vec![Element::zero(); ys.len()];
		let xsl = xs.len();
		for i in 0..xsl {
			let yslice = self.mul(ys[i], invdenoms[i]);
			for j in 0..ys.len() {
				if nums[i][j] != Element::zero() && ys[i] != Element::zero() {
					b[j] ^= self.mul(nums[i][j], yslice);
				}
			}
		}
		b
	}
}

fn _simple_ft(field: &BinaryField, domain: &[Element], poly: &[Element]) -> Vec<Element> {
	domain.into_iter().map(|&item| field.eval_poly_at(poly, item)).collect::<Vec<Element>>()
}

// Returns `evens` && `odds` such that{
// poly(x) = evens(x**2+kx) + x * odds(x**2+kx)
// poly(x+k) = evens(x**2+kx) + (x+k) * odds(x**2+kx)
//
// Note that this satisfies two other invariants{
//
// poly(x+k) - poly(x) = k * odds(x**2+kx)
// poly(x)*(x+k) - poly(x+k)*x = k * evens(x**2+kx)

fn cast(field: &BinaryField, poly: &[Element], k: Element) -> (Vec<Element>, Vec<Element>) {
	println!("cast");
	if poly.len() <= 2 {
		return (vec![poly[0]], vec![if poly.len() == 2 { poly[1] } else { Element::zero() }]);
	}
	assert!(is_power_of_2(dbg!(poly.len()) as u16));

	let mod_power = poly.len() >> 1_usize;
	let half_mod_power = mod_power >> 1_usize;
	let k_to_half_mod_power = field.exp(k, half_mod_power.into());
	dbg!(mod_power);
	dbg!(half_mod_power);
	dbg!(k_to_half_mod_power);

	assert_eq!(mod_power, 2 * half_mod_power); // holds since poly is single val / one bit set / 2^x

	// Calculate low = poly % (x**2 - k*x)**half_mod_power
	// && high = poly // (x**2 - k*x)**half_mod_power
	// Note that (x**2 - k*x)**n = x**2n - k**n * x**n in binary fields
	let mut low_and_high = poly.to_vec();
	dbg!(poly.len());

	{
		let (low, high) = low_and_high.split_at_mut(mod_power + half_mod_power);
		for i in 0..half_mod_power {
			low[mod_power + i] ^= field.mul(high[i], k_to_half_mod_power);
		}
	}
	{
		let (low, high) = low_and_high.split_at_mut(mod_power);
		for i in 0..half_mod_power {
			low[i + half_mod_power] ^= field.mul(high[i], k_to_half_mod_power);
		}
	}
	let (low, high) = low_and_high.split_at(mod_power);
	// Recursively compute two half-size sub-problems, low && high
	let mut low_cast = dbg!(cast(field, low, k));
	let high_cast = dbg!(cast(field, high, k));
	// Combine the results
	(
		{
			low_cast.0.extend(high_cast.0.into_iter());
			low_cast.0
		},
		{
			low_cast.1.extend(high_cast.1.into_iter());
			low_cast.1
		},
	)
}

// Returns a polynomial p2 such that p2(x) = poly(x**2+kx)
fn compose(field: &BinaryField, poly: &[Element], k: Element) -> Vec<Element> {
	if poly.len() == 2 {
		return vec![poly[0], field.mul(poly[1], k), poly[1], Element::zero()];
	}
	if poly.len() == 1 {
		let mut res = poly.to_vec();
		res.push(Element::zero());
		return res;
	}
	// Largest mod_power=2**k such that mod_power >= poly.len()/2
	assert!(is_power_of_2(poly.len() as u16));
	let mod_power: usize = poly.len() >> 1_usize;
	let k_to_mod_power = field.exp(k, mod_power.into());
	// Recursively compute two half-size sub-problems, the bottom && top half
	// of the polynomial
	let low = compose(field, &poly[..mod_power], k);
	let high = compose(field, &poly[mod_power..], k);
	// Combine them together, multiplying the top one by (x**2-k*x)**n
	// Note that (x**2 - k*x)**n = x**2n - k**n * x**n in binary fields
	let mut o = vec![Element::zero(); poly.len() << 1];
	for (i, (&low, &high)) in low.iter().zip(high.iter()).enumerate() {
		o[i] ^= low;
		o[i + mod_power] ^= field.mul(high, k_to_mod_power);
		o[i + 2 * mod_power] ^= high;
	}
	o
}

// Equivalent to [field.eval_poly_at(poly, x) for x in domain]
// Special thanks to www.math.clemson.edu/~sgao/papers/GM10.pdf for insights
// though this algorithm is not exactly identical to any algorithm in the paper
fn fft(field: &BinaryField, domain: &[Element], poly: &[Element]) -> Vec<Element> {
	println!("fft");
	// Base case: constant polynomials
	// if domain.len() == 1{
	//     return [poly[0]]
	dbg!(domain.len());
	dbg!(poly.len());
	if domain.len() <= 8 {
		return _simple_ft(field, domain, poly);
	}
	// Split the domain into two cosets A && B, where for x in A, x+offset is in B
	let offset = domain[1];
	dbg!(offset);
	// Get evens, odds such that{
	// poly(x) = evens(x**2+offset*x) + x * odds(x**2+offset*x)
	// poly(x+k) = evens(x**2+offset*x) + (x+k) * odds(x**2+offset*x)
	let (evens, odds) = cast(field, poly, offset);
	dbg!(evens.len());
	dbg!(odds.len());
	// The smaller domain D = [x**2 - offset*x for x in A] = [x**2 - offset*x for x in B]
	let cast_domain = domain.iter().step_by(2).map(|&x| field.mul(x, offset ^ x)).collect::<Vec<Element>>();
	dbg!(&cast_domain);
	// Two half-size sub-problems over the smaller domain, recovering
	// evaluations of evens && odds over the smaller domain
	let even_points = fft(field, &cast_domain[..], &evens[..]);
	let odd_points = fft(field, &cast_domain[..], &odds[..]);
	// Combine the evaluations of evens && odds into evaluations of poly
	let dl = domain.len() >> 1_usize;
	let mut o = Vec::with_capacity(dl << 1_usize);
	for i in 0..dl {
		o.push(even_points[i] ^ field.mul(domain[i * 2], odd_points[i]));
		o.push(even_points[i] ^ field.mul(domain[i * 2 + 1], odd_points[i]));
	}
	o
}

// The inverse function of fft, does the steps backwards
fn invfft(field: &BinaryField, domain: &[Element], vals: &[Element]) -> Vec<Element> {
	// Base case: constant polynomials
	if domain.len() == 1 {
		return vals.to_vec();
	}
	// if domain.len() <= 4{
	//     return field.lagrange_interp(domain, vals)
	// Split the domain into two cosets A && B, where for x in A, x+offset is in B
	let offset = domain[1];
	// Compute the evaluations of the evens && odds polynomials using the invariants{
	// poly(x+k) - poly(x) = k * odds(x**2+kx)
	// poly(x)*(x+k) - poly(x+k)*x = k * evens(x**2+kx)
	let mut even_points = vec![Element::zero(); vals.len() >> 1_usize];
	let mut odd_points = vec![Element::zero(); vals.len() >> 1_usize];
	let dl = domain.len() >> 1_usize;
	for i in 0..dl {
		let (p_of_x, p_of_x_plus_k) = (vals[i * 2], vals[i * 2 + 1]);
		let x = domain[i * 2];
		even_points[i] = field.div(field.mul(p_of_x, x ^ offset) ^ field.mul(p_of_x_plus_k, x), offset);
		odd_points[i] = field.div(p_of_x ^ p_of_x_plus_k, offset);
	}
	let cast_domain = domain.into_iter().step_by(2).map(|&x| field.mul(x, offset ^ x)).collect::<Vec<Element>>();
	// Two half-size problems over the smaller domains, recovering
	// the polynomials evens && odds
	let evens = invfft(field, &cast_domain[..], &even_points[..]);
	let odds = invfft(field, &cast_domain[..], &odd_points[..]);
	// Given evens && odds where poly(x) = evens(x**2+offset*x) + x * odds(x**2+offset*x),
	// recover poly
	let mut composed_evens = compose(field, &evens[..], offset);
	composed_evens.push(Element::zero());
	let mut composed_odds = vec![Element::zero()];
	composed_odds.extend(compose(field, &odds[..], offset).iter());
	(0..vals.len()).into_iter().map(|i| composed_evens[i] ^ composed_odds[i]).collect::<Vec<_>>()
}

// shift_polys[i][j] is the 2**j degree coefficient of the polynomial that
// evaluates to [1,1...1, 0,0....0] with 2**(i-1) ones && 2**(i-1) zeroes
static SHIFT_POLYS: &[&[usize]] = &[
	&[],
	&[1],
	&[32755, 32755],
	&[52774, 60631, 8945],
	&[38902, 5560, 44524, 12194],
	&[55266, 46488, 60321, 5401, 40130],
	&[21827, 32224, 51565, 15072, 8277, 64379],
	&[59460, 15452, 60370, 24737, 20321, 35516, 39606],
	&[42623, 56997, 25925, 15351, 16625, 47045, 38250, 17462],
	&[7575, 27410, 32434, 22187, 28933, 15447, 37964, 38186, 4776],
	&[39976, 61188, 42456, 2155, 6178, 34033, 52305, 14913, 2896, 48908],
	&[6990, 12021, 36054, 16198, 17011, 14018, 58553, 13272, 25318, 5288, 21429],
	&[16440, 34925, 14360, 22561, 43883, 36645, 7613, 26531, 8597, 59502, 61283, 53412],
];

fn invfft2(field: &BinaryField, vals: &[Element]) -> Vec<Element> {
	if vals.len() == 1 {
		return vals.to_vec();
	}
	let len_half = vals.len() >> 1_usize;
	let left = invfft2(field, &vals[..len_half]);

	let tmp = invfft2(field, &vals[len_half..]);
	let right = shift(field, &tmp[..], len_half.into());

	let mut o = vec![Element::zero(); vals.len()];
	for (j, (left, right)) in left.into_iter().zip(right.into_iter()).enumerate() {
		o[j] ^= left;
		for (i, &coeff) in SHIFT_POLYS[log2(vals.len() as u16) as usize].into_iter().enumerate() {
			o[(1 << i) + j] ^= field.mul(left ^ right, coeff.into());
		}
	}
	o
}

// Multiplies two polynomials using the FFT method
fn mul(field: &BinaryField, domain: &[Element], p1: &[Element], p2: &[Element]) -> Vec<Element> {
	assert!(p1.len() <= domain.len() && p2.len() <= domain.len());
	let values1 = fft(field, domain, p1);
	let values2 = fft(field, domain, p2);
	let values3 =
		values1.into_iter().zip(values2.into_iter()).map(|(v1, v2)| field.mul(v1, v2)).collect::<Vec<Element>>();
	invfft(field, domain, &values3[..])
}

// Generates the polynomial `p(x) = (x - xs[0]) * (x - xs[1]) * ...`
fn zpoly(field: &BinaryField, xs: Vec<Element>) -> Vec<Element> {
	if xs.is_empty() {
		return vec![Element::one()];
	}
	if xs.len() == 1 {
		return vec![xs[0], Element::one()];
	}
	let domain =
		(0_u16..(2_u16 << log2(xs.iter().max().unwrap().0) + 1)).into_iter().map(Element::from).collect::<Vec<_>>();
	let offset = domain[1];
	let z_left = zpoly(field, xs.iter().step_by(2).copied().collect());
	let z_right = zpoly(field, xs.iter().skip(1).step_by(2).copied().collect());
	mul(field, &domain[..], &z_left, &z_right)
}

// Returns q(x) = p(x + k)
fn shift(field: &BinaryField, poly: &[Element], k: Element) -> Vec<Element> {
	if poly.len() == 1 {
		return poly.to_vec();
	}
	// Largest mod_power=2**k such that mod_power >= poly.len()/2
	assert!(is_power_of_2(poly.len() as u16));
	let mod_power = poly.len() >> 1_usize;
	let k_to_mod_power = field.exp(k, Element::from(mod_power));
	// Calculate low = poly % (x+k)**mod_power
	// && high = poly // (x+k)**mod_power
	// Note that (x+k)**n = x**n + k**n for power-of-two powers in binary fields
	let mut low_and_high = poly.to_vec();
	let (low, high) = low_and_high.split_at_mut(mod_power);
	for i in 0..mod_power {
		low[i].0 ^= field.mul(high[i], k_to_mod_power).0;
	}
	[shift(field, low, k), shift(field, high, k)].concat()
}

// Interpolates the polynomial where `p(xs[i]) = vals[i]`
fn interpolate(field: &BinaryField, xs: &[Element], vals: &[Element]) -> Vec<Element> {
	assert!(!xs.is_empty());
	let domain_size = Element::one() << xs.iter().max().unwrap().log2() + Element::one();
	assert!((domain_size << 1_usize) <= (Element::one() << field.height));
	let domain = (0..domain_size.0).into_iter().map(Element::from).collect::<Vec<_>>();
	let big_domain = (0..(domain_size.0 << 1_usize)).into_iter().map(Element::from).collect::<Vec<_>>();
	let z = zpoly(field, domain.iter().filter(|&x| !xs.contains(x)).copied().collect());
	// print("z = ", z)
	let z_values = fft(field, &big_domain[..], &z);
	// print("z_values = ", z_values)
	let mut p_times_z_values = vec![Element::zero(); domain.len()];
	for (&v, &d) in vals.iter().zip(xs.into_iter()) {
		let i = d.0 as usize;
		p_times_z_values[i] = field.mul(v, z_values[i]);
	}
	// print("p_times_z_values = ", p_times_z_values)
	let p_times_z = invfft(field, &domain[..], &p_times_z_values);
	// print("p_times_z = ", p_times_z)
	let shifted_p_times_z_values = fft(field, &big_domain[..], &p_times_z[..]);
	let shifted_p_times_z_values = &shifted_p_times_z_values[domain_size.0 as usize..];
	// print("shifted_p_times_z_values =", shifted_p_times_z_values)
	let shifted_p_values = shifted_p_times_z_values
		.into_iter()
		.zip(z_values[domain_size.0 as usize..].into_iter())
		.map(|(&x, &y)| field.div(x, y))
		.collect::<Vec<Element>>();
	// print("shifted_p_values =", shifted_p_values)
	let shifted_p = invfft(field, &domain, &shifted_p_values);
	shift(field, &shifted_p, domain_size)
}

#[cfg(test)]
mod tests {
	use super::*;

	macro_rules! assert_slice_eq {
		($left:expr, $right:expr) => {
			let left_slice: &[_] = $left;
			let right_slice: &[_] = $right;
			for (idx, (&left, &right)) in left_slice.iter().zip(right_slice.iter()).enumerate()
			{
				assert_eq!(
					left,
					right,
					"Failed at idx = {}, {:?} vs {:?}",
					idx,
					&left_slice[idx.saturating_sub(2)..idx.saturating_add(2)],
					&right_slice[idx.saturating_sub(2)..idx.saturating_add(2)],
				);
			}
			assert_eq!(left_slice.len(), right_slice.len());
		};
	}

	#[test]
	fn test_pow_mod() {
		assert_eq!(Element::from(1).pow_mod(10.into(), 1024), Element::one());
		assert_eq!(Element::from(2).pow_mod(10.into(), 1024), Element::zero());
		assert_eq!(Element::from(3).pow_mod(3.into(), 1024), Element::from(27));
	}

	#[test]
	fn test_mul() {
		let field = BinaryField::new(1033.into()).unwrap();
		assert_eq!(field.mul(128.into(), 128.into()), Element::from(144));
		assert_eq!(field.mul(37.into(), 11.into()), Element::from(327));
		assert_eq!(field.mul(1.into(), 1.into()), Element::one());
		assert_eq!(field.mul(1023.into(), 2.into()), Element::from(1015));


		assert_eq!(field.mul(256.into(), 256.into()), Element::from(576));

		assert_eq!(field.mul(256.into(), Element::zero()), Element::zero());
		assert_eq!(field.mul(Element::zero(), 256.into()), Element::zero());
		assert_eq!(field.mul(256.into(), Element::one()), Element::from(256));
		assert_eq!(field.mul(Element::one(), 256.into()), Element::from(256));
	}

	#[test]
	fn fft_simple_works() {
		let pd = 512;
		let field = BinaryField::new(1033.into()).unwrap();
		let domain = (0_usize..pd).into_iter().map(Element::from).collect::<Vec<_>>();
		let poly = domain.iter().map(|x| x.pow_mod(9.into(), pd as u64)).collect::<Vec<Element>>();
		let z = _simple_ft(&field, &domain[..], &poly[..]);

		let expected_poly: Vec<_> = vec![0, 1, 0, 227, 0, 357, 0, 327, 0, 329, 0, 299, 0, 429, 0, 143, 0, 145, 0, 371, 0, 501, 0, 471, 0, 473, 0, 443, 0, 61, 0, 287, 0, 289, 0, 3, 0, 133, 0, 103, 0, 105, 0, 75, 0, 205, 0, 431, 0, 433, 0, 147, 0, 277, 0, 247, 0, 249, 0, 219, 0, 349, 0, 63, 0, 65, 0, 291, 0, 421, 0, 391, 0, 393, 0, 363, 0, 493, 0, 207, 0, 209, 0, 435, 0, 53, 0, 23, 0, 25, 0, 507, 0, 125, 0, 351, 0, 353, 0, 67, 0, 197, 0, 167, 0, 169, 0, 139, 0, 269, 0, 495, 0, 497, 0, 211, 0, 341, 0, 311, 0, 313, 0, 283, 0, 413, 0, 127, 0, 129, 0, 355, 0, 485, 0, 455, 0, 457, 0, 427, 0, 45, 0, 271, 0, 273, 0, 499, 0, 117, 0, 87, 0, 89, 0, 59, 0, 189, 0, 415, 0, 417, 0, 131, 0, 261, 0, 231, 0, 233, 0, 203, 0, 333, 0, 47, 0, 49, 0, 275, 0, 405, 0, 375, 0, 377, 0, 347, 0, 477, 0, 191, 0, 193, 0, 419, 0, 37, 0, 7, 0, 9, 0, 491, 0, 109, 0, 335, 0, 337, 0, 51, 0, 181, 0, 151, 0, 153, 0, 123, 0, 253, 0, 479, 0, 481, 0, 195, 0, 325, 0, 295, 0, 297, 0, 267, 0, 397, 0, 111, 0, 113, 0, 339, 0, 469, 0, 439, 0, 441, 0, 411, 0, 29, 0, 255, 0, 257, 0, 483, 0, 101, 0, 71, 0, 73, 0, 43, 0, 173, 0, 399, 0, 401, 0, 115, 0, 245, 0, 215, 0, 217, 0, 187, 0, 317, 0, 31, 0, 33, 0, 259, 0, 389, 0, 359, 0, 361, 0, 331, 0, 461, 0, 175, 0, 177, 0, 403, 0, 21, 0, 503, 0, 505, 0, 475, 0, 93, 0, 319, 0, 321, 0, 35, 0, 165, 0, 135, 0, 137, 0, 107, 0, 237, 0, 463, 0, 465, 0, 179, 0, 309, 0, 279, 0, 281, 0, 251, 0, 381, 0, 95, 0, 97, 0, 323, 0, 453, 0, 423, 0, 425, 0, 395, 0, 13, 0, 239, 0, 241, 0, 467, 0, 85, 0, 55, 0, 57, 0, 27, 0, 157, 0, 383, 0, 385, 0, 99, 0, 229, 0, 199, 0, 201, 0, 171, 0, 301, 0, 15, 0, 17, 0, 243, 0, 373, 0, 343, 0, 345, 0, 315, 0, 445, 0, 159, 0, 161, 0, 387, 0, 5, 0, 487, 0, 489, 0, 459, 0, 77, 0, 303, 0, 305, 0, 19, 0, 149, 0, 119, 0, 121, 0, 91, 0, 221, 0, 447, 0, 449, 0, 163, 0, 293, 0, 263, 0, 265, 0, 235, 0, 365, 0, 79, 0, 81, 0, 307, 0, 437, 0, 407, 0, 409, 0, 379, 0, 509, 0, 223, 0, 225, 0, 451, 0, 69, 0, 39, 0, 41, 0, 11, 0, 141, 0, 367, 0, 369, 0, 83, 0, 213, 0, 183, 0, 185, 0, 155, 0, 285, 0, 511].into_iter().map(Element::from).collect();
		assert_eq!(expected_poly.len(), pd);
		assert_eq!(poly, &expected_poly[0..]);

		let expected_z: Vec<_> = vec![0, 0, 429, 108, 891, 245, 857, 420, 735, 4, 232, 954, 189, 1018, 657, 708, 453, 938, 145, 643, 605, 202, 95, 306, 141, 954, 830, 123, 1012, 429, 339, 702, 579, 352, 169, 363, 929, 167, 961, 270, 154, 854, 787, 275, 346, 534, 165, 462, 780, 360, 538, 656, 65, 396, 588, 642, 292, 725, 84, 704, 833, 278, 634, 563, 563, 813, 44, 960, 768, 905, 936, 97, 102, 59, 635, 246, 209, 713, 1011, 431, 344, 966, 602, 812, 363, 609, 611, 155, 200, 714, 794, 138, 452, 593, 551, 83, 1003, 762, 159, 630, 292, 1017, 789, 791, 487, 1007, 720, 657, 19, 41, 510, 422, 358, 516, 286, 1019, 721, 59, 560, 0, 206, 149, 894, 955, 407, 116, 319, 390, 1020, 402, 767, 91, 711, 555, 477, 891, 509, 346, 993, 787, 871, 843, 833, 427, 183, 484, 434, 531, 942, 103, 513, 997, 396, 243, 749, 835, 729, 575, 206, 284, 597, 73, 159, 684, 910, 561, 650, 191, 182, 924, 354, 56, 854, 11, 584, 819, 726, 972, 875, 661, 981, 1013, 9, 24, 867, 1014, 537, 391, 997, 599, 71, 818, 64, 533, 587, 759, 650, 747, 280, 189, 699, 137, 338, 539, 857, 938, 128, 849, 820, 81, 595, 382, 812, 431, 931, 458, 890, 996, 153, 559, 738, 794, 930, 206, 961, 669, 402, 721, 730, 631, 357, 337, 529, 965, 121, 769, 931, 936, 842, 73, 583, 482, 821, 865, 227, 894, 40, 498, 701, 21, 50, 282, 374, 102, 600, 743, 711, 971, 1021, 238, 241, 264, 160, 911, 346, 422, 32, 340, 892, 929, 562, 535, 804, 626, 825, 855, 287, 833, 757, 306, 646, 172, 2, 246, 317, 1001, 472, 553, 972, 383, 15, 884, 164, 251, 641, 27, 764, 118, 19, 927, 579, 34, 701, 984, 567, 980, 933, 580, 622, 415, 845, 972, 93, 750, 136, 620, 4, 930, 797, 373, 353, 458, 33, 144, 916, 538, 444, 26, 302, 83, 328, 622, 104, 1022, 834, 149, 663, 605, 511, 814, 978, 190, 480, 228, 102, 231, 407, 635, 589, 724, 352, 886, 121, 261, 1018, 509, 387, 514, 532, 255, 926, 111, 129, 408, 957, 18, 390, 209, 324, 339, 1002, 449, 873, 879, 33, 633, 539, 816, 450, 870, 114, 520, 665, 910, 147, 127, 577, 594, 97, 279, 926, 894, 584, 830, 26, 591, 593, 226, 57, 365, 922, 319, 744, 2, 58, 538, 64, 723, 709, 762, 727, 3, 611, 680, 948, 70, 920, 542, 854, 593, 831, 853, 473, 932, 342, 184, 330, 804, 433, 84, 549, 622, 172, 646, 776, 648, 63, 816, 219, 260, 417, 774, 517, 928, 664, 152, 101, 910, 118, 384, 628, 726, 589, 118, 556, 5, 665, 290, 366, 41, 60, 581, 780, 69, 216, 426, 64, 561, 658, 656, 417, 375, 882, 825, 561, 103, 892, 218, 91, 144, 91, 160, 132, 895, 78, 464, 146, 153, 117, 456, 810, 80, 822, 388, 834, 1006, 95, 422, 50, 862, 121, 50, 592, 424, 601, 283, 819, 835, 262, 791, 518, 507].into_iter().map(Element::from).collect();
		assert_eq!(expected_z.len(), z.len());
		assert_slice_eq!(&z[..], &expected_z[0..]);
	}

	#[test]
	fn fft_faster_works() {
		let pd = 1024;
		let field = BinaryField::new(1033.into()).unwrap();
		let domain = (0_usize..pd).into_iter().map(Element::from).collect::<Vec<_>>();
		let poly = domain.iter().map(|x| x.pow_mod(9.into(), pd as u64)).collect::<Vec<Element>>();
		let z = fft(&field, &domain[..], &poly[..]);

		let recovered_poly = invfft(&field, &domain[..], &z[..]);
		assert_eq!(&poly[..], &recovered_poly[..]);
	}

	#[test]
	fn fft_encode_and_recover() {
		let pd = 1024;
		let field = BinaryField::new(1033.into()).unwrap();
		let domain = (0_usize..pd).into_iter().map(Element::from).collect::<Vec<_>>();

		let poly3 = (0..25)
			.into_iter()
			.map(Element::from)
			.map(|x: Element| x.pow_mod(9.into(), pd as u64))
			.collect::<Vec<Element>>();

		let xs = (0..25).into_iter().map(Element::from).map(|x| (x * 11) % 32).collect::<Vec<Element>>();

		let ys = xs.iter().map(|&x| field.eval_poly_at(&poly3, x)).collect::<Vec<Element>>();
		let poly4 = interpolate(&field, &xs[..], &ys[..]);

		assert_eq!(&poly4[..poly3.len()], &poly3[..]);
	}

	#[test]
	fn fft_encode_and_recover_subset() {
		let pd = 1024;
		let field = BinaryField::new(1033.into()).unwrap();
		let domain = (0_usize..pd).into_iter().map(Element::from).collect::<Vec<_>>();

		let poly3 = (0..25)
			.into_iter()
			.map(Element::from)
			.map(|x: Element| x.pow_mod(9.into(), pd as u64))
			.collect::<Vec<Element>>();

		// skip the first 1
		let xs = (1..25).into_iter().map(|x| (x * 11) % 32).map(Element::from).collect::<Vec<Element>>();
		let ys = xs.iter().map(|&x| field.eval_poly_at(&poly3[..], x)).collect::<Vec<_>>();
		let poly5 = interpolate(&field, &xs[..], &ys[..]);

		assert_eq!(&poly5[..poly3.len()], &poly3);
	}

	#[test]
	fn fft_edr() {
		// for GF(2^16)
		// let poly = 0x01_02_10 as u32;

		let field = BinaryField::new(1033.into()).unwrap();

		let pd = 1024;
		println!("S1");
		let poly = (0_usize..pd)
			.into_iter()
			.map(Element::from)
			.map(|x| x.pow_mod(9.into(), pd as u64))
			.collect::<Vec<Element>>();

		println!("S2");
		let domain = (0_usize..pd).into_iter().map(Element::from).collect::<Vec<_>>();
		let z = fft(&field, &domain[..], &poly[..]);
		println!("S3");
		let z2 = _simple_ft(&field, &domain[..], &poly[..]);
		for (idx, (&zv, &z2v)) in z.iter().zip(z2.iter()).enumerate() {
			assert_eq!(
				zv,
				z2v,
				"Failed at idx = {}, {:?} vs {:?}",
				idx,
				&z[idx.saturating_sub(2)..idx.saturating_add(2)],
				&z2[idx.saturating_sub(2)..idx.saturating_add(2)],
			);
		}
		assert_eq!(&z[..], &z2[..]);

		println!("S4");
		let poly2 = invfft(&field, &domain[..], &z[..]);
		assert_eq!(&poly2[..], &poly[..]);

		let poly3 = (0..25)
			.into_iter()
			.map(Element::from)
			.map(|x: Element| x.pow_mod(9.into(), pd as u64))
			.collect::<Vec<Element>>();

		let xs = (0..25).into_iter().map(Element::from).map(|x| (x * 11) % 32).collect::<Vec<Element>>();

		let ys = xs.iter().map(|&x| field.eval_poly_at(&poly3, x)).collect::<Vec<Element>>();
		let poly4 = interpolate(&field, &xs[..], &ys[..]);

		assert_eq!(&poly4[..poly3.len()], &poly3[..]);

		let xs = (1..25).into_iter().map(|x| (x * 11) % 32).map(Element::from).collect::<Vec<Element>>();
		let ys = xs.iter().map(|&x| field.eval_poly_at(&poly3[..], x)).collect::<Vec<_>>();
		let poly5 = interpolate(&field, &xs[..], &ys[..]);

		assert_eq!(&poly5[..poly3.len()], &poly3);
	}
}
