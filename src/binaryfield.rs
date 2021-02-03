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

    fn pow<R>(&self, exp: R) -> Self where R: Into<Self> {
        Self(self.0.overflowing_pow(exp.into().0 as u32).0)
    }
}

// #[target_feature = "step_trait"]
// impl Step for Element {
//     fn steps_between(start: &Self, end: &Self) -> Option<usize> {
//         let val = end.0.saturating_sub(start.0);
//         Some(val as usize + 1_usize)
//     }
// }


impl<R> std::ops::BitXor<R> for Element where R: Into<Element> {
    type Output = Self;
    fn bitxor(self, rhs: R) -> Self::Output {
        Self(self.0 ^ rhs.into().0)
    }
}


impl<R> std::ops::BitXorAssign<R> for Element where R: Into<Element> {
    fn bitxor_assign(&mut self, rhs: R) {
        self.0 ^= rhs.into().0
    }
}

impl<R> std::ops::Mul<R> for Element where R: Into<Element> {
    type Output = Self;
    fn mul(self, rhs: R) -> Self::Output {
        let a = self.0;
        let b = rhs.into().0;
        Self(raw_mul(a,b))
    }
}

impl<R> std::ops::Rem<R> for Element  where R: Into<Element>  {
    type Output = Self;
    fn rem(self, rhs: R) -> Self::Output {
        let a = self.0;
        let b = rhs.into().0;
        Self(raw_mod(a,b))
    }
}


impl<R> std::ops::Add<R> for Element where R: Into<Element> {
    type Output = Self;
    fn add(self, rhs: R) -> Self::Output {
        let rhs = rhs.into();
        Self(self.0 + rhs.0)
    }
}


impl<R> std::ops::Sub<R> for Element where R: Into<Element> {
    type Output = Self;
    fn sub(self, rhs: R) -> Self::Output {
        let rhs = rhs.into();
        Self(self.0 - rhs.0)
    }
}

impl<R> std::ops::Shl<R> for Element where R: Into<Element> {
    type Output = Self;
    fn shl(self, rhs: R) -> Self::Output {
        let rhs = rhs.into();
        Self(self.0 << rhs.0)
    }
}
impl<R> std::ops::Shr<R> for Element where R: Into<Element> {
    type Output = Self;
    fn shr(self, rhs: R) -> Self::Output {
        let rhs = rhs.into();
        Self(self.0 >> rhs.0)
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
    return x > 0_u16 && x & (x-1) == 0
}

#[inline(always)]
fn raw_mul(a: u16, b: u16) -> u16 {
    if a.saturating_mul(b) == 0 {
        return 0
    }
    let mut o = 0;
    for i in 0..(log2(b) + 1) {
        if (b & (1<<i)) != 0x0 {
            o ^= a<<i
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
struct BinaryField{
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
                return Ok(())
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
    fn add(&self, x: Element, y: Element) -> Element{
        x ^ y
    }

    fn sub(&self, x: Element, y: Element) -> Element {
        self.add(x,y)
    }

    fn mul(&self, x: Element, y:Element)  -> Element {
        if x * y == 0_u16 {
            Element::zero()
        } else {
            let idx = self.invcache[x.0 as usize].unwrap() + self.invcache[y.0 as usize].unwrap();
            self.cache[idx]
        }
    }

    fn sqr(&self, x: Element)  -> Element {
        if x == 0_usize {
            Element::zero()
        } else {
            let idx = (self.invcache[x.0 as usize].unwrap() << 1_usize) % self.order.0 as usize;
            self.cache[idx]
        }
    }

    fn div(&self, x: Element, y:Element)  -> Element {
        if x == 0 {
            Element::zero()
        } else {
            let idx: Element = self.order + self.invcache[x.0 as usize].unwrap() - self.invcache[y.0 as usize].unwrap();
            self.cache[idx.0 as usize]
        }
    }

    fn inv(&self, x: Element)  -> Element {
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
            outputs[i] = if value.is_some() {
                 self.mul(partials[i], inv)
                } else  {
                    Element::zero()
                };
            inv = self.mul(inv, value.unwrap_or(Element::one()))
        }
        outputs
    }

    // TODO XXX understand why there are two diff impls
    fn div_(&self, x: Element, y:Element) -> Element {
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
        self.add_polys(a,b)
    }

    fn mul_by_const(&self, a: &[Element], c: Element) -> Vec<Element> {
        a.into_iter().map(move |x| self.mul(*x, c)).collect::<Vec<Element>>()
    }

    fn mul_polys(&self, a: Vec<Element>, b: Vec<Element>) -> Vec<Element> {
        let mut o = vec![Element::zero(); a.len() + b.len() - 1];
        for (i, &aval) in a.iter().enumerate() {
            for (j, &bval) in b.iter().enumerate() {
                o[i+j] ^= self.mul(aval, bval)
            }
        }
        o
    }

    fn div_polys(&self, mut a: Vec<Element>, b: Vec<Element>) -> Vec<Element> {
        assert!(a.len() >= b.len());
        let mut o = vec![];
        let mut apos = a.len() - 1_usize;
        let mut bpos = b.len() - 1_usize;
        let mut diff = apos - bpos;
        while diff >= 0_usize {
            let quot = self.div(a[apos], b[bpos]);
            o.insert(0, quot);
            for (i,b) in (0..bpos).into_iter().enumerate().rev() {
                a[diff+i] ^= self.mul(Element::from(b), quot);
            }
            apos -= 1_usize;
            diff -= 1_usize;
        }
        o
    }

    // Build a polynomial that returns 0 at all specified xs
    fn zpoly(&self, xs: &[Element]) -> Vec<Element> {
        let mut roots = vec![Element::one()];
        for &x in xs {
            roots.insert(0, Element::zero());
            let rl = roots.len()-1;
            for j in 0..rl {
                let val = roots[j+1];
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
        let mut nums = xs.iter().map(|&x| self.div_polys(root.clone(), vec![x, Element::one()]) ).collect::<Vec<Vec<Element>>>();
        // Generate denominators by evaluating numerator polys at each x
        let denoms = xs.iter().zip(nums.iter()).map(|(&x, num)| Some(self.eval_poly_at(num, x))).collect::<Vec<Option<Element>>>();
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
    domain.into_iter().map(|&item| field.eval_poly_at(poly, item)).collect::<Vec<_>>()
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
        return (
            vec![poly[0]],
            vec![if poly.len() == 2 { poly[1] } else { Element::zero() } ])
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
    let mut low_cast = dbg!(cast(field, dbg!(low), k));
    let high_cast = dbg!(cast(field, dbg!(high), k));
    // Combine the results
    (
        { low_cast.0.extend(high_cast.0.into_iter()); low_cast.0 },
        { low_cast.1.extend(high_cast.1.into_iter()); low_cast.1 }
    )
}

// Returns a polynomial p2 such that p2(x) = poly(x**2+kx)
fn compose(field: &BinaryField, poly: &[Element], k: Element) -> Vec<Element> {
    if poly.len() == 2 {
        return vec![poly[0], field.mul(poly[1], k), poly[1], Element::zero()]
    }
    if poly.len() == 1 {
        let mut res = poly.to_vec();
        res.push(Element::zero());
        return res
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
        o[i+mod_power] ^= field.mul(high, k_to_mod_power);
        o[i+2*mod_power] ^= high;
    }
    o
}

// Equivalent to [field.eval_poly_at(poly, x) for x in domain]
// Special thanks to www.math.clemson.edu/~sgao/papers/GM10.pdf for insights
// though this algorithm is not exactly identical to any algorithm in the paper
fn fft(field: &BinaryField, domain: &[Element], poly: &[Element]) -> Vec<Element>{
    println!("fft");
    // Base case: constant polynomials
    // if domain.len() == 1{
    //     return [poly[0]]
    dbg!(domain.len());
    dbg!(poly.len());
    if domain.len() <= 8{
        return _simple_ft(field, domain, poly)
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
    dbg!(cast_domain.len());
    // Two half-size sub-problems over the smaller domain, recovering
    // evaluations of evens && odds over the smaller domain
    let even_points = fft(field, &cast_domain[..], &evens[..]);
    let odd_points = fft(field, &cast_domain[..], &odds[..]);
    // Combine the evaluations of evens && odds into evaluations of poly
    let dl = domain.len() >> 1_usize;
    let mut o = Vec::with_capacity(dl << 1_usize);
    for i in 0..dl {
        o.push(even_points[i] ^ field.mul(domain[i*2], odd_points[i]));
        o.push(even_points[i] ^ field.mul(domain[i*2+1], odd_points[i]));
    }
    o
}



// The inverse function of fft, does the steps backwards
fn invfft(field: &BinaryField, domain: &[Element], vals: &[Element]) -> Vec<Element> {
    // Base case: constant polynomials
    if domain.len() == 1 {
        return vals.to_vec()
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
        let (p_of_x, p_of_x_plus_k) = (vals[i*2], vals[i*2+1]);
        let x = domain[i*2];
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
    (0..vals.len()).into_iter().map(|i| { composed_evens[i] ^ composed_odds[i] } ).collect::<Vec<_>>()
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
	&[16440, 34925, 14360, 22561, 43883, 36645, 7613, 26531, 8597, 59502, 61283, 53412]
    ];

fn invfft2(field: &BinaryField, vals: &[Element]) -> Vec<Element> {
    if vals.len() == 1 {
        return vals.to_vec()
    }
    let len_half = vals.len() >> 1_usize;
    let left = invfft2(field, &vals[..len_half]);

    let tmp = invfft2(field, &vals[len_half..]);
    let right = shift(field, &tmp[..], len_half.into());

    let mut o = vec![Element::zero(); vals.len()];
    for (j, (left, right)) in left.into_iter().zip(right.into_iter()).enumerate() {
        o[j] ^= left;
        for (i, &coeff) in SHIFT_POLYS[log2(vals.len() as u16) as usize].into_iter().enumerate() {
            o[(1<<i)+j] ^= field.mul(left ^ right, coeff.into());
        }
    }
    o
}

// Multiplies two polynomials using the FFT method
fn mul(field: &BinaryField, domain: &[Element], p1: &[Element], p2: &[Element]) -> Vec<Element> {
    assert!(p1.len() <= domain.len() && p2.len() <= domain.len());
    let values1 = fft(field, domain, p1);
    let values2 = fft(field, domain, p2);
    let values3 = values1.into_iter().zip(values2.into_iter()).map(|(v1,v2)| field.mul(v1, v2)).collect::<Vec<Element>>();
    invfft(field, domain, &values3[..])
}

// Generates the polynomial `p(x) = (x - xs[0]) * (x - xs[1]) * ...`
fn zpoly(field: &BinaryField, xs: Vec<Element>) -> Vec<Element> {
    if xs.is_empty() {
        return vec![Element::one()]
    }
    if xs.len() == 1 {
        return vec![xs[0], Element::one()]
    }
    let domain = (0_u16..(2_u16 << log2(xs.iter().max().unwrap().0) + 1 )).into_iter().map(Element::from).collect::<Vec<_>>();
    let offset = domain[1];
    let z_left = zpoly(field, xs.iter().step_by(2).copied().collect());
    let z_right = zpoly(field, xs.iter().skip(1).step_by(2).copied().collect());
    mul(field, &domain[..], &z_left, &z_right)
}

// Returns q(x) = p(x + k)
fn shift(field: &BinaryField, poly: &[Element], k: Element) -> Vec<Element> {
    if poly.len() == 1{
        return poly.to_vec()
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
fn interpolate(field: &BinaryField, xs: &[Element], vals: &[Element]) -> Vec<Element>
{
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
    let shifted_p_times_z_values = &shifted_p_times_z_values[domain_size.0 as usize ..];
    // print("shifted_p_times_z_values =", shifted_p_times_z_values)
    let shifted_p_values = shifted_p_times_z_values
        .into_iter()
        .zip(z_values[domain_size.0 as usize ..].into_iter())
        .map(|(&x,&y)| field.div(x, y) )
        .collect::<Vec<Element>>();
    // print("shifted_p_values =", shifted_p_values)
    let shifted_p = invfft(field, &domain, &shifted_p_values);
    shift(field, &shifted_p, domain_size)
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fft_edr() {
        // for GF(2^16)
        let poly = 0x01_02_10 as u32;



        let field = BinaryField::new(1033.into()).unwrap();

        let pd = 1024;
        println!("S1");
        let poly = (0_usize..pd).into_iter().map(Element::from).map(|x| x.pow(9) % pd).collect::<Vec<Element>>();

        println!("S2");
        let domain = (0_usize..pd).into_iter().map(Element::from).collect::<Vec<_>>();
        let z = fft(&field, &domain[..], &poly[..]);
        println!("S3");
        let z2 = _simple_ft(&field, &domain[..], &poly[..]);
        assert_eq!(&z[..], &z2);

        println!("S4");
        let poly2 = invfft(&field, &domain[..], &z[..]);
        assert_eq!(&poly2[..], &poly[..]);


        let poly3 = (0..25).into_iter().map(Element::from).map(|x: Element| x.pow(9) % pd).collect::<Vec<Element>>();

        let xs = (0..25).into_iter().map(Element::from).map(|x| (x * 11) % 32).collect::<Vec<Element>>();

        let ys = xs.iter().map(|&x| {
            field.eval_poly_at(&poly3, x)
        }).collect::<Vec<Element>>();
        let poly4 = interpolate(&field, &xs[..], &ys[..]);

        assert_eq!(&poly4[..poly3.len()], &poly3[..]);

        let xs = (1..25).into_iter().map(|x| (x * 11) % 32).map(Element::from).collect::<Vec<Element>>();
        let ys = xs.iter().map(|&x| { field.eval_poly_at(&poly3[..], x) }).collect::<Vec<_>>();
        let poly5 = interpolate(&field, &xs[..], &ys[..]);

        assert_eq!(&poly5[..poly3.len()], &poly3);
    }
}
