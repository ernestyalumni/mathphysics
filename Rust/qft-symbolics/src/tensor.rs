//! Strict abstract-index algebra for commuting tensors and a constant symmetric
//! Lorentz metric. Independent bounded orbit enumeration, not xPerm or a CAS.
//! Every raw term is validated before zero coefficients or contractions vanish.
//! Spinor epsilon is an ordinary antisymmetric tensor, never a symmetric metric.
use crate::{polynomial::Poly, scalar::Q};
use std::{
    collections::{BTreeMap, BTreeSet, VecDeque},
    fmt,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Space {
    Lorentz,
    Undotted,
    Dotted,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Variance {
    Up,
    Down,
}
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum Label {
    Named(String),
    Bound(u16),
}
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Index {
    space: Space,
    label: Label,
    variance: Variance,
}
impl Index {
    pub fn new(space: Space, name: impl Into<String>, variance: Variance) -> Self {
        Self {
            space,
            label: Label::Named(name.into()),
            variance,
        }
    }
    pub fn up(space: Space, name: impl Into<String>) -> Self {
        Self::new(space, name, Variance::Up)
    }
    pub fn down(space: Space, name: impl Into<String>) -> Self {
        Self::new(space, name, Variance::Down)
    }
    pub fn space(&self) -> Space {
        self.space
    }
    pub fn variance(&self) -> Variance {
        self.variance
    }
    pub fn name(&self) -> Option<&str> {
        if let Label::Named(n) = &self.label {
            Some(n)
        } else {
            None
        }
    }
    fn key(&self) -> (Space, Label) {
        (self.space, self.label.clone())
    }
}
/// Lorentz slots support all metric-related variance views of the same tensor.
/// Spinor slots have fixed variance; epsilon raising/lowering is not implemented.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Slot {
    Lorentz,
    FixedSpinor(Space, Variance),
}
impl Slot {
    fn accepts(&self, i: &Index) -> bool {
        match self {
            Self::Lorentz => i.space == Space::Lorentz,
            Self::FixedSpinor(s, v) => *s != Space::Lorentz && i.space == *s && i.variance == *v,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Dimension {
    Fixed(u16),
    Symbol(String),
}
impl Dimension {
    fn poly(&self) -> Poly {
        match self {
            Self::Fixed(n) => Poly::constant(*n as i64),
            Self::Symbol(n) => Poly::var(n),
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Symmetry {
    None,
    Symmetric,
    Antisymmetric,
    /// output slot i receives input slot permutation[i]; sign must be +/-1.
    Generators(Vec<(Vec<usize>, i8)>),
}
#[derive(Clone, Debug, PartialEq, Eq)]
struct TensorSpec {
    slots: Vec<Slot>,
    group: Vec<(Vec<usize>, i8)>,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Context {
    tensors: BTreeMap<String, TensorSpec>,
    dimensions: BTreeMap<Space, Dimension>,
    budget: usize,
}
impl Default for Context {
    fn default() -> Self {
        Self {
            tensors: BTreeMap::new(),
            dimensions: BTreeMap::from([
                (Space::Lorentz, Dimension::Fixed(4)),
                (Space::Undotted, Dimension::Fixed(2)),
                (Space::Dotted, Dimension::Fixed(2)),
            ]),
            budget: 100_000,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Error {
    Invalid(String),
    SignatureMismatch,
    ContextMismatch,
    Unsupported(String),
    Budget,
}
impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Invalid(s) => write!(f, "invalid indices/declaration: {s}"),
            Self::SignatureMismatch => write!(f, "free-index signatures differ"),
            Self::ContextMismatch => write!(f, "tensor contexts differ"),
            Self::Unsupported(s) => write!(f, "unsupported: {s}"),
            Self::Budget => write!(f, "tensor orbit/work budget exhausted"),
        }
    }
}
impl std::error::Error for Error {}
type Result<T> = std::result::Result<T, Error>;
struct Budget(usize);
impl Budget {
    fn charge(&mut self, n: usize) -> Result<()> {
        self.0 = self.0.checked_sub(n).ok_or(Error::Budget)?;
        Ok(())
    }
    fn step(&mut self) -> Result<()> {
        if self.0 == 0 {
            Err(Error::Budget)
        } else {
            self.0 -= 1;
            Ok(())
        }
    }
}
impl Context {
    pub fn with_budget(mut self, steps: usize) -> Self {
        self.budget = steps;
        self
    }
    pub fn with_lorentz_dimension(mut self, dimension: Dimension) -> Result<Self> {
        if matches!(dimension, Dimension::Fixed(0))
            || matches!(&dimension,Dimension::Symbol(s) if s.is_empty())
        {
            return Err(Error::Invalid("empty/zero dimension".into()));
        }
        self.dimensions.insert(Space::Lorentz, dimension);
        Ok(self)
    }
    pub fn declare(
        &mut self,
        name: impl Into<String>,
        slots: Vec<Slot>,
        symmetry: Symmetry,
    ) -> Result<()> {
        let name = name.into();
        if name.is_empty() || self.tensors.contains_key(&name) {
            return Err(Error::Invalid("empty or duplicate tensor symbol".into()));
        }
        if slots.len() > 8 {
            return Err(Error::Unsupported("tensor rank >8".into()));
        }
        if slots
            .iter()
            .any(|s| matches!(s, Slot::FixedSpinor(Space::Lorentz, _)))
        {
            return Err(Error::Invalid("use metric-flexible Lorentz slots".into()));
        }
        let n = slots.len();
        let identity: Vec<_> = (0..n).collect();
        let generators = match symmetry {
            Symmetry::None => vec![],
            Symmetry::Generators(g) => g,
            Symmetry::Symmetric | Symmetry::Antisymmetric => {
                let sign = if matches!(symmetry, Symmetry::Antisymmetric) {
                    -1
                } else {
                    1
                };
                (1..n)
                    .map(|i| {
                        let mut p = identity.clone();
                        p.swap(i - 1, i);
                        (p, sign)
                    })
                    .collect()
            }
        };
        for (p, sign) in &generators {
            if ![-1, 1].contains(sign)
                || p.len() != n
                || p.iter().copied().collect::<BTreeSet<_>>() != identity.iter().copied().collect()
            {
                return Err(Error::Invalid("not a signed slot permutation".into()));
            }
            if (0..n).any(|i| slots[i] != slots[p[i]]) {
                return Err(Error::Invalid(
                    "symmetry exchanges incompatible slots".into(),
                ));
            }
        }
        let mut seen = BTreeSet::from([(identity.clone(), 1)]);
        let mut queue = VecDeque::from([(identity, 1)]);
        let mut budget = Budget(self.budget);
        while let Some((p, s)) = queue.pop_front() {
            budget.step()?;
            for (g, t) in &generators {
                budget.step()?;
                let q = (0..n).map(|i| p[g[i]]).collect();
                let state = (q, s * t);
                if seen.insert(state.clone()) {
                    queue.push_back(state);
                }
            }
        }
        self.tensors.insert(
            name,
            TensorSpec {
                slots,
                group: seen.into_iter().collect(),
            },
        );
        Ok(())
    }
    /// This boundary is deliberately unavailable for antisymmetric spinor forms.
    pub fn require_symmetric_metric(&self, space: Space) -> Result<()> {
        if space != Space::Lorentz {
            return Err(Error::Unsupported("spinor epsilon raising/lowering".into()));
        }
        Ok(())
    }
}

#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Factor {
    Tensor(String, Vec<Index>),
    Metric([Index; 2]),
    Delta([Index; 2]),
}
impl Factor {
    pub fn tensor(name: impl Into<String>, indices: Vec<Index>) -> Self {
        Self::Tensor(name.into(), indices)
    }
    pub fn metric(a: Index, b: Index) -> Self {
        Self::Metric([a, b])
    }
    /// Slots must be upper, lower in that order, in the same space.
    pub fn delta(upper: Index, lower: Index) -> Self {
        Self::Delta([upper, lower])
    }
    fn indices(&self) -> &[Index] {
        match self {
            Self::Tensor(_, v) => v,
            Self::Metric(v) | Self::Delta(v) => v,
        }
    }
    fn indices_mut(&mut self) -> &mut [Index] {
        match self {
            Self::Tensor(_, v) => v,
            Self::Metric(v) | Self::Delta(v) => v,
        }
    }
    fn validate(&self, ctx: &Context) -> Result<()> {
        for i in self.indices() {
            if matches!(&i.label,Label::Named(s) if s.is_empty()) {
                return Err(Error::Invalid("empty index label".into()));
            }
        }
        match self {
            Self::Tensor(name, indices) => {
                let spec = ctx
                    .tensors
                    .get(name)
                    .ok_or_else(|| Error::Invalid(format!("undeclared tensor {name}")))?;
                if indices.len() != spec.slots.len()
                    || indices.iter().zip(&spec.slots).any(|(i, s)| !s.accepts(i))
                {
                    return Err(Error::Invalid(format!("slot types for {name}")));
                }
            }
            Self::Metric([a, b]) => {
                ctx.require_symmetric_metric(a.space)?;
                if b.space != a.space || a.variance != b.variance {
                    return Err(Error::Invalid(
                        "metric requires two equal-variance Lorentz slots".into(),
                    ));
                }
            }
            Self::Delta([a, b]) => {
                if a.space != b.space || a.variance != Variance::Up || b.variance != Variance::Down
                {
                    return Err(Error::Invalid(
                        "delta requires typed upper/lower slots".into(),
                    ));
                }
            }
        }
        Ok(())
    }
    // A Lorentz contraction can turn g into delta or vice versa.
    fn normalise_identity(&mut self) {
        if let Self::Metric([a, b]) | Self::Delta([a, b]) = self {
            let (a, b) = (a.clone(), b.clone());
            *self = if a.variance == b.variance {
                Self::Metric([a, b])
            } else if a.variance == Variance::Up {
                Self::Delta([a, b])
            } else {
                Self::Delta([b, a])
            };
        }
    }
    fn variants(&self, ctx: &Context, budget: &mut Budget) -> Result<Vec<(Self, i8)>> {
        // Charge before eagerly cloning a potentially large symmetry orbit.
        let n = match self {
            Self::Tensor(name, _) => ctx.tensors[name].group.len(),
            Self::Metric(_) => 2,
            Self::Delta(_) => 1,
        };
        budget.charge(n)?;
        Ok(match self {
            Self::Tensor(name, indices) => ctx.tensors[name]
                .group
                .iter()
                .map(|(p, s)| {
                    (
                        Self::Tensor(
                            name.clone(),
                            p.iter().map(|i| indices[*i].clone()).collect(),
                        ),
                        *s,
                    )
                })
                .collect(),
            Self::Metric([a, b]) => {
                vec![(self.clone(), 1), (Self::Metric([b.clone(), a.clone()]), 1)]
            }
            Self::Delta(_) => vec![(self.clone(), 1)],
        })
    }
}
#[derive(Clone, Debug)]
pub struct Term {
    pub coefficient: Poly,
    pub factors: Vec<Factor>,
}
impl Term {
    pub fn new(coefficient: Poly, factors: Vec<Factor>) -> Self {
        Self {
            coefficient,
            factors,
        }
    }
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Expr {
    context: Context,
    signature: Vec<Index>,
    terms: BTreeMap<Vec<Factor>, Poly>,
}
#[derive(Default)]
struct Classification {
    signature: Vec<Index>,
    dummy: BTreeSet<(Space, Label)>,
}
fn classify(factors: &[Factor], ctx: &Context) -> Result<Classification> {
    if factors.len() > 8 || factors.iter().map(|f| f.indices().len()).sum::<usize>() > 20 {
        return Err(Error::Unsupported(
            "at most 8 commuting factors and 20 slots".into(),
        ));
    }
    let mut occurrences: BTreeMap<_, Vec<_>> = BTreeMap::new();
    for f in factors {
        f.validate(ctx)?;
        for i in f.indices() {
            occurrences.entry(i.key()).or_default().push(i.clone());
        }
    }
    classify_occurrences(occurrences)
}
fn classify_occurrences(
    occurrences: BTreeMap<(Space, Label), Vec<Index>>,
) -> Result<Classification> {
    let mut result = Classification::default();
    for (key, indices) in occurrences {
        match indices.as_slice() {
            [i] => {
                if matches!(i.label, Label::Bound(_)) {
                    return Err(Error::Invalid("dangling internal dummy".into()));
                }
                result.signature.push(i.clone());
            }
            [a, b] if a.variance != b.variance => {
                result.dummy.insert(key);
            }
            [_, _] => return Err(Error::Invalid("same-variance repeated index".into())),
            _ => return Err(Error::Invalid("index occurs more than twice".into())),
        }
    }
    result.signature.sort();
    Ok(result)
}
fn signature_product(a: &[Index], b: &[Index]) -> Result<Vec<Index>> {
    let mut occurrences: BTreeMap<_, Vec<_>> = BTreeMap::new();
    for i in a.iter().chain(b) {
        occurrences.entry(i.key()).or_default().push(i.clone());
    }
    Ok(classify_occurrences(occurrences)?.signature)
}
fn contract(term: &mut Term, ctx: &Context, budget: &mut Budget) -> Result<()> {
    loop {
        let mut change = false;
        'search: for i in 0..term.factors.len() {
            budget.step()?;
            let (a, b) = match &term.factors[i] {
                Factor::Metric([a, b]) | Factor::Delta([a, b]) => (a.clone(), b.clone()),
                _ => continue,
            };
            if matches!(term.factors[i], Factor::Delta(_)) && a.key() == b.key() {
                term.coefficient = term.coefficient.clone() * ctx.dimensions[&a.space].poly();
                term.factors.remove(i);
                change = true;
                break;
            }
            for (old, new) in [(&a, &b), (&b, &a)] {
                for j in 0..term.factors.len() {
                    if i == j {
                        continue;
                    }
                    if let Some(k) = term.factors[j]
                        .indices()
                        .iter()
                        .position(|v| v.key() == old.key())
                    {
                        term.factors[j].indices_mut()[k] = new.clone();
                        term.factors[j].normalise_identity();
                        term.factors.remove(i);
                        change = true;
                        break 'search;
                    }
                }
            }
        }
        if !change {
            return Ok(());
        }
        classify(&term.factors, ctx)?;
    }
}
fn alpha_normalise(factors: &[Factor], dummy: &BTreeSet<(Space, Label)>) -> Vec<Factor> {
    let mut names = BTreeMap::new();
    let mut next: BTreeMap<Space, u16> = BTreeMap::new();
    let mut out = factors.to_vec();
    for f in &mut out {
        for i in f.indices_mut() {
            let key = i.key();
            if dummy.contains(&key) {
                let id = *names.entry(key).or_insert_with(|| {
                    let n = next.entry(i.space).or_default();
                    let id = *n;
                    *n += 1;
                    id
                });
                i.label = Label::Bound(id);
            }
        }
    }
    out
}
struct Orbit<'a> {
    variants: Vec<Vec<(Factor, i8)>>,
    dummy: &'a BTreeSet<(Space, Label)>,
    best: Option<Vec<Factor>>,
    signs: BTreeSet<i8>,
}
impl Orbit<'_> {
    fn enumerate(
        &mut self,
        used: u16,
        ordered: &mut Vec<Factor>,
        sign: i8,
        budget: &mut Budget,
    ) -> Result<()> {
        budget.step()?;
        if ordered.len() == self.variants.len() {
            // Symmetric metric compatibility also identifies simultaneous
            // raising/lowering of the two ends of a Lorentz dummy contraction.
            // Never apply this sign-free operation to spinor epsilon spaces.
            let lorentz: Vec<_> = self
                .dummy
                .iter()
                .filter(|(s, _)| *s == Space::Lorentz)
                .collect();
            for mask in 0usize..(1usize << lorentz.len()) {
                budget.step()?;
                let mut oriented = ordered.clone();
                for f in &mut oriented {
                    for index in f.indices_mut() {
                        if let Some(k) = lorentz.iter().position(|key| **key == index.key()) {
                            if mask & (1 << k) != 0 {
                                index.variance = match index.variance {
                                    Variance::Up => Variance::Down,
                                    Variance::Down => Variance::Up,
                                };
                            }
                        }
                    }
                }
                let key = alpha_normalise(&oriented, self.dummy);
                match &self.best {
                    None => {
                        self.best = Some(key);
                        self.signs = BTreeSet::from([sign]);
                    }
                    Some(best) if key < *best => {
                        self.best = Some(key);
                        self.signs = BTreeSet::from([sign]);
                    }
                    Some(best) if key == *best => {
                        self.signs.insert(sign);
                    }
                    _ => {}
                }
            }
            return Ok(());
        }
        for i in 0..self.variants.len() {
            if used & (1 << i) != 0 {
                continue;
            }
            for j in 0..self.variants[i].len() {
                let (factor, s) = self.variants[i][j].clone();
                ordered.push(factor);
                self.enumerate(used | (1 << i), ordered, sign * s, budget)?;
                ordered.pop();
            }
        }
        Ok(())
    }
}
fn canonicalise(term: &Term, ctx: &Context, budget: &mut Budget) -> Result<(Vec<Factor>, i8)> {
    let classification = classify(&term.factors, ctx)?;
    let mut orbit = Orbit {
        variants: term
            .factors
            .iter()
            .map(|f| f.variants(ctx, budget))
            .collect::<Result<Vec<_>>>()?,
        dummy: &classification.dummy,
        best: None,
        signs: BTreeSet::new(),
    };
    orbit.enumerate(0, &mut vec![], 1, budget)?;
    let sign = if orbit.signs.len() == 2 {
        0
    } else {
        *orbit.signs.first().expect("orbit includes identity")
    };
    Ok((orbit.best.expect("orbit includes identity"), sign))
}
impl Expr {
    pub fn from_terms(ctx: &Context, terms: Vec<Term>) -> Result<Self> {
        let mut signature = None;
        let mut out = BTreeMap::<Vec<Factor>, Poly>::new();
        let mut budget = Budget(ctx.budget);
        for mut term in terms {
            budget.step()?;
            let sig = classify(&term.factors, ctx)?.signature;
            if signature.as_ref().is_some_and(|s| *s != sig) {
                return Err(Error::SignatureMismatch);
            }
            signature = Some(sig.clone());
            // Only after type validation may a coefficient-zero term disappear.
            if term.coefficient.is_zero() {
                continue;
            }
            contract(&mut term, ctx, &mut budget)?;
            if classify(&term.factors, ctx)?.signature != sig {
                return Err(Error::Invalid("contraction changed free signature".into()));
            }
            let (key, sign) = canonicalise(&term, ctx, &mut budget)?;
            let value = out.remove(&key).unwrap_or_default()
                + term.coefficient * Poly::constant(sign as i64);
            if !value.is_zero() {
                out.insert(key, value);
            }
        }
        Ok(Self {
            context: ctx.clone(),
            signature: signature.unwrap_or_default(),
            terms: out,
        })
    }
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }
    pub fn signature(&self) -> &[Index] {
        &self.signature
    }
    pub fn term_count(&self) -> usize {
        self.terms.len()
    }
    /// Read-only canonical monomials. Re-ingestion is validated like raw input.
    pub fn monomials(&self) -> impl Iterator<Item = (&[Factor], &Poly)> {
        self.terms.iter().map(|(f, p)| (f.as_slice(), p))
    }
    pub fn scalar(&self) -> Option<Poly> {
        if !self.signature.is_empty() {
            return None;
        }
        if self.terms.is_empty() {
            Some(Poly::default())
        } else if self.terms.len() == 1 {
            self.terms.get(&vec![]).cloned()
        } else {
            None
        }
    }
    pub fn negate(&self) -> Self {
        let mut out = self.clone();
        for c in out.terms.values_mut() {
            *c = -c.clone();
        }
        out
    }
    pub fn add(&self, b: &Self) -> Result<Self> {
        if self.context != b.context {
            return Err(Error::ContextMismatch);
        }
        if self.signature != b.signature {
            return Err(Error::SignatureMismatch);
        }
        let mut out = self.clone();
        for (key, p) in &b.terms {
            let sum = out.terms.remove(key).unwrap_or_default() + p.clone();
            if !sum.is_zero() {
                out.terms.insert(key.clone(), sum);
            }
        }
        Ok(out)
    }
    pub fn residual(&self, b: &Self) -> Result<Self> {
        self.add(&b.negate())
    }
    pub fn multiply(&self, b: &Self) -> Result<Self> {
        if self.context != b.context {
            return Err(Error::ContextMismatch);
        }
        let signature = signature_product(&self.signature, &b.signature)?;
        if self.is_zero() || b.is_zero() {
            return Ok(Self {
                context: self.context.clone(),
                signature,
                terms: BTreeMap::new(),
            });
        }
        if self
            .terms
            .len()
            .checked_mul(b.terms.len())
            .is_none_or(|n| n > 1024)
        {
            return Err(Error::Budget);
        }
        let mut terms = vec![];
        for (left, p) in &self.terms {
            for (right, q) in &b.terms {
                let mut factors = left.clone();
                let mut fresh = right.clone();
                // Bound and free labels are distinct types. Shift all RHS binders,
                // so independently constructed contractions cannot capture each other.
                let offset = left
                    .iter()
                    .flat_map(Factor::indices)
                    .filter_map(|i| {
                        if let Label::Bound(n) = i.label {
                            Some(n + 1)
                        } else {
                            None
                        }
                    })
                    .max()
                    .unwrap_or(0);
                for f in &mut fresh {
                    for i in f.indices_mut() {
                        if let Label::Bound(n) = i.label {
                            i.label = Label::Bound(n + offset);
                        }
                    }
                }
                factors.extend(fresh);
                terms.push(Term::new(p.clone() * q.clone(), factors));
            }
        }
        let out = Self::from_terms(&self.context, terms)?;
        if out.signature != signature {
            return Err(Error::SignatureMismatch);
        }
        Ok(out)
    }
    pub fn scale(&self, c: Q) -> Self {
        let mut out = self.clone();
        for p in out.terms.values_mut() {
            *p = p.clone() * Poly::constant(c);
        }
        out.terms.retain(|_, p| !p.is_zero());
        out
    }
}
