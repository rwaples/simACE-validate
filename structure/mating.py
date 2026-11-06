"""Layer 3b: simACE's `standard` mating rules, read off one recorded pedigree.

simACE's standard model (simACE methods.md, "Reproduction" and "Monozygotic
twins"; v0.3.0 `simulate.py`, `_mating_standard` and `place_twins`) makes
each generation of N children from the previous one:

- R1: every parent draws ZTP(lambda) mating slots, the sexes independently.
- R2: the sex with more slots keeps a uniform random subset of them, so both
  have T = min(total female, total male) slots.
- R3: female slots in order, male slots shuffled, paired by position into T
  matings. A duplicate couple has its father swapped with a random other
  mating, for at most 5 passes.
- R4: t ~ Bin(N, p_mz / (1 + p_mz)) MZ pairs, capped at N // 2, leaving
  B = N - t births; (b_1..b_T) ~ Multinomial(B, 1/T) births per mating.
- R5: t of the B births, uniformly without replacement, are twin births. A
  mating's children are its births plus one co-twin per twin birth.
- R6: sex is Bernoulli(1/2); an MZ co-twin copies the first twin's sex.
- R7: each couple's children have consecutive ids in birth order, so
  co-twins have consecutive ids, and couples are in mother-id order.

A mating with no births leaves no trace, so a generation shows its couples
(GLOSSARY): K couples with c_k >= 1 children each. Generation 0's parents are
erased (ADR 0001 §3), so only generations g >= 1 are read.

`mating_stats` gives per-generation counts for the SLiM comparison (Part B).
`mating_law` scores the same counts against each rule's own null (Part A):

| row | statistic of generation g | null given everything before g |
|---|---|---|
| sex_males | males, one co-twin of each MZ pair left out | Bin(N - t_g, 1/2), exact (R6) |
| mz_pairs | MZ pairs t_g | Bin(N, p_mz / (1 + p_mz)) capped at N // 2, exact (R4) |
| alloc_sumsq | sum of c_k^2 | B births surjected onto the K couples, t twin births among them (R4, R5) |
| alloc_multi_mother | mean c_k over couples of mothers with >= 2 couples | the same draws (R4, R5) |
| n_couples | K | ZTP slots, trim, matings occupied by B births (R1, R2, R4) |
| n_mothers_multi | mothers with >= 2 couples | the same draws (R1, R2, R4) |
| n_fathers_multi | fathers with >= 2 couples | the same draws (R1, R2, R4) |

Why each null is exact:

- Sex is drawn apart from mating, and a co-twin copies its twin, so the other
  children are iid Bernoulli(1/2). t_g is drawn from N alone, after pairing
  and independently of it.
- Multinomial cells are exchangeable. Given t (so B) and which K matings are
  occupied, their births are B balls on K boxes conditioned on none empty,
  whatever T was: the allocation rows do not depend on lambda. iid ZTP(mu)
  counts conditioned on summing to B have that law for any mu, and
  `surjection_draws` samples it exactly by halving. Twin births are a
  uniform t-subset of the B births, so given the births a couple's twin births are multivariate
  hypergeometric, and c_k = b_k + (its twin births).
- A person's couples are their occupied slots, and a mating is occupied when
  it has a birth (twin births fall on births). The duplicate swap exchanges
  fathers between matings, so it keeps every person's slot count, and who
  pairs with whom does not enter these rows.

Conditioning. A generation's children are the next generation's parents, so
summed statistics would condition on earlier outcomes. Instead each row takes
one randomised PIT u_g per generation, under the law of generation g given
everything before it. Mating depends on the past only through the parents'
sexes, and each row's other conditioning (t_g, the couples) is drawn apart
from the row's own randomness. So u_1..u_{G-1} are iid U(0, 1), and their
Stouffer combination is one exactly uniform PIT per pedigree, scored over
seeds by calibration.calibrate (ADR 0001 §2).

The allocation and partner-count nulls assume no duplicate couple survives
the 5 swap passes; sex_males and mz_pairs do not depend on pairing.
tools/residual_duplicates.py measures how rare a survivor is. A surviving duplicate
is kept as a second mating next to the first, so it shows as one couple whose
children came from two multinomial cells.

R6's shared sex and parents and R7's order are deterministic: `generations`
raises on a violation, which makes the case an ERROR row.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy import optimize, special, stats

from calibration import randomised_pit, randomised_pit_discrete

if TYPE_CHECKING:
    import polars as pl

# Power mutations (power.py): each changes what simACE is compared with, never
# simACE. The values sit here once; the runners and the descriptions read them.
LAMBDA_ALT = 0.6  # lambda_06
PMZ_ALT = 0.03  # pmz_03
P_MALE_ALT = 0.52  # sex_052
DIRICHLET_SHAPE = 10.0  # alloc_overdispersed

# std_mating (oracle mating_law): mutations of the per-rule nulls.
MUTATIONS = {
    "null": "simACE's statistic replaced by one more draw from its null (checks the scoring only)",
    "null_slim": "Part A scored on slim/standard.slim pedigrees instead of simACE's",
    "lambda_06": f"partner-count null drawn with lambda = {LAMBDA_ALT} instead of the case's",
    "pmz_03": f"MZ null with p_mz = {PMZ_ALT} instead of the case's",
    "sex_052": f"sex null with P(male) = {P_MALE_ALT}",
    "alloc_overdispersed": (
        f"allocation null with Dirichlet({DIRICHLET_SHAPE:g})-multinomial births, conditioned on no empty couple"
    ),
    "alloc_by_mother": "allocation null with couple weights 1/(mother's couples)",
}

# std_structure (oracle slim_standard): (effect, standard.slim defines replaced).
SLIM_MUTATIONS: dict[str, tuple[str, dict[str, float]]] = {
    "null": ("simACE replaced by a second standard.slim run (false-positive rate)", {}),
    "lambda_06": (f"SLiM with LAMBDA = {LAMBDA_ALT} instead of the case's", {"LAMBDA": LAMBDA_ALT}),
    "pmz_03": (f"SLiM with PMZ = {PMZ_ALT} instead of the case's", {"PMZ": PMZ_ALT}),
}

# Mutations whose runner never reads simACE's pedigrees, so run_case skips simulating them.
REPLACES_SIMACE = {"mating_law": {"null_slim"}, "slim_standard": {"null"}}


@dataclass(frozen=True)
class Generation:
    """What generation g >= 1 shows of the mating that made it."""

    n_female: int  # parents: generation g - 1
    n_male: int
    n_kids: int
    n_mz: int  # MZ pairs, t_g
    males: int  # male children, one co-twin of each MZ pair left out
    counts: np.ndarray  # children per couple
    couple_mother: np.ndarray
    couple_father: np.ndarray
    mother_couples: np.ndarray  # per couple: how many couples its mother has


def n_multi(owners: np.ndarray) -> int:
    """Return how many values occur at least twice."""
    return int(np.sum(np.unique(owners, return_counts=True)[1] >= 2))


def occurrences(owners: np.ndarray) -> np.ndarray:
    """Return, for each element, how many times its value occurs."""
    _, inverse, counts = np.unique(owners, return_inverse=True, return_counts=True)
    return counts[inverse]


def check_twins(g: int, kids: np.ndarray, gen: np.ndarray, sex: np.ndarray, twin: np.ndarray) -> np.ndarray:
    """Check that generation g's MZ links are reciprocal, consecutive and same-sex; return each pair's lower id."""
    first = kids[twin[kids] > kids]
    second = twin[first]
    if 2 * len(first) != np.sum(twin[kids] >= 0) or np.any(gen[second] != g) or np.any(twin[second] != first):
        raise ValueError(f"generation {g}: an MZ link is not reciprocal within the generation")
    if np.any(second != first + 1):
        raise ValueError(f"generation {g}: co-twins do not have consecutive ids")
    if np.any(sex[first] != sex[second]):
        raise ValueError(f"generation {g}: an MZ pair differs in sex")
    return first


def generations(ped: pl.DataFrame) -> list[Generation]:
    """Read every generation g >= 1, checking parents, the MZ invariants and R7's order; rows must be in id order.

    Recorded generation 0 has no parents, so only its MZ links are checked.
    """
    gen = ped["generation"].to_numpy()
    sex = ped["sex"].to_numpy().astype(np.int64)
    mother = ped["mother"].to_numpy().astype(np.int64)
    father = ped["father"].to_numpy().astype(np.int64)
    twin = ped["twin"].to_numpy().astype(np.int64)
    check_twins(0, np.flatnonzero(gen == 0), gen, sex, twin)
    out = []
    for g in range(1, int(gen.max()) + 1):
        kids = np.flatnonzero(gen == g)
        if len(kids) == 0:
            raise ValueError(f"generation {g}: no children")
        prev = np.flatnonzero(gen == g - 1)
        for label, parents, want in (("mother", mother[kids], 0), ("father", father[kids], 1)):
            if np.any(parents < 0) or np.any(gen[parents] != g - 1) or np.any(sex[parents] != want):
                raise ValueError(f"generation {g}: a {label} is not a generation-{g - 1} individual of that sex")
        mother_steps = np.diff(mother[kids])
        if np.any(mother_steps < 0):
            raise ValueError(f"generation {g}: mothers are not in id order")
        keys, counts = np.unique(np.stack([mother[kids], father[kids]]), axis=1, return_counts=True)
        runs = 1 + np.count_nonzero(mother_steps | np.diff(father[kids]))
        if runs != len(counts):
            raise ValueError(f"generation {g}: a couple's children are not one block of consecutive ids")
        first = check_twins(g, kids, gen, sex, twin)
        second = twin[first]
        if np.any(mother[first] != mother[second]) or np.any(father[first] != father[second]):
            raise ValueError(f"generation {g}: an MZ pair has different parents")
        out.append(
            Generation(
                n_female=int(np.sum(sex[prev] == 0)),
                n_male=int(np.sum(sex[prev] == 1)),
                n_kids=len(kids),
                n_mz=len(first),
                males=int(np.sum(sex[kids] == 1) - np.sum(sex[second] == 1)),
                counts=counts,
                couple_mother=keys[0],
                couple_father=keys[1],
                mother_couples=occurrences(keys[0]),
            )
        )
    return out


def mating_stats(ped: pl.DataFrame) -> dict[str, float]:
    """Return the per-generation couple, partner, MZ and allocation counts of one pedigree."""
    out: dict[str, float] = {}
    for g, gen in enumerate(generations(ped), start=1):
        out[f"n_couples_g{g}"] = len(gen.counts)
        out[f"n_mothers_multi_g{g}"] = n_multi(gen.couple_mother)
        out[f"n_fathers_multi_g{g}"] = n_multi(gen.couple_father)
        out[f"n_mz_pairs_g{g}"] = gen.n_mz
        out[f"alloc_sumsq_g{g}"] = int(np.sum(gen.counts**2))
    return out


def ztp(rng: np.random.Generator, lam: np.ndarray | float, size: int | tuple[int, ...]) -> np.ndarray:
    """Draw zero-truncated Poisson(lam): 1 + Poisson(lam - t), t the first arrival given one in [0, lam]."""
    t = -np.log1p(rng.uniform(size=size) * np.expm1(-lam))
    return 1 + rng.poisson(np.maximum(lam - t, 0.0))


Pmf = tuple[int, np.ndarray]  # (offset, log p): exp(log p[i]) is proportional to P(S = offset + i)


def trimmed(offset: int, log_p: np.ndarray) -> Pmf:
    """Scale log p to a maximum of 0 and drop the ends below double precision (exp underflows to 0) relative to it."""
    log_p = log_p - log_p.max()
    nz = np.flatnonzero(np.exp(log_p) > 0)
    return offset + int(nz[0]), log_p[nz[0] : nz[-1] + 1]


def cell_pmf(n_max: int, lam: float, shape: float | None) -> Pmf:
    """Return one cell's law on c = 1..n_max: ZTP(lam), or with shape the zero-truncated Poisson-Gamma(shape) of mean lam."""
    c = np.arange(1, n_max + 1)
    if shape is None:
        return trimmed(1, c * np.log(lam) - special.gammaln(c + 1))
    return trimmed(1, special.gammaln(c + shape) - special.gammaln(c + 1) + c * np.log(lam / (lam + shape)))


def cell_mean(lam: np.ndarray, shape: float | None) -> np.ndarray:
    """Return the mean of cell_pmf's law before its n_max cap."""
    log_p_zero = -lam if shape is None else -shape * np.log1p(lam / shape)
    return lam / -np.expm1(log_p_zero)


def split_sums(rng: np.random.Generator, s: np.ndarray, left: Pmf, right: Pmf, max_cells: int = 1 << 20) -> np.ndarray:
    """Draw each left-half sum a given the node's sum s, from P(a | s) proportional to left(a) right(s - a).

    The products are formed in logs and scaled to a row maximum of 1, so a
    row of tiny products neither underflows nor rounds its target to 0. Rows
    go in chunks of at most max_cells products; the uniforms are drawn in row
    order either way, so the chunking does not change the draws.
    """
    (o_l, lp_l), (o_r, lp_r) = left, right
    lo = np.maximum(o_l, s - (o_r + len(lp_r) - 1))
    hi = np.minimum(o_l + len(lp_l) - 1, s - o_r)
    if not np.all(lo <= hi):
        raise ValueError("a node's sum lies outside its support")
    out = np.empty_like(s)
    step = max(1, max_cells // (int((hi - lo).max()) + 1))
    for i in range(0, len(s), step):
        rows = slice(i, i + step)
        a = lo[rows, None] + np.arange(int((hi[rows] - lo[rows]).max()) + 1)
        ok = a <= hi[rows, None]
        a = np.where(ok, a, lo[rows, None])
        log_p = np.where(ok, lp_l[a - o_l] + lp_r[s[rows, None] - a - o_r], -np.inf)
        top = log_p.max(axis=1, keepdims=True)
        if not np.all(np.isfinite(top)):
            raise ValueError("a node's sum lies outside its support")
        cum = np.cumsum(np.exp(log_p - top), axis=1)  # each row's largest weight is 1, so its total is >= 1
        target = (1 - rng.uniform(size=len(a))) * cum[:, -1]  # > 0, so the pick has positive weight
        out[rows] = lo[rows] + np.sum(cum < target[:, None], axis=1)
    return out


def surjection_draws(
    rng: np.random.Generator, n: int, weights: np.ndarray, m: int, dirichlet_shape: float | None = None
) -> np.ndarray:
    """Return m exact draws of counts on K couples, each >= 1 and summing to n, shape (m, K).

    The counts are independent ZTP(mu * w_k) conditioned on summing to n: with
    equal weights Multinomial(n, 1/K) conditioned on every count >= 1, with
    weights Multinomial(n, w / sum w) so conditioned, for any mu. With
    dirichlet_shape a, the cells are instead zero-truncated Poisson-Gamma(a),
    so with equal weights the counts are Dirichlet(a)-multinomial conditioned
    on every count >= 1: an overdispersed allocation for power.py.

    The couples, sorted by weight, are split in halves down to single cells.
    Going up, each node's sum law is the direct convolution of its halves'
    (nodes with the same weights share one). Going down from the root's n,
    each node's sum is split between its halves from the exact conditional
    law. mu, set so that the mean sum is n, only centres the laws in floating
    point; the conditionals do not depend on it.
    """
    k = len(weights)
    if not k <= n:
        raise ValueError(f"{n} births cannot occupy {k} couples")
    if dirichlet_shape is not None and np.ptp(weights) > 0:
        raise ValueError("a Dirichlet-multinomial allocation takes equal weights")
    if k == n:
        return np.ones((m, k), dtype=np.int64)
    order = np.argsort(weights, kind="stable")
    w = weights[order] / weights.mean()  # the law does not depend on the weights' scale
    mu = optimize.brentq(lambda s: np.sum(cell_mean(s * w, dirichlet_shape)) - n, 1e-12, 1e6)
    laws: dict[bytes, Pmf] = {}

    def law(lo: int, hi: int) -> bytes:
        key = w[lo:hi].tobytes()
        if key not in laws:
            if hi - lo == 1:
                laws[key] = cell_pmf(n - k + 1, mu * w[lo], dirichlet_shape)
            else:
                mid = (lo + hi) // 2
                (o_a, a), (o_b, b) = laws[law(lo, mid)], laws[law(mid, hi)]
                with np.errstate(divide="ignore"):  # an underflowed sum is log 0 = -inf
                    laws[key] = trimmed(o_a + o_b, np.log(np.convolve(np.exp(a), np.exp(b))))
        return key

    law(0, k)
    draws = np.empty((m, k), dtype=np.int64)
    nodes, sums = [(0, k)], np.full((m, 1), n)
    while nodes:
        groups: dict[tuple[bytes, bytes], list[int]] = {}
        for j, (lo, hi) in enumerate(nodes):
            if hi - lo == 1:
                draws[:, order[lo]] = sums[:, j]
            else:
                mid = (lo + hi) // 2
                groups.setdefault((law(lo, mid), law(mid, hi)), []).append(j)
        next_nodes, next_sums = [], []
        for (left_law, right_law), js in groups.items():  # one call per pair of half laws on this level
            s = sums[:, js]
            left = split_sums(rng, s.ravel(), laws[left_law], laws[right_law]).reshape(s.shape)
            halves = [((lo + hi) // 2, lo, hi) for lo, hi in (nodes[j] for j in js)]
            next_nodes += [(lo, mid) for mid, lo, _ in halves] + [(mid, hi) for mid, _, hi in halves]
            next_sums += [left, s - left]
        nodes = next_nodes
        sums = np.hstack(next_sums) if next_sums else sums
    return draws


def twin_births(rng: np.random.Generator, births: np.ndarray, t: int) -> np.ndarray:
    """Count each couple's twin births: t of each row's births, uniformly without replacement (R5)."""
    m, k = births.shape
    b = int(births.sum()) // m if m else 0
    offsets = b * np.arange(m)[:, None]
    picked = np.array([rng.choice(b, t, replace=False) for _ in range(m)], dtype=np.int64).reshape(m, t) + offsets
    ends = np.cumsum(births, axis=1) + offsets
    return np.bincount(np.searchsorted(ends.ravel(), picked.ravel(), side="right"), minlength=m * k).reshape(m, k)


def trimmed_slots(rng: np.random.Generator, owners: np.ndarray, target: int) -> np.ndarray:
    """Keep a uniform random subset of target slots, in order (R2)."""
    if len(owners) == target:
        return owners
    return owners[np.sort(rng.choice(len(owners), size=target, replace=False))]


def partner_draws(rng: np.random.Generator, gen: Generation, lam: float, m: int) -> np.ndarray:
    """Return m draws of (couples, mothers with >= 2, fathers with >= 2) given the parents' sexes."""
    out = np.empty((m, 3), dtype=np.int64)
    for i in range(m):
        f = np.repeat(np.arange(gen.n_female), ztp(rng, lam, gen.n_female))
        mm = np.repeat(np.arange(gen.n_male), ztp(rng, lam, gen.n_male))
        t = min(len(f), len(mm))
        f = trimmed_slots(rng, f, t)
        mm = rng.permutation(trimmed_slots(rng, mm, t))
        occupied = rng.multinomial(gen.n_kids - gen.n_mz, np.full(t, 1.0 / t)) > 0
        out[i] = occupied.sum(), n_multi(f[occupied]), n_multi(mm[occupied])
    return out


def mz_rate(p_mz: float) -> float:
    """Return the per-child MZ-pair probability p_mz / (1 + p_mz) of R4; without the cap, E[t] / E[B] = p_mz."""
    return p_mz / (1 + p_mz)


def binomial_pit(x: int, n: int, p: float, rng: np.random.Generator, cap: int | None = None) -> tuple[float, float]:
    """Return the randomised PIT of x under min(Bin(n, p), cap), and its z under that same law."""
    at_cap = cap is not None and x >= cap
    pmf = stats.binom.sf(x - 1, n, p) if at_cap else stats.binom.pmf(x, n, p)
    u = float(randomised_pit_discrete(stats.binom.cdf(x - 1, n, p), pmf, rng))
    if cap is None or cap >= n:
        mean, sd = n * p, np.sqrt(n * p * (1 - p))
    else:
        values = np.arange(cap + 1)
        probs = np.append(stats.binom.pmf(values[:-1], n, p), stats.binom.sf(cap - 1, n, p))
        mean = float(probs @ values)
        sd = float(np.sqrt(probs @ (values - mean) ** 2))
    return u, float((x - mean) / sd) if sd > 0 else 0.0


def rank_pit(obs: float, draws: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    """Return the randomised PIT of obs among the draws, and its z."""
    n_less, n_equal = np.atleast_1d(np.sum(draws < obs)), np.atleast_1d(np.sum(draws == obs))
    u = float(randomised_pit(n_less, n_equal, len(draws), rng)[0])
    sd = draws.std(ddof=1)
    return u, float((obs - draws.mean()) / sd) if sd > 0 else 0.0


def stouffer(u: list[float]) -> float:
    """Combine independent U(0, 1) values into one exactly U(0, 1) value."""
    z = stats.norm.ppf(np.clip(u, 1e-300, 1 - 1e-16))
    return float(stats.norm.cdf(z.sum() / np.sqrt(len(z))))


def allocation_draws(
    rng: np.random.Generator, gen: Generation, m: int, mutation: str | None
) -> tuple[np.ndarray, np.ndarray]:
    """Return m null draws of (sum c^2, mean c over multi-couple mothers' couples) on the observed couples.

    Births are surjected onto the couples (R4); the allocation mutations
    perturb that step only. Each draw then places the t twin births uniformly
    among its births (R5) and adds one child per twin birth.
    """
    k = len(gen.counts)
    births = gen.n_kids - gen.n_mz
    if mutation == "alloc_overdispersed":
        draws = surjection_draws(rng, births, np.ones(k), m, dirichlet_shape=DIRICHLET_SHAPE)
    elif mutation == "alloc_by_mother":
        draws = surjection_draws(rng, births, 1.0 / gen.mother_couples, m)
    else:
        draws = surjection_draws(rng, births, np.ones(k), m)
    if gen.n_mz:
        draws = draws + twin_births(rng, draws, gen.n_mz)
    multi = gen.mother_couples >= 2
    mean_multi = draws[:, multi].mean(axis=1) if multi.any() else np.full(m, np.nan)
    return (draws**2).sum(axis=1), mean_multi


def observed_allocation(gen: Generation) -> tuple[int, float]:
    """Return the observed (sum c^2, mean c over multi-couple mothers' couples)."""
    multi = gen.mother_couples >= 2
    return int(np.sum(gen.counts**2)), float(gen.counts[multi].mean()) if multi.any() else np.nan


def mating_law(
    ped: pl.DataFrame, lam: float, p_mz: float, m: int, rng: np.random.Generator, mutation: str | None = None
) -> dict[str, tuple[float, float]]:
    """Return each Part A row's per-pedigree PIT and mean z (module docstring)."""
    lam_null = LAMBDA_ALT if mutation == "lambda_06" else lam
    q_null = mz_rate(PMZ_ALT if mutation == "pmz_03" else p_mz)
    p_male = P_MALE_ALT if mutation == "sex_052" else 0.5
    swap = mutation == "null"  # simACE's value becomes one more null draw
    per_gen: dict[str, list[tuple[float, float]]] = {}

    def add(row: str, pit: tuple[float, float]) -> None:
        per_gen.setdefault(row, []).append(pit)

    def exact(row: str, x: int, n: int, p: float, cap: int | None = None) -> None:
        if swap:
            x = min(int(rng.binomial(n, p)), n if cap is None else cap)
        add(row, binomial_pit(x, n, p, rng, cap))

    def ranked(row: str, obs: float, draws: np.ndarray) -> None:
        if np.isnan(obs):
            return  # no multi-couple mother this generation: the row has no value given the couples
        add(row, rank_pit(draws[-1], draws[:-1], rng) if swap else rank_pit(obs, draws[:m], rng))

    for gen in generations(ped):
        exact("sex_males", gen.males, gen.n_kids - gen.n_mz, p_male)
        exact("mz_pairs", gen.n_mz, gen.n_kids, q_null, cap=gen.n_kids // 2)
        sumsq, mean_multi = allocation_draws(rng, gen, m + 1, mutation)
        obs_sumsq, obs_multi = observed_allocation(gen)
        ranked("alloc_sumsq", obs_sumsq, sumsq)
        ranked("alloc_multi_mother", obs_multi, mean_multi)
        partners = partner_draws(rng, gen, lam_null, m + 1)
        observed = (len(gen.counts), n_multi(gen.couple_mother), n_multi(gen.couple_father))
        for j, row in enumerate(("n_couples", "n_mothers_multi", "n_fathers_multi")):
            ranked(row, observed[j], partners[:, j].astype(float))
    return {
        row: (stouffer([u for u, _ in pits]), float(np.mean([z for _, z in pits]))) for row, pits in per_gen.items()
    }
