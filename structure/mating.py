"""Layer 3b: simACE's `standard` mating rules, read off one recorded pedigree.

simACE's standard model (simACE methods.md, "Reproduction"; v0.2.0
`simulate.py`) makes each generation of N children from the previous one:

- R1: every parent draws ZTP(lambda) mating slots, the sexes independently.
- R2: the sex with more slots keeps a uniform random subset of them, so both
  have T = min(total female, total male) slots.
- R3: female slots in order, male slots shuffled, paired by position into T
  matings. A duplicate couple has its father swapped with a random other
  mating, for at most 5 passes.
- R4: (c_1..c_T) ~ Multinomial(N, 1/T) children per mating.
- R5: each mating with c >= 2 gets one MZ pair with probability p_mz.
- R6: sex is Bernoulli(1/2); an MZ co-twin copies the first twin's sex.

A mating with no children leaves no trace, so a generation shows its couples
(GLOSSARY): K couples with c_k >= 1 children each. Generation 0's parents are
erased (ADR 0001 §3), so only generations g >= 1 are read.

`mating_stats` gives per-generation counts for the SLiM comparison (Part B).
`mating_law` scores the same counts against each rule's own null (Part A):

| row | statistic of generation g | null given everything before g |
|---|---|---|
| sex_males | males, one co-twin of each MZ pair left out | Bin(N - t_g, 1/2), exact (R6) |
| mz_pairs | MZ pairs t_g | Bin(E_g, p_mz), E_g couples with c >= 2, exact (R5) |
| alloc_sumsq | sum of c_k^2 | uniform surjection of N onto the K couples (R4) |
| alloc_multi_mother | mean c_k over couples of mothers with >= 2 couples | the same surjection (R4) |
| n_couples | K | ZTP slots, trim, occupied matings (R1, R2, R4) |
| n_mothers_multi | mothers with >= 2 couples | the same draws (R1, R2) |
| n_fathers_multi | fathers with >= 2 couples | the same draws (R1, R2) |

Why each null is exact:

- Sex is drawn apart from mating, and a co-twin copies its twin, so the other
  children are iid Bernoulli(1/2). Each eligible mating's MZ roll is an
  independent Bernoulli(p_mz) that depends only on c >= 2.
- Multinomial cells are exchangeable. Given which K are occupied, their counts
  are N balls on K boxes conditioned on none empty, whatever T was: the
  allocation rows do not depend on lambda. iid ZTP(mu) counts conditioned on
  summing to N have that law for any mu, which only sets the acceptance rate.
- A person's couples are their occupied slots. The duplicate swap exchanges
  fathers between matings, so it keeps every person's slot count, and who
  pairs with whom does not enter these rows.

Conditioning. A generation's children are the next generation's parents, so
summed statistics would condition on earlier outcomes. Instead each row takes
one randomised PIT u_g per generation, under the law of generation g given
everything before it. Mating depends on the past only through the parents'
sexes, and each row's other conditioning (t_g, E_g, the couples) is drawn
apart from the row's own randomness. So u_1..u_{G-1} are iid U(0, 1), and
their Stouffer combination is one exactly uniform PIT per pedigree, scored
over seeds by calibration.calibrate (ADR 0001 §2).

Every exact claim assumes no duplicate couple survives the 5 swap passes;
tools/residual_duplicates.py measures how rare that is. A surviving duplicate
merges two matings into one couple, which can then carry two MZ pairs: the
invariant check below raises on that, as on any other broken MZ pair.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy import optimize, stats

from calibration import randomised_pit, randomised_pit_discrete

if TYPE_CHECKING:
    import polars as pl

# Power mutations (power.py): each changes what simACE is compared with, never
# simACE. The values sit here once; the runners and the descriptions read them.
LAMBDA_ALT = 0.6  # lambda_06
PMZ_ALT = 0.03  # pmz_03
P_MALE_ALT = 0.52  # sex_052
GAMMA_SHAPE = 10.0  # alloc_overdispersed

# std_mating (oracle mating_law): mutations of the per-rule nulls.
MUTATIONS = {
    "null": "simACE's statistic replaced by one more draw from its null (checks the scoring only)",
    "null_slim": "Part A scored on slim/standard.slim pedigrees instead of simACE's",
    "lambda_06": f"partner-count null drawn with lambda = {LAMBDA_ALT} instead of the case's",
    "pmz_03": f"MZ null with p_mz = {PMZ_ALT} instead of the case's",
    "sex_052": f"sex null with P(male) = {P_MALE_ALT}",
    "alloc_overdispersed": (
        f"allocation null with Gamma({GAMMA_SHAPE:g})-perturbed couple weights; acceptance tilts the weights kept"
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
    """Check that generation g's MZ links are reciprocal and same-sex; return each pair's lower id."""
    first = kids[twin[kids] > kids]
    second = twin[first]
    if 2 * len(first) != np.sum(twin[kids] >= 0) or np.any(gen[second] != g) or np.any(twin[second] != first):
        raise ValueError(f"generation {g}: an MZ link is not reciprocal within the generation")
    if np.any(sex[first] != sex[second]):
        raise ValueError(f"generation {g}: an MZ pair differs in sex")
    return first


def generations(ped: pl.DataFrame) -> list[Generation]:
    """Read every generation g >= 1, checking parents and the MZ invariants; rows must be in id order.

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
        prev = np.flatnonzero(gen == g - 1)
        for label, parents, want in (("mother", mother[kids], 0), ("father", father[kids], 1)):
            if np.any(parents < 0) or np.any(gen[parents] != g - 1) or np.any(sex[parents] != want):
                raise ValueError(f"generation {g}: a {label} is not a generation-{g - 1} individual of that sex")
        keys, couple, counts = np.unique(
            np.stack([mother[kids], father[kids]]), axis=1, return_inverse=True, return_counts=True
        )
        first = check_twins(g, kids, gen, sex, twin)
        second = twin[first]
        if np.any(mother[first] != mother[second]) or np.any(father[first] != father[second]):
            raise ValueError(f"generation {g}: an MZ pair has different parents")
        pair_couple = couple[np.searchsorted(kids, first)]
        if len(np.unique(pair_couple)) < len(pair_couple):
            raise ValueError(f"generation {g}: a couple has two MZ pairs")
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


def surjection_draws(
    rng: np.random.Generator, n: int, weights: np.ndarray, m: int, gamma_shape: float | None = None
) -> np.ndarray:
    """Return m draws of counts on K couples, each >= 1 and summing to n, shape (m, K).

    The counts are independent ZTP(mu * w_k) conditioned on summing to n. With
    equal weights that is Multinomial(n, 1/K) conditioned on every count >= 1,
    for any mu; mu, set so that the mean sum is n, only sets the acceptance
    rate. With gamma_shape, each proposal draws its own weights
    w_k ~ Gamma(gamma_shape) rescaled to mean 1: an overdispersed allocation on
    the same occupied couples, for power.py.

    Rejection sampling: K - 1 counts are drawn, the last is n minus their sum,
    and the proposal is kept with probability p(last) / max p, where p is the
    last cell's ZTP pmf (the largest-weight cell goes last).
    """
    k = len(weights)
    if k == n:
        return np.ones((m, k), dtype=np.int64)
    order = np.argsort(weights, kind="stable")
    w = weights[order]
    mu = optimize.brentq(lambda s: np.sum(s * w / -np.expm1(-s * w)) - n, 1e-12, 1e6)
    batch = max(64, (1 << 22) // k)
    got: list[np.ndarray] = []
    n_got = 0
    while n_got < m:
        lam = np.broadcast_to(mu * w, (batch, k))
        if gamma_shape is not None:
            g = rng.gamma(gamma_shape, size=(batch, k))
            lam = mu * g * (k / g.sum(axis=1, keepdims=True))
        x = ztp(rng, lam[:, :-1], (batch, k - 1))
        last = n - x.sum(axis=1)
        lam_last = lam[:, -1]
        mode = np.maximum(1, np.floor(lam_last))
        log_max = np.maximum(stats.poisson.logpmf(mode, lam_last), stats.poisson.logpmf(mode + 1, lam_last))
        log_p = stats.poisson.logpmf(np.maximum(last, 1), lam_last)
        keep = (last >= 1) & (np.log(rng.uniform(size=batch)) < log_p - log_max)
        got.append(np.column_stack([x[keep], last[keep]]))
        n_got += int(keep.sum())
    draws = np.empty((m, k), dtype=np.int64)
    draws[:, order] = np.concatenate(got)[:m]
    return draws


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
        occupied = rng.multinomial(gen.n_kids, np.full(t, 1.0 / t)) > 0
        out[i] = occupied.sum(), n_multi(f[occupied]), n_multi(mm[occupied])
    return out


def binomial_pit(x: int, n: int, p: float, rng: np.random.Generator) -> tuple[float, float]:
    """Return the randomised PIT of x under Bin(n, p), and its z."""
    u = float(randomised_pit_discrete(stats.binom.cdf(x - 1, n, p), stats.binom.pmf(x, n, p), rng))
    sd = np.sqrt(n * p * (1 - p))
    return u, float((x - n * p) / sd) if sd > 0 else 0.0


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
    """Return m null draws of (sum c^2, mean c over multi-couple mothers' couples) on the observed couples."""
    k = len(gen.counts)
    if mutation == "alloc_overdispersed":
        draws = surjection_draws(rng, gen.n_kids, np.ones(k), m, gamma_shape=GAMMA_SHAPE)
    elif mutation == "alloc_by_mother":
        draws = surjection_draws(rng, gen.n_kids, 1.0 / gen.mother_couples, m)
    else:
        draws = surjection_draws(rng, gen.n_kids, np.ones(k), m)
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
    p_null = PMZ_ALT if mutation == "pmz_03" else p_mz
    p_male = P_MALE_ALT if mutation == "sex_052" else 0.5
    swap = mutation == "null"  # simACE's value becomes one more null draw
    per_gen: dict[str, list[tuple[float, float]]] = {}

    def add(row: str, pit: tuple[float, float]) -> None:
        per_gen.setdefault(row, []).append(pit)

    def exact(row: str, x: int, n: int, p: float) -> None:
        add(row, binomial_pit(int(rng.binomial(n, p)) if swap else x, n, p, rng))

    def ranked(row: str, obs: float, draws: np.ndarray) -> None:
        if np.isnan(obs):
            return  # no multi-couple mother this generation: the row has no value given the couples
        add(row, rank_pit(draws[-1], draws[:-1], rng) if swap else rank_pit(obs, draws[:m], rng))

    for gen in generations(ped):
        exact("sex_males", gen.males, gen.n_kids - gen.n_mz, p_male)
        exact("mz_pairs", gen.n_mz, int(np.sum(gen.counts >= 2)), p_null)
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
