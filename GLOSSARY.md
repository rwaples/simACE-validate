# simACE-validate glossary

**Oracle**:
Code that simACE does not share, which draws the quantity under test on a
pedigree simACE produced. Examples are the R liability oracle, msprime
gene-drop and SLiM. Called at a pinned version, never copied (ADR 0001).
_Avoid_: reference simulator, ground truth

**Case**:
One simACE configuration, run over R seeds and handed to one Oracle. One
entry under `cases:` in cases.yaml.

**Statistic**:
One number computed on a pedigree's values, such as a pair correlation for
one relationship class or a variance in one generation. It is computed by
the same function on simACE's values and on each Oracle draw. Calibrating
one Statistic over a Case's seeds gives one row of results.

**PIT**:
The randomised rank of simACE's Statistic among the Oracle's M draws, scaled
to (0, 1). It is uniform when simACE samples from the Oracle's law.

**Couple**:
A mother and father with at least one child together in the recorded
pedigree. Under simACE's `standard` mating, a mating that gets no children
leaves no trace, so couples are the visible subset of matings.
_Avoid_: mating, pair (for this visible subset)

**Layer**:
A grouping of Cases by what they test:
1. the liability law on a fixed pedigree
2. gene-drop through a fixed pedigree
3. pedigree structure
4. simACE's downstream stages on a gene-drop A (simACE-ts)
5. population moments under assortative mating, against xftsim (ADR 0003,
   proposed)

It is a reporting group, not a Bonferroni family (ADR 0001).

**Realised truth**:
The expectation of a Statistic given the values the Oracle supplied on the
very rows simACE used, with the rest drawn from simACE's stated law. In
simACE-ts, a liability correlation's expectation is computed from the
realised gene-drop A and the configured C and E, not from 2φ.
_Avoid_: nominal truth (for this)

**Part**:
One seed's result for a Case that runs off the laptop: its Statistics and
the provenance they were computed with. The harness scores a Case's Parts
as it scores Statistics it computed itself.

**Matched configuration**:
A simACE configuration and an outside simulator's configuration whose
population laws agree to first order, although their pedigrees differ (for
example, family-size distributions). Only population moments are compared
across a matched configuration, never pedigrees (ADR 0003).

**Reference sample**:
The committed per-seed Statistics of an outside simulator at a Matched
configuration. It is made once per pin and configuration, and simACE's
seeds are scored against it as a two-sample comparison (ADR 0003,
proposed).
_Avoid_: ground truth

**Pilot seed**:
A seed run only to measure a Case's cost before its R is fixed. Its seeds
come from a separate range and are never scored.

**Mutation**:
A deliberate error in an Oracle's law, used only by the power study.
`null` is the Mutation that replaces simACE with one more Oracle draw.

**Status**:
- PASS: both tests are at least α/K.
- FAIL: a test is below α/K.
- INFO: reported, not scored.
- GAP: not checked, with an id.
- ERROR: the run failed.
