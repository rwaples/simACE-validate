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

It is a reporting group, not a Bonferroni family (ADR 0001).

**Mutation**:
A deliberate error in an Oracle's law, used only by the power study.
`null` is the Mutation that replaces simACE with one more Oracle draw.

**Status**:
- PASS: both tests are at least α/K.
- FAIL: a test is below α/K.
- INFO: reported, not scored.
- GAP: not checked, with an id.
- ERROR: the run failed.
