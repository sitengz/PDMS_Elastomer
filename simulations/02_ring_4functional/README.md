# Tetrafunctional rings

This production case uses 750 cyclic 160-bead precursors with four evenly
spaced functional sites at beads 1, 41, 81, and 121. Designed ring cycles must
be reported separately from reaction-induced network loops.

## Matched controls

- 3,000 precursor functional groups;
- 750 four-functional, 32-bead crosslinkers (3,000 crosslinker sites);
- `1:1` functional-group stoichiometry;
- 95% target conversion;
- no moderator or filler;
- 144,000 total beads;
- common density, thermal history, crosslink-site seed, and placement seed.

## Run order

Run `bulk/`, obtain its mean effective-strand `Ree`, and supply that value
as `BULK_REE_ANGSTROM` to the three film runners. Every geometry subfolder
contains its own configuration, runner, documentation, and output ignore rule.
