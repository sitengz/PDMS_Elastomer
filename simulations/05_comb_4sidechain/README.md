# Four-side-chain combs

This production case uses 750 grafted 160-bead precursors. The 96-bead
backbone carries four 16-bead side chains at the generator's regular graft
positions. All four side-chain ends are functional; both backbone ends are
neutral.

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
