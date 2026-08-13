# Bifunctional linear reference

This production case uses 1,500 linear 80-bead precursors with two functional
chain ends. It is the classical bifunctional-chain reference. Its parent-chain
end-to-end distance and its reduced effective-strand end-to-end distance are
identical when both ends react.

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
