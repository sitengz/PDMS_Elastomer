# Short-linear, high-crosslinker-loading reference

This production case uses 3,000 bifunctional linear chains with 40 beads per
chain. Stoichiometry derives 1,500 four-functional, 32-bead crosslinkers.
The system therefore contains 3,000 linear contours, matching the designed
ring-arc and star-arm counts while also matching their 40-bead contour length.

## Composition

- 6,000 precursor functional groups;
- 1,500 four-functional crosslinkers (6,000 crosslinker sites);
- `1:1` functional-group stoichiometry;
- 95% target conversion, corresponding to 5,700 target reaction bonds;
- no moderator or filler;
- 168,000 total beads;
- common density, thermal history, crosslink-site seed, and placement seed.

This case has twice the crosslinker loading relative to the precursor mass,
but its crosslinker number density is 1.714 times that of
`01_linear_reference` because the additional crosslinkers increase the final
material volume.

## Run order

1. Run `bulk/`.
2. Analyze the equilibrated bulk snapshot and record the selected mean
   effective-strand `Ree` in angstrom.
3. Supply that value as `BULK_REE_ANGSTROM` to `film_2Ree/`,
   `film_4Ree/`, and `film_8Ree/`.

Every geometry subfolder contains its own `model.conf`, `run.sh`,
`README.md`, and output ignore rule.
