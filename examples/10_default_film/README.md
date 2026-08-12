# Example 10: Default PDMS elastomer film

This is the film counterpart of `examples/01_default`. It uses the same
chemistry and composition:

- 900 linear strands with 128 beads per strand;
- two functional end groups per strand;
- 450 linear crosslinkers with 32 beads per crosslinker;
- four functional sites per crosslinker;
- `1:1` strand-end-to-crosslinker functional-group stoichiometry;
- no moderators;
- no filler;
- model-fixed 2.801 A bonds and 7.5 A initial spacing at 800 K;
- 129,600 beads in total.

The film thickness is intentionally left as the
`FILM_THICKNESS_ANGSTROM` placeholder. From the repository root, replace
`100` below with the desired nominal wall-free material thickness:

```bash
FILM_THICKNESS_ANGSTROM=100 bash examples/10_default_film/run.sh
```

For a requested thickness `H`, the generator sets the complete box length to
`Lz = H + 2*rc`, where `rc` is the 300 K wall cutoff. The generated case is
named `PDMS_elastomer_film_H<THICKNESS>` and is ignored by Git.

Additional generator options can be appended as command-line overrides:

```bash
FILM_THICKNESS_ANGSTROM=100 bash examples/10_default_film/run.sh \
    --strand-count 800
```
