# Example 03: Tetrafunctional ring strands

This example uses 900 cyclic 128-bead PDMS strands. Four distinct reactive
beads are sampled reproducibly from `strand_reactive_seed`. With
four-functional crosslinkers and `1:1` functional-group stoichiometry, the
generator creates 900 crosslinker molecules and 144,000 total beads.

Run from the repository root:

```bash
bash examples/03_ring_tetrafunctional/run.sh
```

The generated package is written under
`examples/03_ring_tetrafunctional/ring_tetrafunctional/` and ignored by Git.

For four evenly spaced sites at beads 1, 33, 65, and 97, override the
distribution:

```bash
bash examples/03_ring_tetrafunctional/run.sh \
    --strand-reactive-distribution regular \
    --output data.ring_tetrafunctional_regular
```
