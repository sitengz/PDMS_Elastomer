# Example 02: Bifunctional ring strands

This example uses 900 cyclic 128-bead PDMS strands. Each ring has two
regularly spaced reactive beads at sites 1 and 65. With four-functional
crosslinkers and `1:1` functional-group stoichiometry, the generator creates
450 crosslinker molecules and 129,600 total beads.

Run from the repository root:

```bash
bash examples/02_ring_bifunctional/run.sh
```

The generated package is written under
`examples/02_ring_bifunctional/ring_bifunctional/` and ignored by Git.

To sample two distinct reactive sites instead:

```bash
bash examples/02_ring_bifunctional/run.sh \
    --strand-reactive-distribution random \
    --strand-reactive-seed 12345 \
    --output data.ring_bifunctional_random
```
