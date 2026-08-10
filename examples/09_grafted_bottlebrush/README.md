# Example 09: Dense bottlebrush strands

Each strand has a 24-bead backbone and one 5-bead side chain on every
backbone bead. This is selected with graft spacing 0. The requested functional
fraction is 25%, so six evenly distributed side-chain ends are functional.
Backbone ends are neutral.

The example contains 750 bottlebrush strands and 1,125 four-functional
crosslinkers at `1:1` functional-group stoichiometry, for 144,000 beads in
total.

Run from the repository root:

```bash
bash examples/09_grafted_bottlebrush/run.sh
```

The generated package is written under
`examples/09_grafted_bottlebrush/grafted_bottlebrush/` and ignored by Git.
