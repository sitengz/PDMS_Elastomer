# Example 01: Default PDMS elastomer

This example reproduces the generator's current default model:

- 900 linear strands with 128 beads per strand;
- two functional end groups per strand;
- 450 linear crosslinkers with 32 beads per crosslinker;
- four functional sites per crosslinker;
- `1:1` strand-end-to-crosslinker functional-group stoichiometry;
- no moderators;
- no filler;
- 129,600 beads in total.

From the repository root, run:

```bash
bash examples/01_default/run.sh
```

The script builds the generator and creates
`examples/01_default/PDMS_elastomer/` containing:

- `data.PDMS_elastomer`;
- `in.PDMS_elastomer`;
- `submit.PDMS_elastomer.sh`;
- `PDMS_elastomer.info`.

The generated directory is ignored by Git. To vary one setting without
editing `model.conf`, append a command-line override:

```bash
bash examples/01_default/run.sh --strand-count 800 --output data.strands_800
```
