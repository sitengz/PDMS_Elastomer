# Four-side-chain combs: bulk

This bulk reference determines the architecture-specific effective-strand
end-to-end distance used to construct the films.

Run from the repository root:

```bash
bash simulations/04_comb_4sidechain/bulk/run.sh
```

The generator writes `comb_4sidechain_bulk/`. Submit its LAMMPS input, verify the target
conversion and final equilibration, and run `topology_analyzer` and
`basic_network_analyzer` on the final `.npt_eq` snapshot and matching
`.info` file.

Record the mean `Ree` for effective strands in the largest connected
network. Use this same primary definition for all five architectures and pass
the value to the three film runners as `BULK_REE_ANGSTROM`.
