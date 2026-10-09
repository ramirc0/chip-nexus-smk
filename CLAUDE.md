# CLAUDE.md

Skeleton Snakemake workflow. Rules are stubs that demonstrate the directive set.
README.md holds the layout, idioms, and style conventions. Follow it.

## Commands

```bash
pixi run lint                                # ruff check, ruff format --check, numpydoc lint
pixi run snakemake -n -p --profile profiles/local
pixi run snakemake --profile profiles/local --cores 2
```

## Two pixi workspaces

- Root `pixi.toml` is the launcher: Snakemake and its plugins, plus the `lint` env.
- `workflow/envs/pixi.toml` is the rule env. Rule dependencies go there, never in the root.
- `pyproject.toml` is ruff and numpydoc config only. pixi ignores it.
- Rules use `locked=True`. After a manifest edit, run `pixi lock` in `workflow/envs/`.
- Any manifest or lock change reruns every rule.
- `SOFTWARE_ENV` resolves `../envs` relative to the rule file. New rule files MUST sit in `workflow/rules/`.

## Gotchas

- `.smk` files are not linted. Keep helper docstrings numpydoc by hand.
- Plots use poikilos `plain`. `pk.save_figure()` writes SVG and PNG. Rules MUST declare both outputs.
- `MPLBACKEND=Agg` is set in the rule env activation. Scripts rely on it.
- Resources MUST be flat values. Never derive `mem_mb` or `runtime` from input size.
- The `* attempt` retry scaling in `profiles/slurm/config.yaml` is intentional.
- Auto mode may block `pixi lock`, `pixi install`, or workflow runs that fetch the poikilos git dependency. Ask the user to run them.

## Verifying a change

No test suite. Run the workflow end to end on throwaway inputs:

1. Copy `config/config.yaml.template` to `config/config.yaml` and `config/samples.example.tsv` to `config/samples.tsv`.
2. Create `resources/reference.txt` and `resources/sample_{a,b,c}.tsv`, each sample with a `value` column.
3. Run the workflow with the local profile.
4. Delete the inputs, `results/`, `logs/`, and `benchmarks/`. All are gitignored.
