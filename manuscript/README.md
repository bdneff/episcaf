# manuscript/

The living manuscript. Build with:

```bash
latexmk -pdf main.tex        # -> main.pdf   (latexmk -c to clean aux files)
```

Edit `sections/*.tex`; keep results in sync with `results/` and figures in `figures/`
(provenance in `figures/FIGURES.md`). `main.pdf` is committed so the current state is
viewable on GitHub without compiling.
