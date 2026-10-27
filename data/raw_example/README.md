# Raw output example

Raw simulation output of two single runs (K = 0.4, glucose input 7.0 → 15.4 → 4.2 mM over 0–60, 60–120 and 140–200 min):

| Folder | Structure | β-cell fraction | Diameter |
|---|---|---|---|
| `core/` | Core | 60 % | 114 µm |
| `mix/` | Mix | 40 % | 203 µm |

Files (two columns: time in min, value; one line every 0.1 min):

| File | Content |
|---|---|
| `gt_*.dat` | glucose, G/G0 (G0 = 7 mM) |
| `ha_*.dat` | glucagon, mean output per α cell |
| `hb_*.dat` | insulin, mean output per β cell |
| `hd_*.dat` | somatostatin, mean output per δ cell |

Plot them with:

```bash
python scripts/plot_raw_example.py      # -> results/raw_example.png
```

The full raw output of all runs is not included; it can be regenerated with the simulation code (see `docs/REPRODUCE.md`).
