# Constructing Subjectivities Around Reproductive Labour in Spain (ESGE 2024 – ISSP)

Replication code for the Master's thesis *Constructing Subjectivities Around
Reproductive Labour in Spain (ESGE 2024 – ISSP Family and Changing Gender
Roles)*, Master's Degree in Computational Social Sciences, Universidad Carlos
III de Madrid, 2025–2026.

Author: Alicia Mira Guirao
Supervisor: Margarita Torre Fernández

`20260918_script_tfm_alicia_mira.R` reproduces the whole empirical analysis,
from the raw CIS microdata to every table and figure reported in the text and
in the appendices.

## Data

The analysis uses **CIS Study 3506**, *Encuesta Social General Española 2024
(II) / Familia y género (IV) (ISSP)*. The microdata are included in this
repository as `3506_num.csv`, so the script runs as soon as the repository is
cloned; no separate download is needed.

The file was downloaded from the CIS data bank on [DD Month YYYY] and is
redistributed **unmodified**, under the CIS conditions for the reuse of its
data, which authorise reproduction and distribution provided the source is
acknowledged. Source: Centro de Investigaciones Sociológicas (CIS), Study 3506.
The original is available at:

<https://www.cis.es/en/surveys/encuesta-social-general-espanola-2024-ii-esge-/-familia-y-genero-iv-issp->

The file is semicolon-separated, with a comma as the decimal mark and UTF-8
encoding, which is how the CIS distributes it. The reproducibility checks at
the start of the run (1,722 raw rows, 1,714 cases in the working sample, and
so on) confirm that the file has not changed.

## Running the analysis

```r
# from the repository root
source("20260918_script_tfm_alicia_mira.R")
```

The script must be run with the repository root as the working directory: it
reads `3506_num.csv` from the root and writes to `output/` using relative paths. Running time
is roughly fifteen to twenty minutes, dominated by the multiple imputation, the
RIF cluster bootstrap and the Oaxaca–Blinder bootstrap.

Random seeds are fixed at the top of the script, so the run is reproducible.
Block 1 checks every construction step against the benchmark figures
established during the exploratory phase and prints a pass/fail report; the
checks are non-fatal, and a failure signals that the input file or a package
default has changed.

### Requirements

R 4.3 or later and the following packages:

```r
install.packages(c("tidyverse", "survey", "psych", "lavaan", "mice",
                   "survival", "sandwich", "lmtest", "quantreg", "MASS",
                   "stargazer"))
```

**lavaan 0.7-2 or later is required.** Under `parameterization = "theta"` with
ordered data, earlier versions over-constrain the scalar measurement-invariance
model, which produces misfit that looks like threshold non-invariance but is an
artefact of the additional restriction. The script prints the degrees of
freedom of the scalar model as a check: it should report 21, not 24.

## Output

Everything is written to `output/`, which is created on the first run and is
not tracked by git:

```
output/
├── tables/
│   ├── <section>/*.html     one file per table
│   ├── _manifest.html       inventory of every table produced
│   └── _all_tables.html     every table concatenated, for the manuscript
└── figures/
    └── <section>/*.png
```

Sections are `descriptive`, `measurement`, `imputation`, `H1a`, `H1b`, `H2`,
`H3a` and `H3b`. Tables are exported as HTML rather than CSV because the final
deliverable is a Word document: the output can be opened in a browser and
pasted into the manuscript with its structure intact.

## Structure of the script

The file is organised in nine blocks, which must be run in order because each
depends on the objects created by the previous ones:

| Block | Contents |
|---|---|
| 0 | Setup: parameters, exporters, recoding helpers, reproducibility checks |
| 1 | Data preparation and the two long-format files used in H1 |
| 2 | Measurement: ordinal CFA, invariance by sex, scale construction |
| 3 | Multiple imputation of household income and composition |
| 4 | Descriptive analysis, including the exploratory step preceding each hypothesis |
| 5 | H1a: measurability of care relative to domestic work |
| 6 | H1b: self-reports versus proxy-informant reports |
| 7 | H2: family-to-work interference, material position and sex |
| 8 | H3a and H3b: naturalising discourse and perceived proportionality |
| 9 | Table manifest and session information |

## Citation

If you use this code, please cite the thesis and the data source separately.

> Mira Guirao, A. (2026). *Constructing Subjectivities Around Reproductive
> Labour in Spain (ESGE 2024 – ISSP Family and Changing Gender Roles)*.
> Master's thesis, Universidad Carlos III de Madrid.

> Centro de Investigaciones Sociológicas (2025). *Encuesta Social General
> Española 2024 (II) / Familia y género (IV) (ISSP)*. Estudio nº 3506. Madrid:
> CIS.

## Licence

Code released under the MIT Licence (see `LICENSE`). The CIS microdata
(`3506_num.csv`) are not covered by this licence: they are redistributed under,
and remain subject to, the CIS conditions of reuse.
