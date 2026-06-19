# PhyloMLST

PhyloMLST is a command-line tool for clustering MLST profiles into clonal complexes, assigning clonal complexes to international clones (ICs) (or equivalent named clonal complex) using founder sequence types (STs), and optionally correcting those assignments using a phylogenetic tree.

---

## Installation

No installation is required. PhyloMLST is a single Python script with the following dependencies:

```
pandas    >= 3.0.2
networkx  >= 3.6.1
biopython >= 1.87
```

Install dependencies with pip:

```bash
pip install pandas networkx biopython
```

Then place `phylomlst.py` in your working directory and run it directly:

```bash
python phylomlst.py --help
```

---

## Input

PhyloMLST is designed to take the TSV output of [FastMLST](https://github.com/EnzoAndree/FastMLST) as its primary input. When running FastMLST, use the `-to` and `-s` flags to ensure tab-separated output:

```bash
fastmlst --scheme 'abaumannii#2' -s '\t' -to mlst.tsv A_baumannii.fna
```

The MLST input file must contain at minimum a `Genome` column and an `ST` column. All columns other than `Genome`, `ST`, `Scheme`, `clonal_complex`, and `species` are treated as marker genes.

### Founder ST file (`--founder`)

A two-column, headerless TSV mapping founder STs to their IC (or equivalent clonal complex name):

```
ST2    IC1
ST1    IC2
```

### ST-to-IC baseline file (`--st_to_ic`, optional)

A two-column, headerless TSV mapping IC names to a space-separated list of known STs:

```
IC1    1 19
IC2    2 45 187
```

---

## Usage

```bash
python phylomlst.py -m mlst.tsv -f founder_sts.tsv [options]
```

### Arguments

```
required arguments:
  -m, --mlst            TSV file containing MLST data (FastMLST output)
  -f, --founder         TSV file mapping founder STs to ICs

optional arguments:
  -t, --tree            Phylogenetic tree in Newick format for IC correction
                        (default: None)
  -s, --st_to_ic        TSV file mapping ICs to a baseline set of known STs
                        (default: None)
  -o, --output          Output file path
                        (default: out.tsv)
  -r, --report          Write a TSV summary of IC counts to this path
                        (default: None)
  -l, --locus_threshold
                        Number of locus differences allowed when grouping
                        isolates into clonal complexes; accepted values are 1-6
                        (default: 1)
  -n, --not_ic_threshold
                        Number of Not IC isolates permitted within a subtree
                        before it is considered unclean during phylogenetic
                        correction; must be a non-negative integer
                        (default: 0)
  -x, --mask_non_ST
                        Mask CC and IC assignments for isolates with
                        non-numeric STs (e.g. novel allele calls)
```

---

## Examples

**Minimal run** — assign clonal complexes and ICs from MLST data:

```bash
python phylomlst.py \
    -m mlst.tsv \
    -f founder_sts.tsv \
    -o results.tsv
```

**With phylogenetic correction** — correct IC assignments using a Newick tree:

```bash
python phylomlst.py \
    -m mlst.tsv \
    -f founder_sts.tsv \
    -t tree.nwk \
    -o results.tsv
```

**Full run** — all options enabled:

```bash
python phylomlst.py \
    -m mlst.tsv \
    -f founder_sts.tsv \
    -t tree.nwk \
    -s st_to_ic.tsv \
    -o results.tsv \
    -r ic_report.tsv \
    -l 1 \
    -n 2 \
    -x
```

---

## Output

The primary output (`--output`) is a tab-separated file with the original MLST data plus the following appended columns:

| Column | Description |
|---|---|
| `CC` | Clonal complex ID (integer, ranked by size) |
| `Assigned IC` | IC assignment based on founder ST and CC membership |
| `Corrected IC` | IC assignment after phylogenetic correction (only present if `--tree` is provided) |

---

## Terminology

PhyloMLST was designed with *Acinetobacter baumannii* in mind, where globally disseminated lineages are formally named **international clones (ICs)**. However, the tool is applicable to any bacterial species that uses MLST-based typing. The equivalent concept exists across species under different names — **clonal groups (CGs)** in *Klebsiella pneumoniae*, **high-risk clones** in *Pseudomonas aeruginosa*, or simply named STs/clonal complexes in *E. coli*, *S. aureus*, and others. Regardless of the species-specific terminology, PhyloMLST uses the term "IC" internally to refer to any named clonal complex defined by a founder ST.
