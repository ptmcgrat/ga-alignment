# Alignment specification

**Project:** 2a, Alignment
**Repository:** `ga-alignment` *(name assumed)*
**Document version:** 0.1, 2026-10-26 (first draft, from coordinator specification v0.7)
**Owner:** McGrath lab
**Status:** draft for discussion with the owner; no code written

This is the single specification document of project 2a. It is always replaced whole. It records what the
coordinator specification already fixes and leaves everything else open. The methods are not designed
here: the project designs them with the owner.

---

## 1. Purpose and scope

**Does:**
- align every assembly of the panel to the reference *M. zebra* assembly (GCF_041146795.1), pairwise;
- build the pangenome graph of all assemblies (unclipped).

**Does not:** call variants, call arrangements or breakpoints, or annotate genes. Projects 2b (Structure) and
3 (Variant finder) consume its output through the coordinator.

---

## 2. Inputs

| input | version | location | state |
|---|---|---|---|
| assemblies: one genome FASTA per panel assembly | panel P1 (25 assemblies, by accession version) | `<data root>/Genomes/<assembly_id>/<accession_version>/genome.fna.gz` (BGZF) | being downloaded (coordinator task); 2 of 25 were fetched on 2026-10-26 |
| repeat annotation per assembly | from project 1, panel P1 | registered output of project 1 | not yet produced |
| the panel | P1, `config/genome_panel_P1.tsv` | coordinator repository; read with `ga_common.panel` | available |
| sequence names, lengths and NCBI chromosome labels of each assembly | as fetched | `provenance.json` in the same folder | with the download |

- The registry holds no entries yet. Until it does, `ga_common.paths` gives the location of each assembly; the
  project never builds a data path by hand.
- The common library is installed pinned:
  `pip install "git+https://github.com/ptmcgrat/ga-coordinator.git@common-v0.2.0#subdirectory=common"`
  (Python 3.9 or newer).
- Whether the repeat annotation is an input to the alignment at all (soft-masking) is the project's decision
  (section 7).

**Naming caution (observed on the first downloads).** Sequence names are exactly as NCBI delivers them: RefSeq
names (`NC_...`, `NW_...`) for the 8 annotated assemblies, GenBank names (`CM...`, `CP...`) for the others. NCBI's
chromosome labels (`LG2`, ...) are the submitters' and are not guaranteed to be the same chromosome in every
assembly: for example, LG2 of `MZ_MZ4f` is 40,706,781 bases against 48,807,247 for LG2 of the reference
(`NC_135168.1`). The coordinator does not interpret this; how chromosomes are anchored to the reference is for this
project and project 2b to decide.

---

## 2a. Relatedness context (coordinator specification, section 4b)

Context, not results and not a design. Sources: Malinsky et al. 2018, Blumer et al. 2025, Kumar et al. 2025; the
species table is `config/species_context.tsv` in the coordinator repository.

| code | species | ecogroup | expected inversions vs the reference (LG2, 9, 10, 11, 13, 20) |
|---|---|---|---|
| MZ | *Maylandia zebra* (reference species) | mbuna | none |
| AC | *Astatotilapia calliptera* | generalist (stands in for the ancestor) | none |
| RC | *Rhamphochromis longiceps* | Rhamphochromis (pelagic) | LG20 left arm only |
| DT | *Diplotaxodon limnothrissa* | Diplotaxodon (pelagic) | LG9, LG11, LG20 |
| CV | *Copadichromis virginalis* | utaka | LG20; LG9 and LG11 segregating |
| MC | *Mchenga conophoros* | shallow benthic | LG20; LG9 and LG11 segregating |
| YH | *Aulonocara* sp. 'chitande type north' | deep benthic (true deep clade) | all six |
| AS | *Aulonocara stuartgranti* | deep benthic | LG20; LG9 and LG11 segregating; not LG2, 10, 13 |

- **The sequences are very close:** about 0.1–0.25% average divergence between species pairs (SNP level) within a
  radiation roughly 0.4–0.8 million years old. Expect near-identical alignments broken by large rearrangements.
- **Six large inversions**, each roughly 10–20 Mb, on LG2, 9, 10, 11, 13 and 20 (LG11 and LG20 are tandem doubles).
  The table above is a literature expectation, not a result; project 2b's arrangement table holds the results.
  Pairwise alignments and the graph must be able to represent them.
- The reference carries the ancestral state of all six. Kumar et al. note an error or polymorphism in the
  reference on LG2.
- The group tree is not strictly bifurcating (extensive gene flow), so "closely related" is a group-level summary.

---

## 3. Outputs

- Contract `alignment`: PAF with the `cs` tag, each assembly against the reference (owner: this project).
  PAF is 0-based, half-open. Exact fields are to be agreed (section 7).
- Contract `graph`: GFA, unclipped (owner: this project).
- Every output names the panel version it was built from, and is registered through the coordinator.
- Location and version-folder layout of these outputs are not yet defined (coordinator section 13: where
  registered outputs live).

---

## 4. Methods

**Not designed yet.** Fixed by the coordinator specification: pairwise alignment of every assembly to the
reference; an unclipped pangenome graph; whole genomes (decision 2026-10-26); all results in reference
coordinates. Not decided: the graph builder (Minigraph-Cactus or PGGB, coordinator open question), the
pairwise aligner and its settings, and how large inversions are handled in the pairwise step.

The plan (outputs, formats, checks) is agreed with the owner and recorded here **before any code is written**.

---

## 5. Checks (proposed; to be agreed)

- every output validates against its contract version (once defined);
- every panel assembly has a PAF; query names and lengths match the assembly's `provenance.json`;
- target names and lengths match the reference's `provenance.json`;
- outputs name the panel version; commands and tool versions are recorded.

---

## 6. Decisions

| date | decision | reason |
|---|---|---|
| 2026-10 | Pairwise alignment to the reference plus an unclipped graph | coordinator decision |
| 2026-10-26 | Whole genomes | coordinator decision |
| 2026-10-26 | Where the tools run is this project's choice | the owner wants flexibility |

---

## 7. Open questions

| question | blocking |
|---|---|
| Minigraph-Cactus or PGGB? | no |
| Is the repeat annotation (soft-masking) an input to the alignment? | no |
| Exact fields of the `alignment` and `graph` contracts; where registered outputs live | before first integration |
| How chromosomes are anchored across assemblies given mixed naming and submitter-specific LG labels (section 2) | no |
| Is the YH LG10 Y haplotig (`YH_LG10Y`, held beside the panel) aligned too? | no |
| Where the tools run, and how long jobs are kept alive on the shared server | no |

---

## 8. Status

Not started. Specification drafted; waiting for the P1 download (and the repeat annotation, if it is an input)
and for the owner's plan. Nothing is ready to integrate.

---

## 9. Environment and delivery (from coordinator specification, section 3a)

- Repository `ga-alignment`: `~/GenomeAlignments/ga-alignment` on the laptop and on mzebra (outside Dropbox).
  One person works on it at a time (the server account is shared).
- **Code is delivered as numbered zips** (`alignment_001.zip`, ...) with `MANIFEST.txt`, `INSTRUCTIONS.md` and
  `files/`. On the laptop: `git pull`, unzip the exact file, `rsync` `files/` over the repository (excluding `.git`
  and `config/local.yaml`), `git rm` anything under DELETE, check `VERSION`, commit, push. On the server: `git pull`.
- **Data root:** Dropbox `Apps/GenomeAlignments`, mirrored on mzebra at `~/Temp/GenomeAlignments`. Move version
  folders with the coordinator's `scripts/sync.py` (dry run unless `--go`; never overwrites).
- **Throttling:** every script takes `--threads` and `--max-jobs` (`GA_MAX_THREADS`, `GA_MAX_JOBS`); no limit unless asked.
- Install commands are printed in the conversation with every zip, and specifications are shared as files
  readable in the browser.
- Python 3.9 is the minimum for code that uses the common library.

---

## 10. Relatedness context

See section 2a.
