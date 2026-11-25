# Alignment specification

**Project:** 2a, Alignment
**Repository:** `ga-alignment` (`github.com/ptmcgrat/ga-alignment`)
**Document version:** 0.11, 2026-11-25 (0.1: first draft from coordinator specification v0.7; 0.2: P1 downloaded and chromosome labels, from v0.8; 0.3: the first conversation starts at the briefing; working with a student; using the common library, from v0.10; 0.4, 2026-11-18: questions 1 and 2 decided by the student, to whom the owner left them; gate 0 closed; the environment; 0.5, 2026-11-19: Cactus 4.0.0; gate 1 closed; the pairwise run details; the input check; 0.6, 2026-11-20: the input check's results; the B chromosome; Kumar's assemblies by name; gate 2's answers; the pairwise order, threads and version folder; 0.7, 2026-11-24: all 28 pairwise runs, `asm5` confirmed; gate 2 closed; Kumar et al.'s Version of Record (assemblies identified; Supplementary files 2 and 4 as recorded inputs); the scoring rules for the pairwise checks; 0.8, 2026-11-24: column 11 and ambiguous bases; the published tables moved to `published/`; 0.9, 2026-11-24: the first check results; inversion calls from the flanks next to the units; breakpoints scored for YH only; column 10 and N facing N; 0.10, 2026-11-24: the pairwise checks pass; gate 3 closed; the graph settings; the graph runner; 0.11, 2026-11-25: the chromosome graph runs; the graph checks)
**Owner:** McGrath lab
**Status:** questions 1 and 2 decided (sections 4 to 6); environment ready; inputs checked; pairwise alignments done and checked; gates 0 to 3 closed; the graph's chromosome runs done; their checks (zip `011`) run next

This is the single specification document of project 2a. It is always replaced whole. It records what the
coordinator specification already fixes, and the methods the project has decided. Questions 1 and 2 were left to the
student by the owner (2026-11-17); what is still open is in section 7.

---

## 0. Where the first conversation starts, and who is working

**The person working on this project is a student, supervised by the owner (the professor).** The assistant acts as a mentor
and checker: it explains, asks questions, and checks the student's understanding against the documentation and against
evidence the student pastes. The student proposes; **the owner decides.** Where a choice is the owner's, the assistant writes it
into section 7 as an open question and does not settle it.

**Read at the start of every session:** this specification and `docs/alignment_briefing.md`.

**Before the first conversation, the student uploads the reference PDFs** listed in the briefing, Part D (items 1–8 are required,
9 and 10 may follow). The assistant's first action is to list the documents it can see in the project, with the title of each, and
to ask for anything missing or unreadable. It does not go on until items 1–8 are there.

**The first conversation follows this order:**
1. The student reads Kumar et al. 2025 (the six inversions, the reference, and the alignments in their Figures S2–S7) and Quah et al.
   2025. The assistant goes through Parts A and B of the briefing with the student.
2. **Question 1:** how should the pairwise alignments and the graph be built, and how do we compare the options on a subset (briefing, Part B)?
3. **Question 2:** how will we know the alignments and the graph are right (briefing, Part C)?
4. The student writes their proposals for questions 1 and 2 into sections 6 and 7 for the owner to review.

**No installation and no alignment or graph run** until the gate 0 questions (section 4a) are answered and the owner has agreed the
plan for questions 1 and 2, unless the owner says otherwise.

Anything that needs another project goes to the coordinator (blocking items and integration requests only).

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
| assemblies: one genome FASTA per panel assembly | panel P1 (25 assemblies, by accession version) | `<data root>/Genomes/<assembly_id>/<accession_version>/genome.fna.gz` (BGZF) | **downloaded 2026-10-26**: all 25 assemblies (with GFF3 for the 8 annotated), in Dropbox and on mzebra |
| repeat annotation per assembly | from project 1, panel P1 | registered output of project 1 | not yet produced |
| the panel | P1, `config/genome_panel_P1.tsv` | coordinator repository; read with `ga_common.panel` | available |
| sequence names, lengths and NCBI chromosome labels of each assembly | as fetched | `provenance.json` in the same folder | available |
| accessions of LG2, LG9, LG10, LG11, LG13, LG20 per assembly | P1, `config/panel_chromosomes_P1.tsv` | coordinator repository (generated from the downloads) | available |

- The registry holds no entries yet. Until it does, `ga_common.paths` gives the location of each assembly; the
  project never builds a data path by hand.
- The common library is `ga_common`, in the `common/` folder of the coordinator repository (`ga-coordinator`). Install and
  use it as `docs/USING_COMMON.md` in that repository describes: on mzebra,
  `pip install -e ~/GenomeAlignments/ga-coordinator/common`, or pinned to the tag `common-v0.2.1`. It finds the coordinator's
  configuration through `GA_CONFIG_DIR` or the standard location `~/GenomeAlignments/ga-coordinator/config` (Python 3.9 or newer).
- Whether the repeat annotation is an input to the alignment at all (soft-masking) is the project's decision
  (section 7).

**Naming and label cautions (observed on the P1 download; coordinator specification, section 4).**
- Sequence names are exactly as NCBI delivers them: RefSeq names (`NC_...`) for the 8 annotated assemblies,
  GenBank names (`CM...`, `OZ...`) for the others.
- NCBI's chromosome labels are the submitters'. 17 assemblies use `LG1`–`LG20`, `LG22`, `LG23`; 8 use plain
  numbers with the same set (`AC_SAMEA123396440`, `AC_SAMEA123396441`, `AC_fAstCal68`, `AS_SAMEA123396442`,
  `AS_Sanger`, `DT_SAMEA123396444`, `MZ_SAMEA123396453`, `MZ_SAMEA123396454`). That a number is the same
  chromosome as the reference's LG is a **working assumption**, not a result. Accession order is not a guide to
  labels (in `AS_Sanger`, label 9 is `NC_138960.1`, after label 23's `NC_138958.1`).
- **Length observations** (not interpreted by the coordinator; assembly completeness may explain some):
  the reference's LG2 is 48.8 Mb against 38.5–42.3 Mb in the 24 other assemblies; LG3 varies 1.8-fold across
  assemblies (58.7–105.3 Mb; reference 64.4); LG20 is only 17.8 Mb (`MC_MC3m`) and 20.3 Mb (`CV_CV2f`) against
  35–42 Mb in most, and 44.5–46.1 Mb in `DT_004f` and `CV_CV4f`; LG10 of `YH_YH2f` is 28.4 Mb against 34–42 elsewhere;
  `MZ_UMD2025` is short on LG3, LG9, LG11 and LG13. Kumar et al. 2025 describe the reference's left arm of chromosome 2
  and the highly repetitive chromosome 3 as the places where alignments of *M. zebra* genomes to it differ most. How
  chromosomes are anchored across assemblies is for this project and project 2b to decide.
- `YH_LG10Y` (the Y haplotig, held beside the panel) is at `Genomes/YH_LG10Y/h2tg000014l/`, with `genome.fna.gz` and
  `provenance.json` but no `sequence_report.jsonl` (observed 2026-11-18).
- **Input check** (`scripts/check_inputs.py`, 2026-11-19): all 26 genome files match their recorded sha256, and every sequence
  name and length matches `provenance.json`. The 150 labels in `panel_chromosomes_P1.tsv` agree with NCBI's, including the five
  LG10s the panel marks as inferred. The reference has 22 linkage groups (933.9 Mb) and 41 other sequences (28.2 Mb).
- **A B chromosome:** `MZ_UMD2025` carries, besides LG1–LG20, LG22 and LG23, a 26.7 Mb sequence labelled `B` (`CM105978.1`),
  very likely a B chromosome. Its alignments may look like rearrangements.
- **Soft-masking and gaps:** all 25 NCBI files are soft-masked (31–34% lower case); `YH_LG10Y` is not. The Georgia Tech
  assemblies carry 3–37 Mb of N (most in `CV_CV4f`, 37 Mb; `CV_CV2f`, 34 Mb; `DT_004f`, 30 Mb); the Sanger and Cambridge
  ones almost none.
- **Kumar et al.'s assemblies** (Version of Record: Supplementary file 3, Data availability): MZ_GT3.2 is `MZ_NMetZ`
  (959,221,546 bp against our 959,220,039), MZ_GT3.3 is `MZ_MZ4f`, YH_GT1.1 is `YH_YH1m` (male), YH_GT1.2 is `YH_YH2f` (female)
  and YH_GT1.3 is `YH_YH7f`: matched by total size, each within 0.4 Mb while any other pairing differs by more than 10 Mb, and by
  sex in Supplementary file 2. The Version of Record's figure captions swap the YH sexes ("YH_GT1.1f, YH_GT1.2m"); the labels
  inside its plots, version 2, the NCBI names and Supplementary file 2 all agree with the matching here.
- **The reference is GT3a, in Table S4's coordinates:** Kumar et al.'s Data availability names GCF_041146795.1 as GT3a's RefSeq
  annotation, and their Supplementary file 1 gives GT3a chromosome lengths equal to ours (LG1 45,605,991; LG3 64,350,661;
  LG9 37,905,347).
- **`YH_YH1m`'s LG9** (pairwise, 2026-11-20): reference 0–6.2 Mb lies at the start of YH's LG9 on "+", while the rest of the
  chromosome runs reversed; it belongs next to reference 6.8–9.7 Mb, which lies at the other end (YH 34.0–37.2 Mb). Kumar et
  al.'s YH_GT1.1 had it as a separate scaffold; most likely a scaffolding join with one piece reversed (flagged to project 2b).

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

- Contract `alignment`: PAF with the `cs` tag, each assembly against the reference (owner: this project): one PAF per
  assembly, written by minimap2 with `-c --cs` (section 4). PAF is 0-based, half-open. Exact fields are to be agreed (section 7).
- Contract `graph`: GFA, unclipped (owner: this project): Minigraph-Cactus's full graph, `--gfa full` (section 4).
- Every output names the panel version it was built from, and is registered through the coordinator.
- Location and version-folder layout of these outputs are not yet defined (coordinator section 13: where
  registered outputs live). Until the coordinator adds data-root folders for them, results stay in `~/scratch/ga-alignment`.

---

## 4. Methods

Decided by the student, to whom the owner left questions 1 and 2 (2026-11-17); the reasons are in section 6.
Fixed by the coordinator specification: pairwise alignment of every assembly to the reference; an unclipped pangenome
graph; whole genomes (decision 2026-10-26); all results in reference coordinates.

### Pairwise alignment (contract `alignment`)

- **Tool:** minimap2 2.31 (`2.31-r1302`, bioconda), one version for every PAF; the version the binary reports is recorded with
  each output.
- **Target and queries:** the reference `MZ_NMK2024a` against each of the other 24 P1 assemblies; whole genomes, sequence
  names as delivered, inputs as in section 2.
- **Settings:** base-level alignment with the `cs` tag (`-c --cs`), preset `asm5` for all 24.
- **Sensitivity check:** one `asm10` run on `YH_YH1m`, compared with its `asm5` alignment on coverage and identity per
  reference chromosome, the six inversions against Table S4 and Figures S2–S7 of Kumar et al., fragmentation near the
  breakpoints, LG3 coverage, and wall time and peak memory. If the two agree, `asm5` stands and the comparison is reported.
- **Reverse direction** (the reciprocity check): the reference aligned to `YH_YH1m` and to `AC_SAMEA123396441`.
- **Per assembly:** one PAF, plus a log of command, version, input checksums, wall time and peak memory.
- **Records kept:** minimap2's default secondary alignments (up to 5 per primary, per its README); every record is kept in the
  PAF, and statistics use primary alignments only. Whether the registered PAF is filtered is open (section 7).
- **Reference index:** built from the FASTA in each run, minimap2's default. No saved index (`-d`): a saved index fixes the
  indexing parameters (`-k`, `-w`, `-H`, `-I`) of the preset it was built with (minimap2 README), and building it takes minutes.
- **Where the PAFs go:** `~/scratch/ga-alignment/pairwise/PW_v1/<assembly_id>.<preset>.paf` (reverse runs in `reverse/`), with a
  log per run and a `manifest.tsv`, until the coordinator adds a data-root folder for them (section 3).
- **Order** (decided 2026-11-20): `YH_YH1m` with `asm5` and `asm10` side by side; then a review (gate 2, question 10, on its LG9
  records; the two presets compared) that settles the preset; then the other 23 with that preset, and the two reverse runs.
- **Runs, done:** 28 in `PW_v1`, all with exit status 0: the 24 with `asm5`, `YH_YH1m` with `asm10`, the reverse runs for
  `YH_YH1m` and `AC_SAMEA123396441`, and `YH_LG10Y`. Each took 12–24 min and 27–37 GB; a run keeps about four threads busy,
  so the bulk ran 6 at a time with 8 threads each.
- **`asm5` confirmed** for all 24 after the comparison on `YH_YH1m`: `asm10` covered 1.0 point more of the reference (2.8 on
  LG3) with a third fewer records, but its gap-compressed divergence was 0.54% against 0.47% (LG3: 0.97% against 0.72%), so
  what it adds is about 7% divergent: repeat copies rather than the same sequence in both genomes. The `asm10` PAF stays in
  `PW_v1` as the record of the check.
- **Threads:** 16 per run (`-t 16`): two runs at once for `YH_YH1m`, then six (96 threads), under `nice -n 10`, in `screen`.
  One run with 96 threads would leave most of them idle: minimap2 handles each query sequence as a unit, within batches of
  200 Mb (README, steps 2–7).
- **Looks correct**, before the other 23: minimap2 exits 0 and the log records time and peak memory; every record carries `cs`;
  all 22 linkage groups are well covered by primary blocks, except LG2's left arm and parts of LG3; the six inversions show as
  inverted blocks near their Table S4 intervals, LG11 and LG20 as doubles. (The identity criterion first written here, 99.7–99.9%,
  was wrong for the measure computed: column 10 / column 11 counts gap bases. The gap-compressed divergence, the `de` tag, is the
  measure for it: 0.47% on `YH_YH1m`.)
- **Scripts:** `scripts/align_pairwise.py` (the runs, their logs and the manifest), `scripts/paf_overview.py` (per-chromosome
  tables and block lists of one PAF), `scripts/check_pairwise.py` (the exact and accounting checks of section 5) and
  `scripts/summarize_pairwise.py` (the known-feature and measured checks), with `scripts/pafkit.py` shared by the last two.
- **Recorded inputs** (`published/`; not `data/`, which the repository's `.gitignore` keeps out because data lives
  under the data root): `kumar2025_suppfile4_breakpoints.tsv` (Kumar et al.'s Supplementary file 4, the eight
  breakpoint units) and `expected_inversions.tsv` (each assembly's expected state per unit, from their Supplementary file 2),
  each with its source file and that file's sha256.
- **`YH_LG10Y`** (the Y haplotig) is aligned too, after the 24 panel assemblies.

### Graph (contract `graph`)

- **Tool:** `cactus-pangenome` from Cactus v4.0.0 (released 2026-11-09), run under Apptainer from the release image
  `quay.io/comparative-genomics-toolkit/cactus:v4.0.0` (`~/containers/cactus_v4.0.0.sif`, sha256
  `f10a1389628d4f086e60269bf652802c00bb0aac912952742141c5a6a44daae5`). Read from the image itself (2026-11-19): its default
  configuration sets `zipAlleles="1"` and `GAFOverlapFilterExact="1"`, and every planned option is in its help texts. Its
  binaries are built for x86-64-v3, which mzebra's CPUs support. Sample names follow PanSN in 4.0.0 (`SAMPLE#HAPLOTYPE`); the
  assembly ids, none of which ends in `.N`, are used as they are.
- **Reference:** `--reference MZ_NMK2024a`; `--refContigs` lists its 22 linkage-group sequences, from its `provenance.json`;
  `--otherContig` keeps the rest; sample names are the assembly ids.
- **Outputs:** the full graph (`--gfa full`) is the `graph` deliverable. On the subset the clipped graph is also written, to
  measure what clipping would remove. Cactus's stats directory is kept with every run. `--validate` (new in 4.0.0) checks the
  output base for base against the inputs, alongside our own path check (section 5).
- **Settings** (step 4, agreed 2026-11-24): `--reference MZ_NMK2024a`; `--refContigs` with its 22 linkage groups and
  `--otherContig chrOther`, which lumps its 41 other sequences into one; `--gfa full clip` on the subset and `--gfa full` on
  the full graph (`--gfa` alone gives the clip graph); `--validate`; Cactus's default chromosome-assignment thresholds (no
  `--permissiveContigFilter`); the configuration's deletion filters as they are (`delFilter` 10 Mb, `delFilterThreshold`
  0.01, `delFilterMaxRemove` 1 Mb, `delFilterMinSupport` 2, `delFilterQuerySizeThreshold` 2; in the clip graph
  `delEdgeFilter` 10 Mb). The clip graph is written on the subset only, to measure what clipping removes and to check our
  accounting against Cactus's `clipped-by-reference-contig.tsv`.
- **Caps:** `taskset` confines a run to a CPU range; the per-step options (`--mgCores`, `--mapCores`, at most 6, `--consCores`,
  `--indexCores`, and `--maxCores`) size each job; `--maxMemory` is the ceiling `--consRetainPages auto` sizes against. In
  4.0.0, `--maxCores`, `--maxMemory` and `--maxDisk` limit what one job may request, not the run as a whole (its help). Runs
  go under `nice -n 10` in `screen`.
- **Runs:** in `GR_v1`, one folder per run holding its job store, working directory and outputs; after an interruption,
  `--restart` resumes from the job store. Chromosome runs first (LG8, LG9, LG20), from per-chromosome FASTAs chosen by the
  pairwise assignment, then one whole-genome run of the subset; `scripts/build_graph.py`.
- **The full graph** (step 6) is built from scratch, not by adding genomes to the subset graph (`--inGFA`): only the minigraph
  construction would be reused, since the alignment reruns for every genome; reusing the mappings needs `zipAlleles="0"`;
  and construction is iterative in the input genomes, so the full graph would depend on the subset's order. Whether
  `YH_LG10Y` joins the full graph is decided at step 6.
- **Settings:** documented defaults. Two candidates for change, each decided from the pairwise results and recorded in
  section 6 with its reason: `--delFilter` (the reference's LG2) and `--permissiveContigFilter` (LG3). Every option is checked
  against the installed version's own help before it is used.
- **Throttling:** `--maxCores`, `--mgCores`, `--mapCores`, `--consCores`, `--indexCores` and `--maxMemory` set explicitly.
- **Stages:** (1) the subset on LG2, LG3, LG9, LG10, LG11, LG13 and LG20 plus one control chromosome, with each assembly's
  sequences chosen from the pairwise alignments, not from labels; (2) one whole-genome run of the subset, to test Cactus's own
  chromosome assignment and measure cost at a second size; (3) all 25 assemblies in one run, built from scratch rather than
  extended from the subset graph.

### The subset

| role | id | why |
|---|---|---|
| reference | `MZ_NMK2024a` | required |
| no inversions; label test | `AC_SAMEA123396441` | one contig per chromosome; plain-number labels |
| same species | `MZ_NMetZ` and `MZ_MZ4f` | Kumar's two *M. zebra* controls (MZ_GT3.2f and GT3.3f); Figures S2–S7 apply directly |
| all six inversions | `YH_YH1m` | Kumar's YH_GT1.1m (assembly name A_spYH_male_GT1.1); LG10 noted as an X/Y mosaic |
| all six; sequence split across chromosomes | `YH_YH2f` | very likely Kumar's YH_GT1.2f; LG10 X/X; in his Figure S4 its LG10 is split across scaffolds |
| LG20a only | `RC_006pf` | *R. longiceps*, RefSeq names |
| *Diplotaxodon* | `DT_SAMEA123396444` | inverted on 9, 11 and 20, not on 2, 10 or 13 |

The seven are reviewed first against Kumar's figures, and they are the graph's subset.

### What the graph needs from the pairwise alignments

| from the PAFs | what it decides for the graph |
|---|---|
| the reference chromosome each query sequence aligns to, unplaced scaffolds included | the sequences in each chromosome run (not the labels); a prediction of Cactus's own assignment, checked against its split log |
| sequences aligning substantially to two chromosomes | what Cactus will drop at the split; accept and report it, or raise it with 2b |
| size and position of the reference-only stretch on LG2 | whether `--delFilter`'s installed default would drop mappings there, and the value to test |
| how much of each assembly's LG3 aligns, and to what | whether any LG3 risks being called ambiguous (`--permissiveContigFilter`), and what `--otherContig` must keep |
| reference sequences outside the 22 linkage groups that receive alignments | what `--otherContig` keeps |
| each assembly's inverted blocks and large insertions | the yardsticks the graph must reproduce |
| coverage and identity per chromosome | the choice of control chromosome |

The subset's graph runs need only the subset's alignments. The full graph needs all 24 reviewed first: sequences split across
chromosomes, or LG3s at risk, can sit in assemblies outside the subset (`MC_MC3m` and `CV_CV2f` are).

### Order of work

Each step's code comes as the next numbered zip; each step ends with results sent back and reviewed. Gates (section 4a) ask
their questions before their commands.

1. **Environment** (zip `004`): conda env, Cactus image, `scripts/check_env.py`. Done 2026-11-19, with Cactus 4.0.0.
2. **Inputs** (gate 1, zip `005`; done 2026-11-19, every check passed): locate the 25 genomes and `YH_LG10Y` through `ga_common`; check names and lengths against `provenance.json`; one
   sequence table per assembly (name, length, NCBI chromosome label, role); the reference's chromosomes; Kumar's assemblies
   among the panel ids. Decides the final subset.
3. **Pairwise** (gate 2, zip `006`; done): `YH_YH1m` with `asm5` and `asm10`; then the other 23 with `asm5`, the reverse runs
   for `YH_YH1m` and `AC_SAMEA123396441`, and `YH_LG10Y`. The checks (zip `007`) follow, reviewed against Supplementary file 4
   and Figure 3's supplements 2–7, the seven first.
4. **Graph settings**, from step 3, each recorded in section 6 with its reason.
5. **Graph on the subset** (gate 3): the chromosome runs, then one whole-genome run of the subset; checks as in section 5.
6. **Full panel:** the remaining pairwise alignments (they can run during step 5), reviewed the same way; then the full graph,
   all 25 in one run; gate 4 checks; gate 5 registration, once the coordinator has added a data-root folder for the results.

**Prior art, recorded and not adopted:** Kumar et al. 2025 anchored GT3a's own scaffolds to linkage groups with minimap2 2.28
(through D-GENIES); the alignments behind their Figures S2–S7 used minimap2 2.22 with default settings, after which secondary
matches, blocks under 30% identity or under 10 kb, and contigs with under 2 Mb aligned were removed. They confirmed all six
inversions in the YH assemblies. Quah et al. 2025 built a minigraph (structural-variant) graph of 8 assemblies; the Human
Pangenome Reference Consortium used Minigraph-Cactus and PGGB at 90 haplotypes.

The plan (outputs, formats, checks) is recorded here **before any code is written**: v0.4 records it, before the first code (zip `004`).

---

## 4a. Working with a student: how the assistant behaves, and the checkpoint questions

**How the assistant behaves**
1. It treats the user as a student under the owner's supervision (section 0). It asks the question for a gate before giving the
   commands for that gate, and asks the student to paste evidence (the command, the log lines, a version or `ls` output) before
   moving on.
2. If an answer is wrong or vague, it points to the passage of the documentation or the briefing first. After a second miss it
   explains. It never lets a "yes" stand without the evidence.
3. It explains every command before the student runs it: what it does, what the output should look like, how to undo it. Scripts
   take throttling options (`--threads`, `--max-jobs`); long jobs run inside `screen` on the shared server; one person works on a
   repository at a time.
4. It never invents. If it is unsure whether an option exists or what a default is, it says so and has the student check the
   documentation of the **installed** version.
5. It records a two-line summary of the student's answers at the end of each gate, so the owner can review them.
6. Decisions that belong to the owner (the tools, the subset, the checks and when they are good enough, the output layout) go into
   section 7, not into code.

**Checkpoint questions.** The hints are what a good answer contains; the assistant does not show them before the student has tried.

*Gate 0: concepts, before any software*
1. What are the twelve columns of a PAF line, and what does the `cs` tag add? Why must it be written when the alignment is made?
   *(Query and target name, length, start, end; strand; block and mapping quality; the `cs` string lists matches, substitutions and
   indels; a tag that was not written cannot be rebuilt, and variant calling from assembly alignments needs it.)*
2. What does a pangenome graph represent that a set of pairwise PAFs against one reference does not? *(Relations between the
   non-reference assemblies, and sequence that the reference lacks, such as private insertions.)*
3. Quah et al. built their graph with minigraph. What level of variation does that represent, and what do Cactus (in Minigraph-Cactus)
   and PGGB add? *(Structural variants only; base-level alignment, so SNVs and indels too.)*
4. What does "clipped" mean in Minigraph-Cactus, what is the default, and why does our specification say unclipped? *(Unaligned
   stretches over 10 kb are removed by default; our large private insertions run to about 250 kb; the reference is never clipped.)*
5. How does an inversion appear in a PAF, and why can a single collinear alignment not describe it? *(Alignment blocks on the
   opposite strand, split at the breakpoints.)*

*Gate 1: inputs*
6. Which files do you use, where are they, and how do you locate them? *(`genome.fna.gz` under the Genomes folder by assembly id and
   accession version, through `ga_common.paths`; never a hand-built path.)*
7. What are the names of the sequences the tools will see, and how will you keep track of which reference chromosome each belongs to?
   *(RefSeq or GenBank names as received; two label styles; `provenance.json` and `config/panel_chromosomes_P1.tsv`; the
   plain-number labels are an assumption to be tested.)*
8. Which assembly is the reference, what is its accession, and what is unusual about its LG2 and LG3? *(`MZ_NMK2024a`,
   GCF_041146795.1; LG2 is 48.8 Mb against 38.5–42.3 Mb elsewhere; LG3 is repetitive and varies in length.)*

*Gate 2: running the pairwise alignment*
9. Which minimap2 preset do you use, and why? What does the documentation of the installed version say about divergence? *(Check
   the installed version; the figures differ between versions.)*
10. What does your PAF contain where an assembly carries an inversion? *(Look at it; do not assume.)*
11. How many threads, and why? How do you limit them on the shared server, and how does the job survive a lost connection? Where do
    the outputs go while the job runs? *(`-t`, the throttling flag, `screen`; scratch first, then the data root.)*

*Gate 3: the graph*
12. Which graph builder, and why? What did your comparison on the subset show? *(The numbers from the comparison, not preference.)*
13. What do `--clip` and `--filter` do (Minigraph-Cactus), or `-p`, `-s` and `-n` (PGGB)? *(Check the documentation.)*
14. How do you make the run restartable, and how do you keep memory and disk under control on a shared server? *(The builder's job
    store or per-chromosome runs; throttling.)*

*Gate 4: checking*
15. How do you check that each assembly's path in the graph reproduces its sequence? *(Compare the path-spelled sequence with the
    input; report any difference.)*
16. Which known features do you check, and what do you expect to see? *(Briefing C1: the six inversions per species, the reference's
    LG2 and LG3, the short LG20, the plain-number labels.)*
17. Is a discrepancy an assembly error or a real inversion, and who decides? *(Project 2b decides; this project reports and flags
    through the coordinator.)*

*Gate 5: before registering*
18. What goes in the manifest, and why is a version folder never overwritten? *(Panel version, assemblies and reference, tool
    versions, commands, parameters, date; a version is a record of what was done.)*
19. Which fields do you propose for the `alignment` and `graph` contracts, and what do you need from the coordinator? *(A proposal;
    the coordinator adds the output layout and the validators to the common library once it is agreed.)*

---

## 5. Checks

Written once and applied to every run, in four kinds. Thresholds for the measured kind are set from the subset's numbers
before the full run, relative to the same chromosome in the other assemblies.

1. **Exact: 100%, or the output is not accepted.** Query and target names and lengths match their `provenance.json`; every
   panel assembly has a PAF; every `cs` string reproduces its aligned query segment and agrees with the CIGAR, its matched
   bases less N facing N equal column 10, and its columns minus the ambiguous bases (the `nn` tag) equal column 11, which leaves them out; every graph path
   spells its input (our own script, case-insensitive, differences reported with positions); one complete, forward, acyclic
   reference path per reference sequence; every output validates against its contract version (once defined).
2. **Accounting:** every input base is aligned or in the graph, or listed as missing with a reason; our count from the GFA
   matches Cactus's `clipped-by-reference-contig.tsv` exactly.
3. **Known features** (decided 2026-11-24). Units and breakpoint intervals: Kumar et al.'s Supplementary file 4 (2, 9, 10, 11a,
   11b, 13, 20a, 20b). Expected states: their Supplementary file 2, in `data/expected_inversions.tsv` (an assembly's own genotype
   where its individual was genotyped, otherwise its species', same sex first; the single LG11 genotype applies to 11a and 11b).
   Blocks: Kumar et al.'s four filters (primary; identity at least 30%; at least 10,000 bp; from sequences with at least
   2,000,000 bp of such blocks). A unit is inverted when at least 80% of its interior's aligned bases on its main sequence lie on
   the strand opposite to that sequence's strand in the flanks next to the units (within 2 Mb before the first unit and after
   the last on the chromosome, with at least 0.5 Mb aligned there), normal when at most 20% do, and otherwise unresolved;
   unresolved too when the interior and the flanks lie on different sequences. A segment belongs to a unit's inverted stretch
   only if most of it lies inside the unit. For YH, whose alignments Supplementary file 4 came from (its caption), an inverted
   unit's observed breakpoint intervals must overlap the published ones; for the other assemblies they are reported, with the
   distance of the block ends from the published intervals. Scored on the subset, 55 units: a disagreement
   sends the plan back to question 1; an unresolved unit is reported, not failed. Heterozygous and unknown states are not scored
   (`YH_YH1m`'s LG10, heterozygous in this male). Reported, not scored: the other 17 assemblies against the same expectations;
   `YH_LG10Y` on LG10 (expected normal: the X inverted, the Y normal); the reference's LG2 left arm; LG3; the short chromosomes
   (LG20 in `MC_MC3m` and `CV_CV2f`, LG10 in `YH_YH2f`, and `MZ_UMD2025`); label assignments; the B chromosome in `MZ_UMD2025`.
4. **Measured:** coverage and identity per chromosome (both the BLAST-like identity, column 10 / column 11, and the
   gap-compressed divergence, the `de` tag); reciprocity; each assembly's share of the graph per chromosome; graph length per
   chromosome; time, memory and disk. Flagged, not scored: a chromosome whose coverage is more than 5 points below the subset's
   median, or whose divergence is more than twice it.
5. **Graph checks** (agreed 2026-11-25, before any graph result was read; `scripts/check_graph.py`). Exact: every walk in the
   full and clip graphs spells its input contig over its stated range, compared in upper case, with the first differing
   position; the reference has one walk per reference contig, complete, forward only, and visiting no node twice.
   Accounting: per input contig, its bases in each graph plus Cactus's absent intervals for that graph equal its length, and
   our counts per contig and per reference contig and genome equal Cactus's `clipped-by-input-contig.tsv` and
   `clipped-by-reference-contig.tsv` exactly. Known features: each walk is projected onto the reference walk (steps on
   reference nodes, merged into blocks where the strand holds and both coordinates continue within 10 kb), then called by
   the PAF check's rules, and scored on the subset, in the **full** graph: the clip phase severs some inversion junctions
   (LG20: 4 of 16 reverse-strand runs). Reported: the graph's call and breakpoints against the PAF's; reference-only
   stretches (at least 10 kb) against those no subset PAF covers, flagged below 90% either way, with Cactus's `refgaps.bed`
   alongside; and the measured figures (graph size, bases per genome and graph, clipping by reason, time, memory, disk).

Graph and PAF must agree on the inversion intervals, the reference-only stretches and the contig-to-chromosome assignment.
Outputs name the panel version; commands, tool versions and image checksums are recorded. Disagreements are flagged to
project 2b through the coordinator; deciding between assembly error and biology is 2b's. What could not be checked is reported.

---

## 6. Decisions

| date | decision | reason |
|---|---|---|
| 2026-10 | Pairwise alignment to the reference plus an unclipped graph | coordinator decision |
| 2026-10-26 | Whole genomes | coordinator decision |
| 2026-10-26 | Where the tools run is this project's choice | the owner wants flexibility |
| 2026-10-26 | The first conversation starts at the briefing (prior work, then questions 1 and 2); no installation or run before the owner agrees the plan | the owner wants the design questions settled before compute |
| 2026-10-26 | The assistant works with a student under the owner's supervision, with checkpoint questions (section 4a) | the owner's request |
| 2026-10-26 | The student uploads PDFs of the reference papers and pages (briefing Part D); the assistant confirms it can see items 1–8 before going on | uploaded copies are always present and fixed to a version |
| 2026-11-04 | Graph builder: Minigraph-Cactus; the deliverable is the unclipped (full) graph | the owner with the student: the design is anchored on the reference, matching "all results in reference coordinates", and the pipeline reports exactly what is missing from each graph |
| 2026-11-04 | Pairwise aligner: minimap2 | the owner with the student: the only documented option that writes the `cs` tag the `alignment` contract requires |
| 2026-11-04 | The comparison of builders in briefing Part B3 is dropped | the owner's decision |
| 2026-11-17 | The owner leaves the decisions on questions 1 and 2 to the student | the owner's decision |
| 2026-11-17 | Methods, subset and checks as in sections 4 and 5 | the student's decision |
| 2026-11-17 | Preset `asm5` for all 24, with one `asm10` run on `YH_YH1m` as a sensitivity check | `asm5` fits the panel's low divergence; the check covers the breakpoints and LG3, where divergence is local and higher |
| 2026-11-17 | Environment: conda env `ga-alignment` with Apptainer 1.5.3 and squashfuse 0.6.1, Python 3.11, minimap2 2.31 and `ga_common` at `common-v0.2.1`; Cactus v3.3.0 image in `~/containers` | project 1's tested recipe (base's Apptainer 1.2.2 cannot mount images); one pinned version of each tool for all outputs |
| 2026-11-18 | Pairwise alignment of all 24 assemblies can start | the owner's confirmation |
| 2026-11-19 | Cactus v4.0.0 instead of v3.3.0 | v3.3.0's default configuration has neither `zipAlleles` nor the exact overlap filter, which the inversion checks rely on; 4.0.0 has both on, matches the documentation the plan was written from, fixes a minigraph mapping bug present since v2.1.0, gives reproducible output and adds `--validate` |
| 2026-11-19 | `YH_LG10Y` is aligned too, after the 24 panel assemblies | the student's decision |
| 2026-11-20 | Pairwise order: `YH_YH1m` with `asm5` and `asm10` first, then the other 23 with the chosen preset, and the reverse runs | the student's decision; settles the preset before the bulk |
| 2026-11-20 | 16 threads per run, two runs at once and then six, under `nice -n 10`; version folder `PW_v1` | the student's decision; minimap2 works per query sequence, so more runs use the machine better than more threads per run |
| 2026-11-20 | The seven assemblies of section 4, including both of Kumar's *M. zebra* controls and `DT_SAMEA123396444`, are reviewed first and form the graph's subset | the student's decision |
| 2026-11-24 | `asm5` confirmed for all 24, after the comparison on `YH_YH1m` (decided before the bulk runs; recorded here) | `asm10`'s additional alignment is about 7% divergent (section 4) |
| 2026-11-24 | Known features scored as in section 5: Kumar et al.'s four block filters; the 80% strand rule; breakpoint intervals scored by overlap with Supplementary file 4 | the method of Kumar et al.'s Version of Record; real blocks end short of breakpoints, which lie in repeats (`YH_YH1m`, LG9) |
| 2026-11-24 | Expected states from Supplementary file 2; `YH_YH1m`'s LG10 not scored; DT's LG11 scored as both 11a and 11b | genotypes for six of our individuals and all eight species; LG10 is XY in YH; LG11 is one genotype in Supplementary file 2 |
| 2026-11-25 | The graph checks as in section 5 (item 5); the whole-genome subset run waits for their results on LG9 and LG20 | agreed with the student; a settings problem found on two chromosomes costs less than one found on the whole genome |
| 2026-11-24 | Graph settings and caps as in section 4 (step 4) | agreed with the student after gate 3; options and defaults read from the 4.0.0 help and configuration |
| 2026-11-24 | The full graph is built from scratch at step 6 | section 4: little is reused, reusing mappings needs `zipAlleles` off, and the result would depend on the subset |
| 2026-11-24 | Inversion calls take the background strand from the flanks next to the units (2 Mb on each side, at least 0.5 Mb aligned there) instead of the whole chromosome; a segment joins a unit's inverted stretch only if most of it lies inside the unit. Changed after the first results (section 8) | the whole-chromosome rule let misplaced terminal pieces of up to 17.4 Mb outvote the true neighbours (YH and `MC_MC2f` on LG20, `YH_YH7f` on LG9); reading the immediate flanks is how a dot plot is read. Previewed on the real segments: 55 of 55 scored units agree, at 1, 2 and 5 Mb windows alike |
| 2026-11-24 | Breakpoints scored against Supplementary file 4 only for YH; reported for the other assemblies. Changed after the first results | the table's caption: the intervals were calculated from alignments of YH to GT3a; DT's first LG9 breakpoint lies 0.37 Mb from YH's in all three DT assemblies |
| 2026-11-24 | Breakpoint bounds stay as written in zip `007`; a proposed clipping of flank segments that overlap an inverted stretch was withdrawn | previewed on the real segments, it moved `YH_YH2f`'s first LG9 breakpoint out of the published interval, which it had matched |
| 2026-11-24 | Supplementary files 4 and 2 kept in `published/` as recorded inputs, each with its source and sha256, and pinned by `tests/test_published.py` | the checks read recorded inputs, not numbers in code; `data/` is kept out of the repository by its `.gitignore` (data lives under the data root) |
| 2026-11-24 | Soft-masking is not an input | Minigraph-Cactus ignores it (its FAQ); minimap2 aligns masked sequence (89.5% of the reference covered on `YH_YH1m`, of which only 68% is upper case) |
| 2026-11-19 | Pairwise: minimap2's default secondary alignments, every record kept; the reference indexed in each run; PAFs in scratch until a data-root folder exists | the student's decision; a saved index would fix one preset's indexing parameters (minimap2 README) |

---

## 7. Open questions

| question | blocking |
|---|---|
| Exact fields of the `alignment` and `graph` contracts: short or long `cs`; all PAF records or a filtered set; sorted or not (coordinator, with project 3) | before first integration |
| **Output layout:** folder, version naming and manifest for the alignment and graph outputs (the common library knows only `Genomes/` and `RepeatInformation/`); the student proposes, then the coordinator adds the paths | before the first output is registered |
| Data-root folders for the results, added to `ga_common.paths` and `sync.py` by the coordinator; until then results stay in `~/scratch/ga-alignment` | before results leave scratch |
| Readers and validators for PAF and GFA, and coordinate-frame conversions, to be added to the common library by the coordinator once the fields are agreed | before first integration |
| How chromosomes are anchored across assemblies given mixed naming and submitter-specific LG labels (section 2): tested by the pairwise alignments (section 4), with project 2b | no |
| A shared place for the published tables, if project 2b needs them too (coordinator) | no |
| Whether `YH_LG10Y` joins the full graph | step 6 |
| Misplaced terminal pieces, each at the start of a chromosome on "+" although it belongs at the far end next to its reference neighbour: `YH_YH1m` LG9 (reference 0–6.2 Mb), `YH_YH7f` LG9 (0–9.7 Mb), and LG20 of `YH_YH1m`, `YH_YH2f`, `YH_YH7f` and `MC_MC2f` (0–17.4 Mb). The identical placement in four Georgia Tech chromosome builds suggests how they were built, not biology: to project 2b through the coordinator | no |
| DT's first LG9 breakpoint, 0.37 Mb from YH's in all three DT assemblies; `YH_YH7f`'s first LG9 breakpoint, where flank and inverted segments overlap (14.71–14.94 Mb) around the published interval; `CV_CV4f`'s LG20, normal on 20a and unresolved on 20b against CV's 26 of 26 inverted: to project 2b | no |
| Note for the owner: the Version of Record's figure captions swap the YH sexes (section 2) | no |
| `MZ_UMD2025`'s B chromosome (section 2): its alignments may look like rearrangements; projects 2b and 3 should know (coordinator) | no |
| The briefing says mzebra has no docker or singularity; Apptainer runs there from a project conda env (project 1's recipe). A correction for the coordinator | no |

---

## 8. Status

- Questions 1 and 2 decided (sections 4 to 6) by the student, to whom the owner left them (2026-11-17).
- **Gate 0** (2026-11-12): the PAF columns and the difference between minigraph and Cactus were correct first time. What a
  graph adds over pairwise PAFs, and what clipping removes, were correct on the second try. Explained after two attempts: what
  `cs` records beyond the CIGAR and why it must be written at alignment time; that clipping cannot remove deletions (that is
  `--delFilter`); and how an inversion appears in PAF records, to be rechecked on real output at gate 2, question 10.
- **Gate 1** (2026-11-19): where the genomes are, and that only the first word of each header reaches the tools, were right
  first time; locating files through `ga_common.paths` was right on the second try. Explained after two attempts: that chromosome
  labels are the submitters', in two styles, and are confirmed by alignment rather than taken from the records; that the
  reference's file carries RefSeq names under a folder named by its GenBank accession; and the LG2 and LG3 length observations.
- The environment is ready (2026-11-19): Apptainer 1.5.3, Python 3.11.18, minimap2 2.31-r1302, `ga_common` 0.2.1, and the Cactus
  4.0.0 image (section 4); mzebra has 96 CPU threads with x86-64-v3 and 3.7 TB free in home.
- **Gate 2** (closed after `YH_YH1m`'s runs of 2026-11-20, before the other 23): question 9 right in substance after two corrections (the divergence figure is 0.1–0.25%;
  the `asm10` run tests the alignment around structural variants). Question 11 explained: the throttle flags, `nice`, checking
  the load, `screen`, scratch first, and why several 16-thread runs use the machine better than one with 96. Question 10, on
  real output: the LG9 inversion and both breakpoint intervals found and matched to Table S4. Corrected: an inversion is marked
  by its strand relative to its flanks. Explained: the misplaced first segment of `YH_YH1m`'s LG9, a likely scaffolding error.
- Inputs checked (2026-11-19): every check passed (section 2).
- Pairwise alignments done (section 4): 28 runs in `PW_v1`.
- Zip `007`'s first check run (`checks/20261125-024645`) failed 27 of the 28 runs on one comparison only. Its rule took column
  11 to count every alignment column, but column 11 leaves out ambiguous bases, which minimap2 counts in the `nn` tag: the
  differences were 100, 200 or 300 bases, NCBI's 100-N runs for assembly gaps. Every `cs` string rebuilt both genomes base
  for base and agreed with its CIGAR; the reference has no letters other than A, C, G, T and N. The same zip's two tables
  were left out of git by `.gitignore` (`data/`), so the summary could not start. Zip `008` fixes both; the checks rerun.
- Zip `008`'s check run (`checks/20261125-030421`): 24 of 28 runs passed every exact check; the ambiguous bases counted by
  minimap2's `nn` tags (2,474,206) equal the N columns counted from the genomes. The four failures were N facing N inside
  matched runs, which column 10 does not count as residue matches, in three Georgia Tech assemblies aligned to the Georgia
  Tech reference.
- Zip `008`'s summary (`summary/20261125-030423`): 50 of 55 scored units agreed. The five disagreements: `YH_YH1m` and
  `YH_YH2f` on 20a and 20b, from the whole-chromosome background rule; and `DT_SAMEA123396444`'s first LG9 breakpoint.
  `YH_YH1m`'s breakpoint intervals equal Supplementary file 4's to the base on 2, 11a, and the first of 9 and 13. The Y
  haplotig `YH_LG10Y` is normal on LG10 (other strand 0.009), as the genotypes imply (the X inverted, the Y normal). Zip
  `009` changes the rules (section 6, owner to be told); the checks rerun.
- Zip `009`'s checks pass (`checks/20261125-033145`, `summary/20261125-033146`): 28 of 28 runs pass every exact check; 55 of
  55 scored units agree; 131 of 138 reported units agree. Of the 30 scored breakpoints (all YH), 8 observed intervals equal
  Supplementary file 4's to the base and 25 inverted-stretch ends lie inside the published intervals; the other 5 lie up to
  570 kb away. The findings for project 2b are drafted for the coordinator (section 7).
- **Gate 3** (2026-11-24): question 12 named minigraph; corrected to Minigraph-Cactus, with the reasons recorded in section
  6 (no builder comparison: the owner dropped it); completeness and accuracy mapped onto section 5's checks. Question 13:
  the subgraph relation was right; what `--clip` and `--filter` remove was explained; `--filter` is not used. Question 14
  explained: the job store and `--restart`, per-chromosome steps, and the resource options, confirmed in the 4.0.0 help.
- Chromosome graph runs (2026-11-25), each the reference and the seven, all exit 0 with `--validate` passed: LG8 24 min, 2
  CPU-hours, largest process 8.8 GB; LG9 43 min, 5 CPU-hours, 8.7 GB; LG20 33 min, 3 CPU-hours, 7.0 GB; 1.0-1.3 GB of disk
  each. Every chromosome sequence is complete in its full graph. Left out of every graph, as unassigned or ambiguous, were
  only the extra unplaced scaffolds that the pairwise assignment had added (for example 4 of `YH_YH1m`'s, 2.75 Mb, on LG20):
  Cactus's thresholds are stricter than ours. On LG20, the deletion filter removed 589 mapping lines (2.48 Mb of matches)
  and the clip phase severed 4 of 16 reverse-strand runs. Zip `011` checks the three graphs.
- 2026-11-18: an environment zip mistakenly numbered `alignment_001` was committed here (5d9f9e2) and reverted (f60fbea); the
  same content comes as `004`. The same zip had been pushed to `ga-repeats` by mistake (f15e6dc) and was reverted there.
- The repeat annotation (if it is an input) is not yet produced. Nothing is ready to integrate.

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
- Using the coordinator's code: `docs/USING_COMMON.md` in the coordinator repository (install, configuration, examples, what is missing).
- **Documents in the Claude project:** this specification, `docs/alignment_briefing.md`, `genome_panel_P1.tsv`, `species_context.tsv`
  and `panel_chromosomes_P1.tsv` from the coordinator repository, and the papers and pages of the briefing, Part D, **uploaded by the
  student as PDFs** (items 1–8 required, 9 and 10 may follow).
- Python 3.9 is the minimum for code that uses the common library.
- **`VERSION`** holds the number of the last zip applied, three digits (`011`). Every zip's instructions start by checking `pwd`,
  `git remote -v` and `VERSION` before anything is copied.
- **This project on mzebra:** conda env `ga-alignment` (section 6); Cactus image `~/containers/cactus_v4.0.0.sif` (the v3.3.0 image may be removed), with the
  pull cache and tmp in `~/containers/cache` and `~/containers/tmp`, set per session (`.bashrc` is not edited); working space
  `~/scratch/ga-alignment`; long jobs in `screen`; nobody commits or pushes on mzebra.

---

## 10. Relatedness context

See section 2a.
