# Alignment specification

**Project:** 2a, Alignment
**Repository:** `ga-alignment` *(name assumed)*
**Document version:** 0.3, 2026-10-26 (0.1: first draft from coordinator specification v0.7; 0.2: P1 downloaded and chromosome labels, from v0.8; 0.3: the first conversation starts at the briefing; working with a student; using the common library, from v0.10)
**Owner:** McGrath lab
**Status:** draft for discussion with the owner; no code written; **start at section 0 and `docs/alignment_briefing.md`**

This is the single specification document of project 2a. It is always replaced whole. It records what the
coordinator specification already fixes and leaves everything else open. The methods are not designed
here: the project designs them with the owner.

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

**Not designed yet.** The two questions to settle first are in `docs/alignment_briefing.md` (Parts B and C).
Fixed by the coordinator specification: pairwise alignment of every assembly to the
reference; an unclipped pangenome graph; whole genomes (decision 2026-10-26); all results in reference
coordinates. Not decided: the graph builder (Minigraph-Cactus or PGGB, coordinator open question), the
pairwise aligner and its settings, and how large inversions are handled in the pairwise step.

**Prior art, recorded and not adopted:** Kumar et al. 2025 aligned new assemblies to the reference with minimap2 2.28 (through
D-GENIES) and confirmed all six inversions in the YH assemblies; Quah et al. 2025 built a minigraph (structural-variant) graph of
8 assemblies; the Human Pangenome Reference Consortium used Minigraph-Cactus and PGGB at 90 haplotypes.

The plan (outputs, formats, checks) is agreed with the owner and recorded here **before any code is written**.

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

## 5. Checks (proposed; to be agreed)

- every output validates against its contract version (once defined);
- every panel assembly has a PAF; query names and lengths match the assembly's `provenance.json`;
- target names and lengths match the reference's `provenance.json`;
- outputs name the panel version; commands and tool versions are recorded;
- the checks of the briefing (Part C2) and the known features (Part C1) are reported, with what could not be checked.

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

---

## 7. Open questions

| question | blocking |
|---|---|
| **Question 1:** which pairwise settings and which graph builder, and the comparison on a subset (briefing, Part B)? The student proposes; the owner agrees | **yes, before any run** |
| **Question 2:** what checks, against which known features, and what counts as right (briefing, Part C)? | **yes, before any output is accepted** |
| Is the repeat annotation (soft-masking) an input to the alignment? | no |
| Exact fields of the `alignment` and `graph` contracts | before first integration |
| **Output layout:** folder, version naming and manifest for the alignment and graph outputs (the common library knows only `Genomes/` and `RepeatInformation/`); the student proposes with the owner, then the coordinator adds the paths | before the first output is registered |
| Readers and validators for PAF and GFA, and coordinate-frame conversions, to be added to the common library by the coordinator once the fields are agreed | before first integration |
| Is `M_zebra_GT3a` (Kumar et al.; WGS project JBEVYI000000000) the same assembly as GCF_041146795.1? (check NCBI's assembly report) | no |
| Which builders can run on mzebra (no docker or singularity): binary release or conda | no |
| How chromosomes are anchored across assemblies given mixed naming and submitter-specific LG labels (section 2) | no |
| Is the YH LG10 Y haplotig (`YH_LG10Y`, held beside the panel) aligned too? | no |
| Where the tools run, and how long jobs are kept alive on the shared server | no |

---

## 8. Status

Not started. Specification and briefing drafted; the P1 genomes are available. The first conversation starts at section 0; the
student's proposals for questions 1 and 2 are awaited. The repeat annotation (if it is an input) is not yet produced. Nothing is ready to integrate.

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

---

## 10. Relatedness context

See section 2a.
