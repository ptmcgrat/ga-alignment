# Alignment project: briefing for the first conversation

**For:** the student and the assistant of project 2a (Alignment). **Owner:** McGrath lab.
**Version:** 1, 2026-10-26. Written by the coordinator; the project decides its own methods.

This briefing is where the first conversation starts. It has two questions that have to be settled with the
owner before any alignment or graph run, and the prior work to read first.

1. **How should the pairwise alignments and the pangenome graph be built?** Which tools, with which settings, and
   how do we compare the options before spending the compute on all 25 assemblies?
2. **How will we know the alignments and the graph are right?** What do we check them against, given six large
   inversions, a reference with known problems, and assemblies named in two different ways?

Statements below are of three kinds: **what a source says** (with where), **my reading** (the coordinator's inference,
marked as such), and **what is not known**. Nothing has been run yet.

---

## Part A. What the project has to deliver, and the facts it works with

### A1. The deliverables (coordinator specification, sections 2 and 6)

- **Contract `alignment`:** a PAF file per assembly with the `cs` tag, each assembly against the reference
  (M. zebra, `MZ_NMK2024a`, GCF_041146795.1). PAF is 0-based, half-open.
- **Contract `graph`:** one pangenome graph of all assemblies, **unclipped**, in GFA.
- Whole genomes. All results in reference coordinates. Every output names the panel version (P1).
- Consumers: project 2b (breakpoints, arrangements, assembly-error flags) and project 3 (SNVs, indels, transposon
  insertions, small copy-number variants and inversions, polarized). The exact fields of both contracts are for
  this project to propose, with the owner.

### A2. The data (all downloaded, 2026-10-26)

- 25 assemblies of 8 species, all seven ecogroups; contig N50 22–41 Mb. Plus the YH LG10 Y haplotig (`YH_LG10Y`), which
  is held beside the panel.
- **Sequence names are as NCBI delivered them:** RefSeq (`NC_...`) for the 8 annotated assemblies, GenBank (`CM...`,
  `OZ...`) for the rest. NCBI's chromosome labels are the submitters', in two styles: `LG2` (17 assemblies) or a plain `2`
  (8 assemblies). That a plain number is the same chromosome as the reference's LG is a working assumption.
- `config/panel_chromosomes_P1.tsv` gives, per assembly, the accessions labelled LG2, 9, 10, 11, 13 and 20.
- Genome files are BGZF-compressed (`genome.fna.gz`); BGZF is valid gzip, so tools that read `.gz` read it. Whether each
  tool you choose does so has to be checked.

### A3. The biology (coordinator specification, section 4b)

- Species differ by only 0.1–0.25% at SNPs, but 4.7–9.9% of base pairs differ structurally (Quah et al. 2025), and a graph
  of 8 genomes held 33.1% more sequence than one genome. About 75% of the structural-variant sequence is transposon-derived,
  with the largest insertions 17–252 kb.
- **Six large inversions**, each roughly 10–20 Mb, on LG2, 9, 10, 11, 13 and 20 (LG11 and LG20 are tandem doubles; on LG20
  *Rhamphochromis* carries only the first arm). Expected states per species are in `config/species_context.tsv`; they are a
  literature prior, not a result. The reference carries the ancestral state of all six.

### A4. The reference has known problems

- Kumar et al. 2025 (eLife) describe the reference they use, `M_zebra_GT3a`, and say that alignments of *M. zebra* genomes to it
  differ **most prominently in the left arm of chromosome 2 and in the highly repetitive chromosome 3**, which likely indicates some
  errors in GT3a or a large polymorphism within the species.
- That matches what the P1 download shows: the reference's LG2 is 48.8 Mb against 38.5–42.3 Mb in the 24 other assemblies, and
  LG3 varies 1.8-fold across assemblies (58.7–105.3 Mb; reference 64.4).
- **Not yet verified:** that `M_zebra_GT3a` (WGS project JBEVYI000000000 in Kumar et al.) is the same assembly as
  GCF_041146795.1 / GCA_041146795.1. It looks likely, but the student should confirm it from NCBI's assembly report.
- Other length observations to carry in mind (coordinator specification, section 4): LG20 is 17.8 Mb in `MC_MC3m` and 20.3 Mb in
  `CV_CV2f` against 35–42 Mb in most; `YH_YH2f` LG10 is 28.4 Mb against 34–42 elsewhere; `MZ_UMD2025` is short on LG3, LG9, LG11 and LG13.

### A5. The closest prior work, to read first

- **Kumar et al. 2025** aligned five new assemblies to GT3a and found that the alignments support all six inversions in the YH
  assemblies and none in the control *M. zebra* assemblies (their Figures S2–S7). Single scaffolds often spanned a whole inversion.
  They could not place breakpoints at base-pair resolution because of repetitive DNA near them. They used minimap2 2.28, through
  D-GENIES, to anchor scaffolds. This is the result our alignments should reproduce.
- **Quah et al. 2025** built a graph of 8 long-read assemblies with minigraph, with *A. calliptera* as the backbone. That is a
  structural-variant graph, not a base-level graph.

---

## Part B. Question 1: how should the alignments and the graph be built?

### B1. Pairwise alignment to the reference

- **minimap2** is the usual tool for assembly-to-reference alignment. The documented form is
  `minimap2 -cx asm5 --cs ref.fa asm.fa > asm.paf`. The presets `asm5`, `asm10` and `asm20` are for increasing divergence;
  the documentation is not consistent between versions on the figures (below 1% in the cookbook, about 0.1% in a help text), so
  the student reads the documentation of the **installed** version.
- **The `cs` tag must be written at alignment time.** A tag that was not written cannot be rebuilt later, and `paftools.js call`
  (variant calling from assembly alignments) needs it, along with sorted PAF. Project 3 will care.
- **Inversions:** a collinear aligner reports an inversion as separate alignment blocks on the opposite strand. What exactly the PAF
  contains there, and what the aligner does at the boundaries, is something the student has to look at, not assume.
- **wfmash** (the aligner inside pggb) is a second option: it is described as aligning gigabase genomes in minutes to hours, with a
  segment length (`-s`) that makes alignments more collinear when it is set longer.

### B2. The graph: the two builders named in the specification

| | Minigraph-Cactus (`cactus-pangenome`) | PGGB |
|---|---|---|
| principle | reference-first: minigraph builds a structural-variant backbone from the reference, the other assemblies are mapped to it, and Cactus adds base-level alignment | reference-free: all-to-all alignment with wfmash, graph induced with seqwish, normalised with smoothxg and gfaffix |
| fit to our case | documented as suitable for similar samples, such as the same species; needs a reference (we have one) | designed for any set; key settings are the minimum identity `-p`, the segment length `-s` and the number of haplotypes `-n` (defaults `-p 95 -s 10k`) |
| clipping | by default it removes unaligned stretches over 10 kb (`--clip 10000`) and low-support nodes (`--filter 2`); it emits clipped, filtered and full graphs; the full, unclipped graph is made with `--clip 0 --filter 0`; the reference is never clipped | not applicable in the same way; the graph follows the alignments |
| scale evidence | 90 human haplotypes (HPRC) | 20–90 haplotypes (HPRC, mouse) used about 150 GB at peak on a 386 GB machine (a developer's answer in a GitHub issue); the same thread expected about a day for a two-genome case, and the project's issue threads discuss running by chromosome |
| how it runs | a Toil workflow; docker or singularity image, or a binary release | conda (bioconda) or containers |
| outputs | GFA, vg, HAL, VCF | GFA, odgi |

Readings (the coordinator's, to be challenged):
1. **Clipping matters to us.** The specification asks for an unclipped graph. Quah et al. found large private insertions, up to
   252 kb. Minigraph-Cactus's default clipping would remove unaligned stretches over 10 kb, so the default graph would lose exactly
   those. The unclipped (`--clip 0`) graph, or the "full" output, is what the specification means. The HPRC graphs also excluded highly
   repetitive sequence such as centromeres; the student checks what that implies for the repeat-rich chromosome 3.
2. **A base-level graph is needed.** Project 3 calls SNVs, indels and small variants, so a structural-variant-only (minigraph) graph is
   not enough on its own.
3. **Chromosome anchoring is tied to the builder.** Both approaches work per chromosome at this scale, which brings back the
   naming problem in A2. How each tool assigns assembly sequences to reference chromosomes has to be read in its documentation.
4. **Environment.** mzebra has conda, 96 cores and about 500 GiB of memory, shared, with no docker or singularity. Which builders can
   run there (binary release, conda) is for the student to check.
5. **Inversions in the graph.** Settings such as the segment length change how long collinear blocks are represented, so how the
   six inversions appear in each graph has to be tested, not assumed.

### B3. How to decide (to be designed with the owner, measured not assumed)

Run both builders (and the pairwise aligner) on a small, deliberate subset before committing the full panel. The subset is the student's
to propose and the owner's to agree. It should include the reference, an assembly with no inversions (for example *A. calliptera*), one
with several (for example YH, which carries all six), and one with the LG20 left-arm-only state (*Rhamphochromis*), and probably only the
chromosomes that carry inversions. Report at least: time, memory and disk; whether each known inversion is represented as expected; the
share of each assembly that is in the graph; and what happens to large insertions and to the repeat-rich chromosome 3. Then decide, and
record the decision with the numbers in the specification (section 6).

---

## Part C. Question 2: how will we know the results are right?

**Short answer: against known features, with measurements, and with an honest account of what could not be checked.** Nothing has been run.

### C1. Known features to check against

| feature | what to look for |
|---|---|
| six inversions (LG2, 9, 10, 11, 13, 20) | present where the literature prior says (YH: all six; *Diplotaxodon*: LG9, 11, 20; *Rhamphochromis*: LG20 left arm only), absent in the ancestral set (MZ, AC). Kumar et al.'s Figures S2–S7 are the comparison. |
| reference LG2 (48.8 Mb against 38.5–42.3 Mb elsewhere) | the alignments should show where the reference carries extra or rearranged sequence; do not read it as biology before 2b has looked |
| reference LG3 (repetitive; length varies 58.7–105.3 Mb across assemblies) | alignment coverage and identity there, and what is clipped or missing from the graph |
| LG20 short in `MC_MC3m` (17.8 Mb) and `CV_CV2f` (20.3 Mb) | whether the rest of LG20 sits elsewhere in those assemblies |
| `YH_YH2f` LG10 (28.4 Mb) | against the male's 36.5 Mb |
| mixed names and label styles | every assembly sequence ends up on the right reference chromosome; the plain-number labels are confirmed or refuted |
| `MZ_UMD2025` (short on LG3, 9, 11, 13) | whether this reflects its assembly |

### C2. Checks that apply to every output

1. **Names and lengths:** the query and target names and lengths in every PAF match `provenance.json` of the assembly and of the reference.
2. **Coverage and identity by chromosome:** reported per assembly, with a plot or table, not a single number.
3. **`cs` consistency:** the `cs` strings reproduce the aligned sequence.
4. **Reciprocity:** an alignment of the reference to an assembly agrees with the assembly to the reference, at least on a subset.
5. **The graph:** each assembly's path through the graph spells out its sequence; the reference is complete and unclipped.
6. **Reproducibility:** commands, versions and parameters are recorded so the run can be repeated.

### C3. What is not the alignment project's call

Telling an assembly error from a real inversion belongs to project 2b. Alignment reports what the data show and flags what it
cannot explain; the flags travel through the coordinator.

---

## Part D. Sources, in reading order

1. **Kumar et al. 2025**, eLife (doi 10.7554/eLife.104923). The lab's copy is `KumarElife.pdf` in the coordinator project.
   *Read first: the six inversions, the reference, and the alignments in Figures S2–S7.*
2. **Quah et al. 2025**, Genome Research 35:1094. https://genome.cshlp.org/content/35/5/1094 (open access).
3. **Blumer et al. 2025**, Science (adr9961). The lab's copy is `science_adr9961.pdf` in the coordinator project. Inversions and TE
   insertions typed in the radiation.
4. **Hickey et al. 2024**, the Minigraph-Cactus paper. https://pmc.ncbi.nlm.nih.gov/articles/PMC10638906
5. **Minigraph-Cactus documentation**: https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md
6. **PGGB documentation**, essential parameters: https://pggb.readthedocs.io/en/stable/rst/essential_parameters.html ,
   and the repository https://github.com/pangenome/pggb
7. **wfmash README**: https://github.com/waveygang/wfmash
8. **minimap2 README and cookbook**: https://github.com/lh3/minimap2
9. **HPRC Minigraph-Cactus resources** (what clipping removed): https://github.com/human-pangenomics/hpp_pangenome_resources/blob/main/hprc-v1.0-mc.md
10. **The pggb paper**, bioRxiv 2023.04.05.535718 (https://doi.org/10.1101/2023.04.05.535718), optional.

### Uploading the sources (the student does this, before the first conversation)

The student saves each source **as a PDF** and uploads it to the Alignment project, so the assistant has the same fixed text at the start
of every session. Items 1–8 are needed before the first conversation continues; 9 and 10 may follow.

| # | how to get the PDF | suggested file name |
|---|---|---|
| 1 | the lab's `KumarElife.pdf` (ask Professor McGrath), or the journal's PDF | `Kumar2025_eLife_inversions.pdf` |
| 2 | on the journal page, "Download PDF" | `Quah2025_GenomeRes_pangenome.pdf` |
| 3 | the lab's `science_adr9961.pdf` (ask Professor McGrath) | `Blumer2025_Science.pdf` |
| 4 | the PDF link on the PubMed Central page | `Hickey2024_MinigraphCactus.pdf` |
| 5 | browser: Print, then Save as PDF | `Cactus_pangenome_doc.pdf` |
| 6 | browser: Print, then Save as PDF (the page, and the repository's README) | `PGGB_parameters.pdf`, `PGGB_README.pdf` |
| 7 | browser: Print, then Save as PDF | `wfmash_README.pdf` |
| 8 | browser: Print, then Save as PDF (README and cookbook) | `minimap2_README.pdf`, `minimap2_cookbook.pdf` |
| 9 | browser: Print, then Save as PDF | `HPRC_MC_resources.pdf` |
| 10 | the bioRxiv page, "Download PDF" | `pggb_paper.pdf` |

How:
- **Papers:** use the publisher's or PubMed Central's own PDF, not a print of the web page.
- **Web pages:** Print, then Save as PDF. Open the result and check that the text can be selected (it should not be a picture of the page) and that the whole page is there.
- **Upload:** add the PDFs to the Alignment project's files.
- **Check:** ask the assistant, "List the documents in this project and give the title of each." Re-upload anything missing or unreadable.
- **Versions:** when the tools are available, run `minimap2 --version` and the equivalent for the graph builder, and tell the assistant whether the uploaded documentation matches the installed version.

---

## Part E. Questions this raises for other projects (relayed through the coordinator)

- **For the coordinator:** the folder layout, version naming and manifest for the alignment and graph outputs (the library knows only
  `Genomes/` and `RepeatInformation/`); readers and validators for PAF and GFA; coordinate-frame conversions. The student proposes
  with the owner; the coordinator then adds them to the common library.
- **For project 2b:** which flags it needs from the alignments, and which chromosome anchoring it expects.
- **For project 1:** whether soft-masked genomes are an input to the alignment at all (the specification says so; this project decides).
- **For the owner:** confirm that `M_zebra_GT3a` and GCF_041146795.1 are the same assembly.
