# CLAUDE.md — latinoRDF

Reference for contributors (human and AI) working on **latinoRDF**, the CMS
RDataFrame-based analysis framework for histogram production, plotting and
datacard creation.

Initialize the environment first so the scripts are discoverable:

```bash
source setup.sh
```

---

## 1. Repository layout

| Path | Contents |
|---|---|
| `scripts/` | The executable framework: `mkShapes.py`, `mkSnapshot.py`, `analysisFactory.py`, `mkPlot.py`, `mkDatacards.py`, `mkMergeYears.py`, `haddfast.py`, `haddnano.py`, `easyDescription.py` |
| `lib/` | Plot styling helpers (`tdrStyle.py`, `CMS_lumi.py`) and `search_files.py` |
| `plotsConfigurationsRDF/` | **Separate git repository** (`acmbulla/plotsConfigurationsRDF`), ignored by this repo via `.gitignore:209`. Holds the per-analysis, per-year configuration files. |

> **Important:** `plotsConfigurationsRDF/` is a *different repo*. Changes to
> analysis configs (including the `gInterpreter` removal, see §7) are committed
> and reviewed there, **not** in a latinoRDF branch.

### Searching the tree

The configured shell wraps `grep` with `--ignore-files`, which honours
`.gitignore`. Because `plotsConfigurationsRDF/` is gitignored, a plain
`grep -r <pattern> .` from the repo root **silently finds nothing there**.
Use `command grep -r ...` (or `git grep` inside the nested repo) when searching
configuration files.

---

## 2. Scripts and their purpose

### `analysisFactory.py` — shared base class
Defines `AnalysisFactory`, the base class for both producers. It owns everything
common up to and including `current_node` + cuts:

- the ordered aliases/functions pipeline and C++ code generation
  (`_build_aliases_ordered`, `_build_cpp_common`)
- weights and `supercut` handling
- cuts / categories
- job chunking: `MaxEventsPerJob` / `FilesPerJob`, with parallel event counting
  and a JSON cache (`samples_metadata.json`)
- `split` / `era` event slicing across eras via RDataFrame `Range()`
- `submitBatch` (proxy check, pre-compiled vs compile-on-node, single `.sub` file)
- `checkBatch` and `--resubmit`
- `parallelCompile`

It also provides the config loaders used by every script:
`load_config_files`, `read_exec`, `defaultParser`, `loadDefaultOptions`.

Subclasses must implement `job_prefix()`, `output_file_name()`, `makeScripts()`,
`generate_makefile()`, `_build_cpp_output()`, `_build_cpp_header()`,
`_check_output_exists()`, `_build_bash_precompiled()`,
`_build_bash_compile_on_node()` and `_extra_transfer_files_for_submit()`.

### `mkShapes.py` — histogram production (`ShapeFactory`)
Generates and runs the C++ that fills histograms for every
sample × cut × variable × nuisance. Reads nanoAOD directly, or existing
snapshots via `makeFromSnapshot()` (`--fromSnapshot`).

Options: `--parallelCompile`, `--localCompile`, `--fromSnapshot`,
`--submitBatch`, `--hadd`, `--haddBatch`, `--checkBatch`, `--resubmit`,
`--silentMode`.

### `mkSnapshot.py` — skimmed-tree production (`SnapshotFactory`)
Writes reduced nanoAOD-style trees (snapshots) so later shape production can
skip the expensive full-nanoAOD pass. Tracks what was produced via
`SnapshotMetadata`, and supports `--patch` to top up an incomplete snapshot set.

Options: `--localCompile`, `--submitBatch`, `--checkBatch`, `--resubmit`,
`--silentMode`, `--haddBatch`, `--patch`, `--shapeNuisances`, `--force`.

### `mkMergeYears.py`
Merges per-year shape files into a combined multi-year result, handling
correlated vs uncorrelated nuisances across years. Splits by cut for batch
submission. See §6 for the workflow, the measurements behind its shape, and the
copying limit that caps how far it can be split.

### `haddfast.py`
Fast parallel merger for **shape/histogram** ROOT files. Used by
`mkShapes.py --haddBatch`.

### `haddnano.py`
nanoAOD-aware merger that preserves branch structure and zero-fills missing
branches. Used by `mkSnapshot.py --haddBatch`.

### `mkPlot.py`, `mkDatacards.py`, `easyDescription.py`
Plot rendering, datacard writing and configuration summaries. Not touched by
the current refactor.

---

## 3. Configuration file conventions

A configuration lives in one directory per analysis and year, e.g.
`plotsConfigurationsRDF/SSWW/2022/`. `configuration.py` is the entry point and
names all the others.

| File | Variable defined | Purpose |
|---|---|---|
| `configuration.py` | module-level names | Entry point: `tag`, `tagger`, `lumi`, output dirs, batch paths, and the filenames below |
| `samples.py` | `samples` | Sample definitions: file lists, weights, per-sample job splitting |
| `cuts.py` | `cuts`, `supercut` | Analysis regions/categories; `supercut` is applied to every sample first |
| `aliases.py` | `aliases` | Derived quantities and C++ helpers injected into the generated code |
| `variables.py` | `variables` | Histogram definitions (binning, expression) |
| `nuisances.py` | `nuisances` | Systematic variations |
| `plot.py` | `plot` | Cosmetics: colours, grouping, legend |
| `structure.py` | `structure` | Datacard process structure |
| `snapshot.py` | — | Snapshot configuration, when used |

Config files are executed with `exec()` by `read_exec` / `load_config_files`,
with `os`, `sys`, `re`, `glob` pre-injected, plus `samples`, `cuts` and
`nuisances` where relevant. Pre-seeded dicts mean a file may write
`cuts['SR'] = '...'` without first creating the dict.

### Per-sample job splitting keys (in `samples.py`)

| Key | Effect |
|---|---|
| `FilesPerJob` | N input files per batch job (default 1) |
| `MaxEventsPerJob` | Upper bound on events per job. It is a cap, not a fixed size: a job may hold fewer. Triggers parallel counting + the `samples_metadata.json` cache |
| `split` | List of percentages, one entry per era, which must sum to 100 |
| `era` | Which era this sample instance uses, as a **1-based** index into `split` |

An era split takes priority over `MaxEventsPerJob`, which takes priority over
`FilesPerJob`. When `split` and `MaxEventsPerJob` are both given, the era is
selected first and `MaxEventsPerJob` then caps the jobs inside it.

#### Splitting one MC production across eras

`split` divides a sample's events into contiguous, non-overlapping slices and
`era` picks one. Sharing a production between three years looks like this, with
the same `split` repeated in each year's `samples.py`:

```python
# SSWW/2024/samples.py        SSWW/2025/samples.py        SSWW/2026/samples.py
'split': [40, 40, 20],        'split': [40, 40, 20],      'split': [40, 40, 20],
'era'  : 1,                   'era'  : 2,                 'era'  : 3,
```

Era 1 takes the first 40% of the events, era 2 the next 40%, era 3 the last
20%. Non-overlap is structural — an era is one contiguous range, so it cannot
be made to overlap another — and the last era always ends on the last event, so
rounding can neither drop nor duplicate events. Files are sorted before slicing
so the partition is reproducible, and `Range()` is applied to the files that
straddle an era boundary.

A malformed `split`/`era` raises `ValueError` rather than silently producing the
wrong events: percentages that do not sum to 100, a missing `era`, a
non-integer `era`, and an `era` outside `1..len(split)` are all rejected.

### Alias keys recognised by the code generator

| Key | Meaning | Where it lands in the generated `.cpp` |
|---|---|---|
| `expr` | RDataFrame expression | `SafeDefine(...)` inside `main()` |
| `samples` | Restrict the alias to these samples | — |
| `external` | Path to a `.cc`/`.cpp`/`.h` file; its contents are spliced in verbatim, **deduplicated by path** | before `main()` |
| `linesToAdd` | Literal lines (typically `#include`), deduplicated | top of file, before `external` content |
| `linesToGlobal` | Literal lines at global scope, deduplicated — **the place for global object instances** | before `main()`, *after* `external` content |
| `linesToProcess` | Literal lines emitted inside `main()`, **not** deduplicated | first statements in `main()` |
| `function` + `variables` | Call a C++ functor over the listed columns | `.Define(...)` chain |
| `xmlfile` + `variables` | TMVA `RReader` model evaluation | `.Define(...)` chain |

Emission order before `main()` is:
`linesToAdd` → `external` → `linesToGlobal`.
That ordering is what makes §7 work: class definitions arrive before the global
instances that use them.

> **Gotcha:** `_build_aliases_ordered` routes an alias into the "function"
> branch only if it has `function`, `linesToAdd` **or** `linesToGlobal`. An alias
> carrying only `external` + `expr` is classified as a plain `expr` and **its
> `external` file is never read.** Always pair `external` with `linesToGlobal`
> (or `linesToAdd`).

---

## 4. Batch workflow

```
mkShapes.py --submitBatch     # generate .cpp/.sh, compile, submit to HTCondor
mkShapes.py --checkBatch      # report which jobs are missing or failed
mkShapes.py --checkBatch --resubmit   # rewrite submit_all.sub with only the missing/failed jobs
mkShapes.py --haddBatch       # submit the merge step as HTCondor jobs
```

Notes:

- `submitBatch` writes one `submit_all.sub` describing all jobs, and runs
  `condor_submit` with `cwd=script_batch_location`.
- Batch script names are derived from the generated `.cpp` filename (not the
  legacy generic `my_script_*`), so jobs are identifiable in the queue.
- **The generation step is not a separate invocation.** `--submitBatch`,
  `--localCompile` and `--parallelCompile` all run it first, so
  `mkShapes.py --submitBatch` on an empty job directory works. Only
  `--checkBatch` and the hadd steps skip it, because they read what is already
  on disk. The `.cpp` are written through `_write_if_changed`, so regenerating
  is cheap and leaves an unchanged file's mtime alone.
- The two compilation modes are chosen **per job**: a job ships its binary if
  that binary exists *and is newer than its `.cpp`*, otherwise it ships the
  `.cpp` and compiles on the node (a minute or two). So nothing has to be
  compiled before submitting, a partially built set is fine, and a binary left
  over from an older configuration is recompiled instead of being shipped
  silently. `submitBatch` prints the tally and aborts only when a job has
  neither, which means its scripts were never generated.
  `transfer_input_files` is one line for the whole cluster, so when any job
  compiles on the node the list carries what that needs.
- EOS copies use `xrdcp -f`; the job creates its output directory with
  `mkdir -p` before the executable runs.
- **The X509 proxy is optional.** `submitBatch` looks for one
  (`$X509_USER_PROXY`, `~/.proxy`, `/tmp/x509up_u<uid>`), ships it with the job
  and warns if there is none, but never aborts: CERN EOS over
  `root://eosuser.cern.ch` / `root://eoscms.cern.ch` authenticates with
  Kerberos too, so a user with EOS access and no grid certificate can submit.
  A proxy *is* needed to read inputs through a CMS redirector
  (`root://xrootd-cms...`, `cms-xrd-global`), which only accepts GSI — without
  one those jobs fail on the node with an xrootd authentication error. Relying
  on Kerberos also means the job needs the credential forwarded
  (`MY.SendCredential = True` in the submit file, with a credential stored on
  the schedd); that is not set by default.
- `--haddBatch` uses `haddfast.py` for shapes and `haddnano.py` for snapshots.
- For shapes, both `--hadd` and `--haddBatch` write `root_file_joined.root`
  into `outputDir`. That name is not free: `mkPlot.py` takes it as the default
  `--inputFile`, `mkDatacards.py` hardcodes it, and `mkMergeYears.py` reads
  `<folder>/rootFile/<tag>/root_file_joined.root` for each year. Changing it
  means changing all four. The per-job inputs are named
  `root_file___<sample>_<subname>_<i>.root`, which keeps them distinct from the
  merged file so re-running cannot fold the merged output back into itself.
- **The merge inputs come from `samples.py`, not from the output directory.**
  `selectOutputFiles()` matches the files on disk against the
  `(sample, subname)` pairs the configuration declares; `checkBatch` and the
  pre-hadd completeness check walk the declared subnames too, not
  `os.listdir`. Anything left over is listed and skipped, never merged. See
  "Stale per-job files" below for why this is not a convenience.
- For snapshots the merged file is `snapshot_<tag>.root`, deliberately *not*
  matching `nanoLatino_snapshot____*.root`, since `makeFromSnapshot` globs that
  pattern per sample and would otherwise read the merged file alongside its
  own parts.
- Compilation links ROOT (`root-config --cflags/--libs`) plus
  `-lcorrectionlib`, so `correction.h`-dependent helpers compile and link
  without extra configuration.

### Stale per-job files

A sample's output files are keyed by the **sample name plus the dataset
subname**, but the histograms *inside* them are keyed by the sample name alone.
So after swapping a sample's dataset — say `WZJJ_QCD` from the madgraph
`WZJJ_QCD` to the POWHEG `WZTo3LNu` — the directory holds

```
root_file___WZJJ_QCD_WZJJ_QCD_0.root     <- previous configuration
root_file___WZJJ_QCD_WZTo3LNu_0.root     <- current one
```

and both contain `WZ_SR/histo_WZJJ_QCD`. A merge driven by a glob therefore
**adds the old sample to the new one** instead of replacing it: same folder,
same histogram name, so to `hadd` they are two chunks of one sample. No error,
no warning. Measured on a one-bin test: 100 + 42 came out as 142.

The same applies with no dataset swap at all. Commenting one
`nanoGetSampleFiles` line out of a multi-dataset sample, or commenting out a
whole sample, leaves every file that sample ever produced behind, and a glob
puts them all back in. That is why the configuration is the authority.

Two consequences worth knowing:

- The leftover files are *reported*, not deleted. Deleting them is safe and
  recommended, but nothing in the framework does it for you — an interrupted
  reconfiguration should not destroy results.
- Dropping a subname also leaves its job scripts under
  `script_batch_location/<sample>/<subname>/`. Those are ignored as well, so a
  deleted output file does not make `--checkBatch` or the pre-hadd check report
  phantom missing jobs and refuse to merge.

## 5. Snapshot workflow

```
mkSnapshot.py --submitBatch           # produce skimmed trees
mkSnapshot.py --checkBatch [--resubmit]
mkSnapshot.py --haddBatch             # merge with haddnano.py
mkSnapshot.py --patch                 # top up an incomplete snapshot set
mkShapes.py --fromSnapshot --submitBatch   # fill histograms from snapshots
```

`mkShapes.py --fromSnapshot` routes through `makeFromSnapshot()`, which reads
the snapshot trees instead of nanoAOD. `SnapshotMetadata` records what exists so
`--patch` can regenerate only what is missing.

### Shape nuisances in a snapshot

**There is no flag.** A `suffix` nuisance left uncommented in `nuisances.py` *is*
the request, so the varied trees are part of every step by themselves:

```
mkSnapshot.py --submitBatch           # generates, then queues nominal + varied
mkSnapshot.py --checkBatch [--resubmit]
mkSnapshot.py --haddBatch             # nominal only, by design (see below)

mkSnapshot.py --localCompile          # optional: build locally first
mkSnapshot.py --noShapeNuisances ...  # escape hatch: nominal only
```

`--submitBatch` is enough on its own: it generates the `.cpp` (nominal and
varied), finds no binaries and ships the sources, and the nodes compile.
`--localCompile` only trades node time for local time.

`--shapeNuisances` is still accepted and ignored, with a note. It used to be
required on every step, which made it possible to generate the varied `.cpp`
and then submit only the nominal jobs — a half-run with nothing to say so.
`makeScripts()` now writes both makefiles, `--localCompile` builds both, and
`SnapshotFactory.doShapeNuisances()` is the single predicate.

Only the `suffix` / `suffixFly` nuisances need any of this — they read a varied
copy of the nanoAOD, so they need their own pass. The `weight` ones do not:
`_build_full_variables_to_dump` already adds a `my_sample_weight_<name>Up` and
`Down` branch for each of them, so the varied weight is in the nominal tree. At
148 nuisances that is 258 extra branches, 206 of them PDF, and it dominates the
file size.

Naming, and why it is in one place. A varied job is
`my_snapshot_<nuisance><Up|Down>_<sample>_<subname>_<i>` and writes
`nanoLatino_snapshot____<sample>_<subname>_<i>_<nuisance><Up|Down>.root` — the
nuisance leads in the job name and trails in the file name. `checkBatch` finds
jobs by listing the `.sh` on disk and strips `job_prefix()_`, so it has to undo
that reordering; `SnapshotFactory._nuisance_names()` and
`_output_file_from_suffix()` are the only two places that know the layout.

The varied files are deliberately *not* merged by `--haddBatch`:
`selectOutputFiles` requires what follows `<sample>_<subname>_` to be the
integer job index alone, so `..._0_JERUp.root` is reported as an orphan and
skipped. They are inputs for `--fromSnapshot`, not chunks of the nominal tree.

---

## 6. Year merging (`mkMergeYears.py`)

Sums the per-year shape files into one multi-year file. The sum is **over years
only**: for every cut/variable folder and every sample,
`histo_DY = histo_DY(2022) + histo_DY(2022EE) + ...`, and likewise for every
*shape* nuisance variation. Nothing is summed across samples, regions or
variables. `lnN` nuisances write no histogram; the one case where a `lnN` enters
the file is a year that lacks a shape variation another year has, where its own
nominal scaled by its `lnN` value stands in. `DATA` gets no variations.

```
mkMergeYears.py --pycfg merge_configuration.py     # everything locally
```

### Output naming

The output is `root_file_joined.root`, the same name `mkShapes` produces per
year, in `outputDir`. **An existing file is never overwritten**: the first merge
writes the base name, every later one is suffixed `_v1`, `_v2`, … so the suffix
lands on the *new* file. `versionedOutputPath()` owns this.

> **Gotcha:** `mkPlot.py` defaults `--inputFile` to `root_file_joined.root` and
> `mkDatacards.py` hardcodes it (`mkDatacards.py:909`). From the second merge on
> they still read the *first* file, so pass the versioned path explicitly.

### Batch splitting

A folder of one cut never depends on another cut, so the work splits by cut.
Each job writes `root_file_joined___<cut>.root`; the partials hold **disjoint**
folders, so the hadd is a union, not an addition, and cannot double count.

```
mkMergeYears.py --submitBatch     # one condor job per cut
mkMergeYears.py --checkBatch      # which partials are missing
mkMergeYears.py --hadd            # union into root_file_joined.root (haddfast -C)
```

- `--cut` takes a **comma separated list, all of which run in one process** on
  one node, sequentially, into one output file. The list makes jobs *fewer and
  bigger*, not more parallel. `--submitBatch` queues one cut per job.
- `--output` sets an explicit file name and bypasses the `_vN` versioning. The
  generated job is exactly one `--cut` plus one `--output`.
- A cut with `categories` contributes its unrolled sub-regions, so the job count
  is the unrolled count, not `len(cuts)`.

> **Do not scale this to hundreds of jobs as it stands.** Every job copies all
> of each year's `root_file_joined.root`, whole, even though it reads one region.
> Ten jobs over seven years is 70 copies, which is free. Splitting further —
> region × variable would be 380 jobs and 2660 copies — makes EOS bandwidth the
> bottleneck and the parallelism stops paying. Before going there, change the
> reading: either stream over xrootd with no local copy, or have each job fetch
> only the folders it needs. There is no `--variable` option today; adding one
> is small, but it is the copying that has to be fixed first.

### Why it is shaped this way — measured, not assumed

Benchmarked on 7 synthetic years, 10 cuts × 8 variables × 25 samples × 2 shape
nuisances = 9680 output histograms, and a full-scale variant with 38 variables
and 10 nuisances = 191900 output histograms per year.

| | |
|---|---|
| load-everything version | 352 s |
| one folder at a time | 157 s, output identical bin by bin |
| peak RSS | 392 MB → 401 MB, i.e. **unchanged** |

The merge reads one cut/variable folder from every year, sums, writes and frees,
rather than loading every year's whole file first. Two things to know:

1. **The speed came from `SetDirectory(0)`, not from the streaming itself.** The
   old code called `rootFileNew.cd()` and then `Clone()`d every histogram of
   every year, so all of them were appended to one directory's in-memory object
   list, which makes the cost quadratic. Measured on 40k clones: 4.7 s attached
   with distinct names, **70 s attached with repeating names** — which is this
   case, ~500 names reused in every folder of every year — against 0.17 s
   detached.
2. **Memory was not the win.** The peak is ROOT plus the open input files, not
   the histograms. What streaming removes is the *growth* with the number of
   variables and shape nuisances, which the old shape could not survive.

The remaining floor is per-histogram object lifetime: creating and freeing one
histogram read from a file costs ~0.17 ms, and an explicit `Delete()` and simply
dropping the reference measured **the same**, so there is no trick there. The
merge pays that once per histogram per year, which is why it splits by cut
instead of being optimised further in Python. For reference, the useful work
alone — `Clone` + 6 `Add` + `Write` — is 0.026 ms, and reading all of one year's
histograms is 2.9 s.

### uproot — tried, measured, deliberately not adopted

| | per output histogram |
|---|---|
| uproot read 7 years + sum | 2.33 ms |
| uproot write, with `sumw2` | 0.32 ms |
| **uproot total** | **~2.65 ms** |
| current PyROOT | ~10.2 ms |

So uproot is ~3.8× faster, and `sumw2` survives correctly (written via
`uproot.writing.identify.to_TH1x`, read back with the right errors). It was
**not** adopted because the per-cut split already gives 10× and is verified
identical, while uproot means rewriting `merge()` with different I/O semantics.
The two compose, so uproot is a possible second step, not the first.

If it is ever picked up:

- Use `uproot.recreate()` only — write once, never update or delete in place.
  uproot's update path is where the surprises are.
- `fTsumwx` and `fTsumwx2` carry the mean and rms. The benchmark set them to
  zero; real code must read them from the inputs and sum them, or `GetMean()` on
  the merged file is wrong even though every bin is right.
- Decide what happens to objects that are not `TH1`. uproot will not pass them
  through on its own.

### Duplicate key cycles

`haddfast.py` merges by calling `TH1::Write()` without `kOverwrite`, so a
histogram that is *added* to an existing one gains a key cycle instead of
replacing it: `;1` holds the first chunk alone and the highest cycle the running
sum. `TFile::Get()` resolves a bare name to the highest cycle, so `mkPlot` and
`mkDatacards` were always fine, but **`GetListOfKeys()` returns every cycle, in
descending order**. Code that keys a dict by name and lets the last iteration
win therefore keeps the *partial* histogram. That was a real bug in this script
and it under-counted every sample. Both ends are fixed: `mkShapes.py --haddBatch`
passes `-C` so new files carry one cycle, and `readFolder()` dedups on the
highest cycle like `haddfast.merge()` already did.

---

## 7. gInterpreter policy — **forbidden**

Runtime JIT declaration has been **rejected by reviewers** and must not be
reintroduced.

### Why it was there in the first place

This is the key point, and it is not obvious. Two different `Define` overloads
are in play:

1. **String-expression `Define`** — `SafeDefine(df, name, "some expression")`
   resolves to `df.Define(name, std::string expr)` (`mkShapes.py:135-138`).
   RDataFrame hands that string to **Cling at runtime**. Cling only sees what
   has been declared *into the interpreter*; it does **not** see types or
   objects defined in the compiled translation unit of the generated binary.
2. **Typed callable `Define`** — `df.Define(name, callable, {"col", ...})`.
   Fully compiled, no interpreter involved.

The fake-rate and b-tag aliases used form 1 (`'expr'`). For a JIT-compiled
expression string to resolve `fr_reader_nominal(...)`, that class *and* that
instance had to exist **inside Cling** — which is exactly what
`gInterpreter->Declare` / `ProcessLine` did. So the original pattern was not
gratuitous: with `'expr'`, dropping `Declare` genuinely breaks.

It follows that **`external` + `linesToGlobal` alone is not a fix.** A class
spliced before `main()` is invisible to the JIT that compiles the `expr`
string. The alias has to stop using `expr` at all.

### Required instead: the VBS_functions form

Use the typed-callable form already established by `VBS_functions.cpp`:
`function` + `variables` + `external`, and **no `expr`**.

```python
aliases['zstar_max'] = {
    'function' : 'computeZstar',
    'variables': ['Lepton_eta', 'CleanJet_eta'],
    'external' : '../utils/macros/VBS_functions.cpp',
    'samples'  : ALL,
}
```

For a class that needs construction arguments, add the instance at global scope
with `linesToGlobal` and name that instance as the `function`:

```python
aliases['fakeW'] = {
    'external'     : f'{macros}fake_rate_reader_class.cc',
    'linesToGlobal': [
        f'fake_rate_reader fr_reader_nominal("{eleWP}", "{muWP}", "nominal", 2, "std", "{fakerates}", "2022_v12_pt");'
    ],
    'function'     : 'fr_reader_nominal',
    'variables'    : ['Lepton_pdgId', 'Lepton_pt', 'Lepton_eta',
                      'Lepton_isTightMuon', 'Lepton_isTightElectron',
                      'Lepton_muonIdx', 'CleanJet_pt', 'nCleanJet'],
    'samples'      : ['Fake'],
}
```

Emission order before `main()` (§3) is `linesToAdd` → `external` →
`linesToGlobal`, so the class is always defined before the instance.

### Three constraints the typed `Define` imposes

These bit once each during the conversion. All three are properties of
`RDataFrame::Define(name, callable, columns)`, not of this framework.

**1. No literals — configuration goes in the constructor.** `variables`
entries become columns, so a configuration constant cannot be an `operator()`
argument. `btagSFbc` / `btagSFlight` took 6 columns plus 4 strings (`WP`,
`shift`, `tagger`, `year`); those four moved into the constructor and became
members. `fake_rate_reader` already took only columns. A useful side effect:
`cset->at(tagger + ...)` now resolves once per job instead of once per event.

**2. The callable must be copyable.** `Define` takes it **by value**. A
`std::unique_ptr` member deletes the copy constructor, and a user-declared
destructor suppresses the implicit move constructor — with both, the class is
neither copyable nor movable and `Define` will not compile. Use `shared_ptr`
for owned state and let the destructor go, so copies share it. Check with
`static_assert(std::is_copy_constructible_v<T>)`.

**3. Column types must match exactly.** This is the one that costs the most
time. A string `expr` is JIT-compiled by Cling, which converts silently; the
typed `Define` compares types and throws:

```
RDataFrame: type mismatch: column "_var_btagSFbc_3" is being used as
unsigned int but the Define or Vary node advertises it as int
```

Branch types in these nanoAODs are not what the names suggest, so **pin them
in `variables`** rather than in the C++ — the branch type belongs to the input
data and varies by nanoAOD version, while the helper signature is shared by
every year:

| Branch | Actual type | Note |
|---|---|---|
| `nCleanJet` | `Int_t` | *not* `UInt_t`; cast for a helper taking `unsigned int` |
| `Jet_hadronFlavour` | `RVec<unsigned char>` | `UChar_t` from NanoAODv12 on |
| `Lepton_isTightMuon_<WP>`, `Lepton_isTightElectron_<WP>` | `RVec<bool>` | cast to `RVec<int>` |
| `CleanJet_jetIdx` | `RVec<ULong64_t>` | genuinely unsigned 64-bit |
| `Lepton_muonIdx` | `RVec<int>` | cast kept anyway: only verified on one year |
| `*_pt`, `*_eta`, `Jet_btag*` | `RVec<float>` | |
| `Lepton_pdgId` | `RVec<int>` | |

`(unsigned int)nCleanJet` and `ROOT::VecOps::RVec<int>(branch)` both work as
`variables` expressions; the latter accepts any convertible element type.

To check a file's types before submitting, rather than finding out from a
crashed job:

```python
import ROOT
t = ROOT.TFile.Open(path).Get("Events")
b = t.GetBranch("Jet_hadronFlavour")
print(b.GetLeaf("Jet_hadronFlavour").GetTypeName())
```

### Forbidden patterns

```python
'linesToAdd'    : ['#include "TInterpreter.h"'],
'linesToProcess': [f'gInterpreter->Declare(R"FAKE({macro_source})FAKE");'],
```

```python
'linesToProcess': [f'gInterpreter->ProcessLine("btagSFbc inst(\\"{path}\\", \\"{year}\\");");'],
```

Also forbidden: `gInterpreter->Declare`, `gInterpreter->ProcessLine`,
`gROOT->ProcessLine` and `TInterpreter.h` anywhere in generated C++ or in
`aliases.py` / `nuisances.py`.

Rationale beyond the reviewers' objection: JIT-ing a class at runtime inside an
already-compiled binary costs Cling startup on every batch job, defers all type
errors to runtime, and hides the source from the compiler.

### Notes

- `external` deduplicates by path, so every alias may safely declare its own
  dependency.
- Helper sources must be include-guard-safe and self-contained. The macros in
  `plotsConfigurationsRDF/SSWW/utils/macros/` already are:
  `fake_rate_reader_class.cc` (`FAKERATEREADER`), `evaluate_btagSFbc.cc`
  (`BTAGSF`), `evaluate_btagSFlight.cc` (`BTAGLIGHTSF`), `VBS_functions.cpp`.
- The separate `_fakeRateDecl` / `_btagSF{bc,light}Decl` aliases carrying
  `'expr': '1.'` existed only to host the one-off `Declare`. They are no longer
  needed and should be deleted.

## 8. Commit conventions

Conventional-commit prefixes, imperative mood, lower-case subject:

```
feat:     new user-visible capability
fix:      bug fix or correction of a rejected pattern
refactor: behaviour-preserving restructuring
docs:     documentation only
chore:    tooling, cleanup, no behaviour change
```

Rules for this repo:

- **One feature per commit**, self-contained and reviewable on its own.
- **Every commit must leave the code in a working state** — the branch is
  reviewed commit by commit upstream.
- Framework changes (`scripts/`, `lib/`) go to **latinoRDF**; analysis
  configuration changes go to **plotsConfigurationsRDF**. Never mix the two in
  one commit — they are different repositories.
- Feature work happens on a branch (e.g. `feature/clean-history`), never
  directly on `main`.
