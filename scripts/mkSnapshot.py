#!/usr/bin/env python3

import sys
argv = sys.argv
sys.argv = argv[:1]

import os
import re
import json
import glob
import subprocess
import ROOT

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysisFactory import (AnalysisFactory, defaultParser, loadDefaultOptions,
                             read_exec, load_config_files, _base_ns, _SCRIPTS_DIR)


# -----------------------------------------------------------------------
#  Helpers
# -----------------------------------------------------------------------

def _extract_branch_names(expression):
    """Extract likely branch/variable names from a C++ expression string."""
    # remove strings, numbers, operators — keep identifiers
    tokens = re.findall(r'\b[A-Za-z_][A-Za-z0-9_]*\b', expression)
    # filter out C++ keywords and ROOT types
    _skip = {'auto','int','float','double','bool','true','false','if','else',
              'for','while','return','const','void','std','ROOT','RVec',
              'abs','sqrt','pow','sin','cos','log','exp','min','max'}
    return [t for t in tokens if t not in _skip]


def _collect_cut_variables(cuts):
    """Collect all variable names referenced in cut expressions."""
    vars_used = set()
    for cutName, cut in cuts.items():
        expr = cut['expr'] if isinstance(cut, dict) else cut
        for v in _extract_branch_names(expr):
            vars_used.add(v)
    return vars_used


def _collect_alias_dependencies(aliases_ordered):
    """Collect all input variables needed by aliases."""
    deps = {}
    for aliasName, (kind, payload) in aliases_ordered.items():
        if kind == 'expr':
            deps[aliasName] = _extract_branch_names(payload)
        elif kind == 'function' and 'variables' in payload:
            deps[aliasName] = list(payload['variables'])
    return deps


# -----------------------------------------------------------------------
#  SnapshotMetadata
# -----------------------------------------------------------------------

class SnapshotMetadata:
    """
    Tracks what has been saved in the snapshot and what still needs to be done.
    Saved as snapshot_metadata.json in folder_where_to_save_trees.
    """

    def __init__(self, path):
        self._path = path
        self._data = {
            'cut_to_dump'        : '',
            'cut_expression'     : '',
            'cut_variables'      : [],
            'variables_saved'    : [],
            'alias_dependencies' : {},
            'weight_nuisances'   : [],
            'shape_nuisances'    : [],
            'patches'            : [],
            'completed_files'    : {},
        }
        if os.path.exists(path):
            with open(path, 'r') as f:
                self._data.update(json.load(f))

    def save(self):
        os.makedirs(os.path.dirname(self._path), exist_ok=True)
        with open(self._path, 'w') as f:
            json.dump(self._data, f, indent=2)

    def get(self, key, default=None):
        return self._data.get(key, default)

    def set(self, key, value):
        self._data[key] = value

    def mark_completed(self, job_key, path):
        try:
            mtime = os.path.getmtime(path)
        except Exception:
            mtime = 0
        self._data['completed_files'][job_key] = {'path': path, 'mtime': mtime}

    def is_completed(self, job_key, path):
        entry = self._data['completed_files'].get(job_key)
        if not entry:
            return False
        if not os.path.exists(path):
            return False
        try:
            mtime = os.path.getmtime(path)
            return abs(mtime - entry.get('mtime', -1)) < 2.0
        except Exception:
            return False

    def missing_variables(self, requested):
        saved = set(self._data.get('variables_saved', []))
        return [v for v in requested if v not in saved]

    def can_patch_from_snapshot(self, missing_vars):
        """Check if missing variables can be derived from already-saved branches."""
        saved = set(self._data.get('variables_saved', []))
        alias_deps = self._data.get('alias_dependencies', {})
        for var in missing_vars:
            deps = alias_deps.get(var, [var])
            if not all(d in saved for d in deps):
                return False
        return True


# -----------------------------------------------------------------------
#  SnapshotFactory
# -----------------------------------------------------------------------

class SnapshotFactory(AnalysisFactory):
    """Produces ROOT files with a TTree (Events) snapshot per sample/chunk."""

    def __init__(self):
        super().__init__()
        self._variables_to_dump          = []
        self._cut_to_dump                = ''
        self._folder_where_to_save_trees = ''
        self._scripts_run_folder         = 'scripts_snapshot_run'
        self._script_batch_location      = './scripts_snapshot_batch/'
        self._metadata                   = None
        self._skipShapeNuisances         = False
        # set only while generating a varied job's .cpp (see makeShapeNuisances)
        self._suffix_variation           = None

    # _____________________________________________________________________________
    def setSkipShapeNuisances(self, skip):
        """Opt out of the shape-nuisance jobs even though nuisances.py declares them."""
        self._skipShapeNuisances = bool(skip)

    # _____________________________________________________________________________
    def doShapeNuisances(self):
        """Whether the varied trees are part of this run.

        Derived from the configuration, not from a command-line flag: a
        'suffix' nuisance left uncommented in nuisances.py IS the request. The
        flag only existed because submitBatch could not see those jobs, and
        having to repeat it on every step was a way to generate the varied
        .cpp and then submit only the nominal ones -- a silent half-run.

        setSkipShapeNuisances() is the escape hatch, for a quick nominal-only
        snapshot without editing nuisances.py.
        """
        return bool(self._suffixNuisances()) and not self._skipShapeNuisances

    # _____________________________________________________________________________
    def _suffixNuisances(self):
        """The shape nuisances that need a separate pass over a varied tree."""
        return [(n, nu) for n, nu in self._nuisances.items()
                if nu.get('type') == 'shape'
                and nu.get('kind') in ('suffix', 'suffixFly')]

    def _nuisance_names(self, nuisanceName, direction, sampleName, subname, i):
        """(job_name, out_file) of one varied job -- the ONLY place they are spelled.

        makeShapeNuisances, _batch_jobs_for_chunk and _output_file_from_suffix
        all go through this, because the job name and the output name order
        their fields differently (nuisance first in one, last in the other) and
        hand-writing that twice is how checkBatch came to look for a file no
        job ever wrote.
        """
        return (f"{self.job_prefix()}_{nuisanceName}{direction}_{sampleName}_{subname}_{i}",
                f"nanoLatino_snapshot____{sampleName}_{subname}_{i}_{nuisanceName}{direction}.root")

    # _____________________________________________________________________________
    def _batch_jobs_for_chunk(self, sampleName, subname, i):
        jobs = super()._batch_jobs_for_chunk(sampleName, subname, i)
        if not self.doShapeNuisances():
            return jobs
        for nuisanceName, nuisance in self._suffixNuisances():
            if sampleName not in nuisance.get('samples', {}):
                continue
            for direction in ('Up', 'Down'):
                jobs.append(self._nuisance_names(nuisanceName, direction,
                                                 sampleName, subname, i))
        return jobs

    # _____________________________________________________________________________
    def setSnapshotValues(self, variables_to_dump, cut_to_dump, folder_where_to_save_trees):
        self._variables_to_dump          = variables_to_dump
        self._cut_to_dump                = cut_to_dump
        self._folder_where_to_save_trees = os.path.abspath(
            os.path.expandvars(folder_where_to_save_trees))
        meta_path = os.path.join(self._folder_where_to_save_trees, 'snapshot_metadata.json')
        self._metadata = SnapshotMetadata(meta_path)

    # _____________________________________________________________________________
    def job_prefix(self):
        return "my_snapshot"

    def output_file_name(self, sampleName, subname, i):
        return f"nanoLatino_snapshot____{sampleName}_{subname}_{i}.root"

    def _output_file_from_suffix(self, sampleName, subname, suffix):
        """The output file of the job whose name is job_prefix() + '_' + suffix.

        checkBatch strips job_prefix() and hands us the rest, so for a nominal
        job the suffix IS the output stem. For a varied one it is not: the job
        name starts with <nuisance><Up|Down>_ while the file ends with it, so
        the field has to be moved back. Without this every shape-nuisance job
        was reported missing even after it had run.
        """
        for nuisanceName, _ in self._suffixNuisances():
            for direction in ('Up', 'Down'):
                head = f"{nuisanceName}{direction}_"
                if suffix.startswith(head):
                    stem = suffix[len(head):]
                    return os.path.join(
                        self._folder_where_to_save_trees,
                        f"nanoLatino_snapshot____{stem}_{nuisanceName}{direction}.root")
        return os.path.join(self._folder_where_to_save_trees,
                            f"nanoLatino_snapshot____{suffix}.root")

    # _____________________________________________________________________________
    #  Variables to dump — with auto-collection of cut vars, alias vars, weight nuisances
    # _____________________________________________________________________________

    def _build_full_variables_to_dump(self, aliases_ordered, nuisances):
        """
        Build the complete list of variables to save in the snapshot:
        - explicit variables_to_dump from snapshot.py
        - my_sample_weight (always)
        - weight nuisance branches (my_sample_weight_<name>Up/Down)
        - cut variables (branches used in cut expressions)
        - simple expr aliases (auto-added)
        """
        vtd = list(self._variables_to_dump)

        # always include my_sample_weight
        if 'my_sample_weight' not in vtd:
            vtd = ['my_sample_weight'] + vtd

        # weight nuisances as extra branches
        weight_nuisance_branches = []
        for nuisanceName, nuisance in nuisances.items():
            if nuisance.get('type') == 'shape' and nuisance.get('kind') == 'weight':
                wUp   = f'my_sample_weight_{nuisanceName}Up'
                wDown = f'my_sample_weight_{nuisanceName}Down'
                if wUp   not in vtd: vtd.append(wUp);   weight_nuisance_branches.append(wUp)
                if wDown not in vtd: vtd.append(wDown);  weight_nuisance_branches.append(wDown)

        # cut variables
        cut_vars = _collect_cut_variables(self._cuts)
        for v in cut_vars:
            if v not in vtd:
                vtd.append(v)

        # simple expr aliases
        for aliasName, (kind, payload) in aliases_ordered.items():
            if kind == 'expr' and aliasName not in vtd:
                vtd.append(aliasName)

        return vtd, weight_nuisance_branches

    # _____________________________________________________________________________
    #  C++ weight nuisance defines
    # _____________________________________________________________________________

    def _build_weight_nuisance_defines(self, sampleName, weight, nuisances):
        """Generate SafeDefine calls for weight nuisance branches."""
        out = ""
        for nuisanceName, nuisance in nuisances.items():
            if (nuisance.get('type') == 'shape'
                    and nuisance.get('kind') == 'weight'
                    and sampleName in nuisance.get('samples', {})):
                samples_entry = nuisance['samples'][sampleName]
                if len(samples_entry) == 3:
                    if samples_entry[0] in weight:
                        var_up = weight.replace(samples_entry[0], samples_entry[1])
                        var_do = weight.replace(samples_entry[0], samples_entry[2])
                    else:
                        continue
                else:
                    var_up = f"({weight}) * {samples_entry[0]}"
                    var_do = f"({weight}) * {samples_entry[1]}"
                out += f'    varied_df = SafeDefine(varied_df, "my_sample_weight_{nuisanceName}Up",   "{var_up}");\n'
                out += f'    varied_df = SafeDefine(varied_df, "my_sample_weight_{nuisanceName}Down", "{var_do}");\n'
        return out

    # _____________________________________________________________________________
    #  Varied input — the nominal tree, with the suffix folder as a friend
    # _____________________________________________________________________________

    def _build_cpp_input_files_extra(self, tree_name, file_path, sampleName):
        """Attach the variation's files as a friend of the nominal chain.

        A '*_suffix' production is not a full nanoAOD: it holds only the
        branches the variation touches, each named '<nominal branch>_<map>'.
        Read on its own it defines almost nothing, so the first alias or cut
        that mentions an ordinary branch dies in the JIT with 'use of
        undeclared identifier Lepton_pt'. The varied job therefore reads the
        nominal files and adds the varied ones as a friend, the way mkShapes
        already does, and _build_cpp_after_df() swaps the columns over.
        """
        var = self._suffix_variation
        if not var:
            return ""
        out  = "//  Variation input files (friend of the nominal chain)\n"
        out += f'    auto* varied_friend = new TChain("{tree_name}");\n'
        for file_name in var['files']:
            out += f'    varied_friend->Add("{file_name}");\n'
        out += f'    nominal->AddFriend(varied_friend, "{var["map"]}");\n'
        out += '    auto nominalBranches = getBranchNames(nominal);\n'
        out += '    auto variedBranches  = getBranchNames(varied_friend);\n'
        out += '    \n'
        return out

    def _build_cpp_after_df(self, sampleName):
        """Redefine every nominal column the variation provides.

        Which columns those are is a property of the input files, not of the
        configuration, so the list is walked on the node at runtime: every
        friend branch ending in '_<map>' whose nominal counterpart exists is
        redefined from it. Everything downstream -- aliases, variables, the
        cut -- then sees varied values under the nominal names, which is what
        makes the varied snapshot a drop-in friend for mkShapes
        --fromSnapshot: it reads '<nuisance><Up|Down>.<nominal column>'.
        """
        var = self._suffix_variation
        if not var:
            return ""
        out  = "\n//  Replace the nominal columns with their varied counterparts\n"
        out += f'    const std::string varied_suffix = "_{var["map"]}";\n'
        out += '    for (const auto& branch : variedBranches) {\n'
        out += '      if (branch.size() <= varied_suffix.size()) continue;\n'
        out += '      if (branch.compare(branch.size() - varied_suffix.size(),\n'
        out += '                         varied_suffix.size(), varied_suffix) != 0) continue;\n'
        out += '      const std::string nomCol = branch.substr(0, branch.size() - varied_suffix.size());\n'
        out += '      if (std::find(nominalBranches.begin(), nominalBranches.end(), nomCol)\n'
        out += '          == nominalBranches.end()) continue;\n'
        out += '      varied_df = varied_df.Redefine(nomCol, branch);\n'
        out += '    }\n'
        out += '    \n'
        return out

    # _____________________________________________________________________________
    #  Cuts + snapshot — the mkSnapshot-specific output fragment
    # _____________________________________________________________________________

    def _unroll_cuts(self):
        temp_cuts = {}
        for cutName, cut in self._cuts.items():
            if isinstance(cut, dict):
                expression = cut['expr']
                categories = cut.get('categories', {})
                if categories:
                    temp_cuts[cutName] = "(" + expression + ")"
                    for cat_name, cat in categories.items():
                        temp_cuts[cutName + "_" + cat_name] = "(" + expression + ") && (" + cat + ")"
                else:
                    temp_cuts[cutName] = expression
            else:
                temp_cuts[cutName] = cut
        self._cuts = temp_cuts

    def _build_cpp_cuts_and_output(self, sampleName):
        self._unroll_cuts()
        out = ""
        if self._cut_to_dump not in self._cuts:
            print(f"  WARNING: cut_to_dump='{self._cut_to_dump}' not found in cuts.")
            return out

        this_cut = self._cuts[self._cut_to_dump]
        out += f'\n    auto node_{self._cut_to_dump} = current_node.Filter("{this_cut}", "{self._cut_to_dump}");\n\n'

        columns_str = ", ".join([f'"{v}"' for v in self._current_variables_to_dump])
        out += f'    std::vector<std::string> requested = {{{columns_str}}};\n'
        out += f'    std::vector<std::string> existing = node_{self._cut_to_dump}.GetColumnNames();\n'
        out += f'    std::vector<std::string> columns_to_save;\n'
        out += f'    for (const auto& col : requested) {{\n'
        out += f'      if (std::find(existing.begin(), existing.end(), col) != existing.end())\n'
        out += f'        columns_to_save.push_back(col);\n'
        out += f'    }}\n\n'
        return out

    # _____________________________________________________________________________
    #  create_cpp_source_list_of_files
    # _____________________________________________________________________________

    def create_cpp_source_list_of_files(self, output_name, tree_name, file_path,
                                         output_root_file_name, sampleName, weight,
                                         aliases_to_be_defined, functions_to_be_defined,
                                         aliases_ordered=None,
                                         range_start=None, range_end=None,
                                         extra_weight_defines=""):
        frags = self._build_cpp_common(tree_name, file_path, sampleName, weight,
                                       aliases_ordered or {}, range_start, range_end)

        # A Range() is single-thread by construction. A varied job is run
        # single-thread too: its input is the nominal chain with a friend
        # chain attached, and entry-by-entry alignment between the two is
        # what makes the variation correct -- not worth risking for threads
        # the batch node does not give us anyway (condor allots 1 CPU).
        single_thread = ((range_start is not None and range_end is not None)
                         or self._suffix_variation is not None)
        mt_setting = "ROOT::DisableImplicitMT();" if single_thread else "ROOT::EnableImplicitMT();"

        snapshot_call = f'node_{self._cut_to_dump}.Snapshot("Events", "{output_root_file_name}", columns_to_save);\n'

        cpp_code = f"""
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RVec.hxx"
#include "TFile.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include <algorithm>
#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <typeinfo>
#include "TMVA/RReader.hxx"
#include "TMVA/RInferenceUtils.hxx"

ROOT::RDF::RNode SafeDefine(ROOT::RDF::RNode df, std::string name, std::string expr) {{
    auto colNames = df.GetColumnNames();
    if (std::find(colNames.begin(), colNames.end(), name) == colNames.end())
        return df.Define(name, expr);
    return df;
}}

std::vector<std::string> getBranchNames(TTree* tree) {{
    std::vector<std::string> names;
    TObjArray* branches = tree->GetListOfBranches();
    if (!branches) return names;
    for (int i = 0; i < branches->GetEntries(); ++i)
        names.push_back(branches->At(i)->GetName());
    return names;
}}

{frags['lines_to_add']}

{frags['code_of_function_to_include']}

{frags['lines_to_global']}

int main() {{
    {mt_setting}

    {frags['lines_to_process']}

    {frags['define_input_files_logic']}

    {frags['define_aliases']}

    {frags['define_weights']}

    // --- weight nuisance branches ---
    {extra_weight_defines}

    {frags['booking_logic']}

    {frags['define_variables_logic']}

    {snapshot_call}

    return 0;
}}
"""
        self._write_if_changed(f"{output_name}.cpp", cpp_code)

    # _____________________________________________________________________________
    #  generate_makefile (standalone — no liblibrary_utils)
    # _____________________________________________________________________________

    def generate_makefile(self, file_paths, makefile_name="MakefileSnapshot"):
        cpp_flags = "$(shell root-config --cflags)"
        # Same link line the batch scripts use. correctionlib is not optional:
        # evaluate_btagSF{bc,light}.cc call correction::CorrectionSet, so every
        # configuration with a b-tag scale factor fails at link time without it.
        ld_flags  = ("$(shell root-config --libs) -lTMVA -lXMLIO "
                     "-I/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/include "
                     "-L/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/lib "
                     "-lcorrectionlib")
        targets   = [os.path.splitext(f)[0] for f in file_paths]
        with open(makefile_name, "w") as f:
            f.write("# Generated MakefileSnapshot\nCXX = g++\n")
            f.write(f"CXXFLAGS = -O2 -Wall {cpp_flags}\n")
            f.write(f"LDFLAGS = {ld_flags}\n\n")
            all_targets = " ".join(targets)
            f.write(f"all: {all_targets}\n\n")
            f.write("%: %.cpp\n\t$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)\n\n")
            f.write("clean:\n")
            f.write(f"\tfind {self._scripts_run_folder} -type f ! -name '*.cpp' ! -name '*.h' -delete\n")

    # _____________________________________________________________________________
    #  makeScripts — nominal run
    # _____________________________________________________________________________

    def makeScripts(self, patch_mode=False, force=False):
        print("=====================\n==== makeScripts ====\n=====================")
        if patch_mode:
            print("  -> PATCH MODE: adding missing variables to existing snapshots")

        os.makedirs(self._folder_where_to_save_trees, exist_ok=True)
        os.makedirs(self._scripts_run_folder, exist_ok=True)

        self._pre_count_events()

        list_of_files_to_compile = []

        for sampleName, sample in self._samples.items():
            os.makedirs(self._scripts_run_folder + "/" + sampleName, exist_ok=True)

            global_weight   = sample.get('weight', "1.")
            weights_subname = sample.get('weights', {})
            if 'isData' not in sample:
                global_weight = f"({global_weight}) * {self._lumi}"

            _, _, aliases_ordered = self._build_aliases_ordered(sampleName)

            # build full variables list (auto-adds cut vars, alias vars, weight nuisances)
            full_vtd, weight_nuisance_branches = self._build_full_variables_to_dump(
                aliases_ordered, self._nuisances)
            self._current_variables_to_dump = full_vtd

            # weight nuisance defines for this sample
            extra_weight_defines = self._build_weight_nuisance_defines(
                sampleName, global_weight, self._nuisances)

            for subname, list_root_files in sample['name'].items():
                print(f"  {sampleName}/{subname}: {len(list_root_files)} files")
                os.makedirs(self._scripts_run_folder + "/" + sampleName + "/" + subname, exist_ok=True)

                weight = (f"({global_weight}) * ({weights_subname[subname]})"
                          if subname in weights_subname else global_weight)

                chunks = self._get_chunks(sample, list_root_files)
                print(f"  -> {len(chunks)} jobs")

                for i, chunk in enumerate(chunks):
                    job_key  = f"{sampleName}_{subname}_{i}"
                    out_file = self.output_file_name(sampleName, subname, i)
                    out_path = os.path.join(self._folder_where_to_save_trees, out_file)

                    if not force and not patch_mode and self._metadata.is_completed(job_key, out_path):
                        print(f"    -> {job_key} already done, skipping")
                        continue

                    name_code = (self._scripts_run_folder + "/" + sampleName + "/" + subname
                                 + "/" + self.job_prefix() + "_" + sampleName + "_" + subname + "_" + str(i))

                    if patch_mode:
                        # determine missing variables
                        missing = self._metadata.missing_variables(full_vtd)
                        if not missing:
                            print(f"    -> {job_key} nothing to patch")
                            continue
                        # can we patch from snapshot or need nanoAOD?
                        if self._metadata.can_patch_from_snapshot(missing):
                            print(f"    -> {job_key} patch from snapshot ({len(missing)} new vars)")
                            input_files = [out_path]  # read existing snapshot
                            patch_suffix = f"_patch{len(self._metadata.get('patches', []))}"
                            patch_out    = self.output_file_name(sampleName, subname, i).replace('.root', f'{patch_suffix}.root')
                            self._current_variables_to_dump = missing
                        else:
                            print(f"    -> {job_key} patch requires nanoAOD (missing deps)")
                            input_files = chunk['files']
                            patch_out   = out_file
                    else:
                        input_files = chunk['files']
                        patch_out   = out_file

                    self.create_cpp_source_list_of_files(
                        name_code, "Events", input_files,
                        patch_out if patch_mode else out_file,
                        sampleName, weight,
                        {}, {}, aliases_ordered=aliases_ordered,
                        range_start=chunk.get('range_start'), range_end=chunk.get('range_end'),
                        extra_weight_defines=extra_weight_defines)

                    list_of_files_to_compile.append(name_code + ".cpp")

        # update metadata
        self._metadata.set('cut_to_dump',     self._cut_to_dump)
        self._metadata.set('cut_expression',  self._cuts.get(self._cut_to_dump, ''))
        self._metadata.set('cut_variables',   list(_collect_cut_variables(self._cuts)))
        self._metadata.set('variables_saved', self._current_variables_to_dump)
        self._metadata.set('alias_dependencies',
                           _collect_alias_dependencies(aliases_ordered if 'aliases_ordered' in dir() else {}))
        # track weight nuisances
        wn = [n for n, nu in self._nuisances.items()
              if nu.get('type') == 'shape' and nu.get('kind') == 'weight']
        self._metadata.set('weight_nuisances', wn)
        # track shape nuisances
        sn = [n for n, nu in self._nuisances.items()
              if nu.get('type') == 'shape' and nu.get('kind') in ('suffix', 'suffixFly')]
        self._metadata.set('shape_nuisances', sn)
        self._metadata.save()
        print(f"  -> Metadata saved: {self._metadata._path}")

        self.generate_makefile(list_of_files_to_compile)

        # The varied trees belong to the same run, so generate them here rather
        # than behind a flag -- otherwise submitBatch queues jobs whose .cpp
        # does not exist.
        if self.doShapeNuisances() and not patch_mode:
            self.makeShapeNuisances()
        subprocess.run(["make", "-f", "MakefileSnapshot", "clean"])
        print("Scripts generated. To compile: make -f MakefileSnapshot -j8")
        print("Or use --localCompile.")

    # _____________________________________________________________________________
    #  makeShapeNuisances — generate alternative TTree files for shape nuisances
    # _____________________________________________________________________________

    def makeShapeNuisances(self):
        """
        For each shape nuisance of kind 'suffix', generate a separate snapshot
        reading the alternative TTree (friend tree in nanoAOD).
        Output: nanoLatino_snapshot____<sample>_<subname>_<i>_<nuisanceName>Up.root
        """
        print("===============================")
        print("==== makeShapeNuisances    ====")
        print("===============================")

        os.makedirs(self._folder_where_to_save_trees, exist_ok=True)
        os.makedirs(self._scripts_run_folder, exist_ok=True)

        self._pre_count_events()

        list_of_files_to_compile = []

        for nuisanceName, nuisance in self._nuisances.items():
            if nuisance.get('type') != 'shape':
                continue
            if nuisance.get('kind') not in ('suffix', 'suffixFly'):
                continue

            print(f"  Nuisance: {nuisanceName}")

            for sampleName, sample in self._samples.items():
                if sampleName not in nuisance.get('samples', {}):
                    continue

                _, _, aliases_ordered = self._build_aliases_ordered(sampleName)
                full_vtd, _ = self._build_full_variables_to_dump(aliases_ordered, {})
                self._current_variables_to_dump = full_vtd

                global_weight   = sample.get('weight', "1.")
                weights_subname = sample.get('weights', {})
                if 'isData' not in sample:
                    global_weight = f"({global_weight}) * {self._lumi}"

                for subname, list_root_files in sample['name'].items():
                    chunks = self._get_chunks(sample, list_root_files)

                    for i, chunk in enumerate(chunks):
                        for direction in ['Up', 'Down']:
                            map_name = nuisance[f'map{direction}']
                            folder   = nuisance[f'folder{direction}']

                            # the varied files are a friend of the nominal ones,
                            # not a replacement for them: see
                            # _build_cpp_input_files_extra
                            varied_files = [os.path.join(folder, os.path.basename(f))
                                            for f in chunk['files']]

                            job_name, out_file = self._nuisance_names(
                                nuisanceName, direction, sampleName, subname, i)
                            name_code = os.path.join(self._scripts_run_folder, sampleName,
                                                     subname, job_name)

                            os.makedirs(os.path.dirname(name_code), exist_ok=True)

                            weight = (f"({global_weight}) * ({weights_subname[subname]})"
                                      if subname in weights_subname else global_weight)

                            self._suffix_variation = {
                                'nuisance' : nuisanceName,
                                'direction': direction,
                                'map'      : map_name,
                                'files'    : varied_files,
                            }
                            try:
                                self.create_cpp_source_list_of_files(
                                    name_code, "Events", chunk['files'], out_file,
                                    sampleName, weight,
                                    {}, {}, aliases_ordered=aliases_ordered,
                                    range_start=chunk.get('range_start'),
                                    range_end=chunk.get('range_end'))
                            finally:
                                self._suffix_variation = None

                            list_of_files_to_compile.append(name_code + ".cpp")

        # update metadata
        sn = [n for n, nu in self._nuisances.items()
              if nu.get('type') == 'shape' and nu.get('kind') in ('suffix', 'suffixFly')]
        self._metadata.set('shape_nuisances', sn)
        self._metadata.save()

        self.generate_makefile(list_of_files_to_compile, "MakefileSnapshotNuisances")
        print("Scripts generated. To compile: make -f MakefileSnapshotNuisances -j8")

    # _____________________________________________________________________________
    #  haddBatch — hadd all snapshot files into one per tag
    # _____________________________________________________________________________

    def haddBatch(self, tag, haddnano_path=None):
        if haddnano_path is None:
            haddnano_path = os.path.join(_SCRIPTS_DIR, "haddnano.py")
        super().haddBatch(
            tag          = tag,
            output_dir   = self._folder_where_to_save_trees,
            name_prefix  = "nanoLatino_snapshot____",
            output_name  = f"snapshot_{tag}.root",
            hadd_tool_path = haddnano_path,
            hadd_tool_args = "",   # haddnano doesn't need extra args
        )

    # _____________________________________________________________________________
    #  submitBatch hooks
    # _____________________________________________________________________________

    def _extra_transfer_files_for_submit(self, precompiled, credentials, extra_files_str):
        return f"{credentials}{extra_files_str}"

    def _build_bash_precompiled(self, job_name, out_file, sampleName, subname,
                                proxy_lines, submission_dir, root_cflags, root_libs):
        cp_cmd = self._build_cp_cmd(out_file, self._folder_where_to_save_trees)
        return f"""#!/bin/bash
set -e
echo "Job started at $(date)"; echo "Running on node $(hostname)"
source /cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/setup.sh
{proxy_lines}./{job_name}
{cp_cmd}
rm -f {out_file} {job_name}
echo "Done at $(date)"; ls -lh
"""

    def _build_bash_compile_on_node(self, job_name, out_file, sampleName, subname,
                                    proxy_lines, submission_dir, root_cflags, root_libs):
        cp_cmd = self._build_cp_cmd(out_file, self._folder_where_to_save_trees)
        return f"""#!/bin/bash
set -e
echo "Job started at $(date)"; echo "Running on node $(hostname)"
source /cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/setup.sh
{proxy_lines}echo "Compiling..."
g++ -O2 -Wall {root_cflags} \\
    {job_name}.cpp -o {job_name} \\
    -lTMVA -lXMLIO \\
    -I/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/include \\
    -L/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/lib \\
    -lcorrectionlib {root_libs}
echo "Compilation done at $(date)"
./{job_name}
{cp_cmd}
rm -f {out_file} {job_name}
echo "Done at $(date)"; ls -lh
"""

    # _____________________________________________________________________________
    #  checkBatch — output file is in folder_where_to_save_trees, not outputDir
    # _____________________________________________________________________________

    def checkBatch(self, resubmit=False):
        """Override to also mark completed jobs in metadata."""
        super().checkBatch(resubmit=resubmit)
        # update metadata with completed files
        for sampleName in self._samples:
            sample_dir = os.path.join(self._script_batch_location, sampleName)
            if not os.path.isdir(sample_dir):
                continue
            for subname in os.listdir(sample_dir):
                subname_dir = os.path.join(sample_dir, subname)
                if not os.path.isdir(subname_dir):
                    continue
                sh_files = sorted(f for f in os.listdir(subname_dir) if f.endswith('.sh'))
                for sh_file in sh_files:
                    job_name = sh_file[:-3]
                    suffix   = job_name.replace(self.job_prefix() + "_", "", 1)
                    out_path = self._output_file_from_suffix(sampleName, subname, suffix)
                    job_key  = suffix
                    if os.path.exists(out_path):
                        self._metadata.mark_completed(job_key, out_path)
        self._metadata.save()


# -----------------------------------------------------------------------
#  Main
# -----------------------------------------------------------------------

if __name__ == '__main__':
    sys.argv = argv

    header = """
         --------------------------------------------------------------------------------------------------
         '      ___|                              |             |         \\  |         |                  '
         '    \\___ \\   __ \\    _` |  __ \\    __|  __ \\    _ \\   __|      |\\/ |   _` |  |  /   _ \\   __|   '
         '          |  |   |  (   |  |   | \\__ \\  | | |  (   |  |        |   |  (   |    <    __/  |      '
         '    _____/  _|  _| \\__,_|  .__/  ____/ _| |_| \\___/  \\__|     _|  _| \\__,_| _|\\_\\ \\___| _|      '
         --------------------------------------------------------------------------------------------------
    """
    print(header)

    parser = defaultParser()
    parser.add_argument("--localCompile",    action='store_true', dest="localCompile")
    parser.add_argument("--submitBatch",     action='store_true', dest="submitBatch")
    parser.add_argument("--checkBatch",      action='store_true', dest="checkBatch")
    parser.add_argument("--resubmit",        action='store_true', dest="resubmit")
    parser.add_argument("--silentMode",      action='store_true', dest="silentMode")
    parser.add_argument("--haddBatch",       action='store_true', dest="haddBatch",
                        help="Submit a single condor job to hadd all snapshot files")
    parser.add_argument("--patch",           action='store_true', dest="patch",
                        help="Add missing variables to existing snapshots")
    parser.add_argument("--shapeNuisances",  action='store_true', dest="shapeNuisances",
                        help="No longer needed: the varied trees are produced "
                             "whenever nuisances.py declares a 'suffix' nuisance. "
                             "Accepted and ignored.")
    parser.add_argument("--noShapeNuisances", action='store_true', dest="noShapeNuisances",
                        help="Nominal snapshot only, ignoring the 'suffix' "
                             "nuisances declared in nuisances.py")
    parser.add_argument("--force",           action='store_true', dest="force",
                        help="Force regeneration even if already done")

    opt = parser.parse_args()
    loadDefaultOptions(parser, opt.pycfg)
    opt = parser.parse_args()

    cfg = load_config_files(opt)
    samples   = cfg['samples']
    variables = cfg['variables']
    cuts      = cfg['cuts']
    supercut  = cfg['supercut']
    aliases   = cfg['aliases']
    nuisances = cfg['nuisances']

    # snapshot.py
    snap_ns = {**_base_ns(), 'samples': samples, 'cuts': cuts, 'aliases': aliases,
               'variables_to_dump': [], 'cut_to_dump': '', 'folder_where_to_save_trees': ''}
    if hasattr(opt, 'snapshotFile') and os.path.exists(opt.snapshotFile):
        with open(opt.snapshotFile, 'r') as fh:
            exec(fh.read(), snap_ns)
    variables_to_dump = snap_ns.get('variables_to_dump', [])
    cut_to_dump       = snap_ns.get('cut_to_dump', '')

    folder_where_to_save_trees = snap_ns.get('folder_where_to_save_trees', '') or \
                                 getattr(opt, 'folder_where_to_save_trees', '') or \
                                 os.path.join(os.path.expandvars(opt.outputDir), 'snapshot')
    folder_where_to_save_trees = os.path.expandvars(folder_where_to_save_trees)

    factory = SnapshotFactory()
    factory.setValues(
        folder_where_to_save_trees,
        variables, cuts, supercut, samples, nuisances, aliases, getattr(opt, 'lumi', 1),
        scripts_run_folder    = getattr(opt, 'snapshotRunFolder',   None),
        script_batch_location = getattr(opt, 'batchScriptLocation', None),
        extra_transfer_files  = getattr(opt, 'extra_transfer_files', None),
    )
    factory.setSnapshotValues(variables_to_dump, cut_to_dump, folder_where_to_save_trees)
    factory.setConditions(opt.silentMode)

    factory.setSkipShapeNuisances(opt.noShapeNuisances)
    if opt.shapeNuisances:
        print("  NOTE: --shapeNuisances is no longer needed; the varied trees are\n"
              "        produced whenever nuisances.py declares a 'suffix' nuisance.\n"
              "        Use --noShapeNuisances to opt out.")

    n_suffix = len(factory._suffixNuisances())
    if n_suffix:
        if factory.doShapeNuisances():
            print(f"  -> {n_suffix} shape nuisance(s): {2 * n_suffix} varied tree(s) per job.")
        else:
            print(f"  -> {n_suffix} shape nuisance(s) declared, skipped (--noShapeNuisances).")

    # Generating the scripts is the prerequisite of compiling and of submitting,
    # and it is cheap and idempotent (_write_if_changed), so it runs for those
    # too. It used to be a separate bare invocation, so --submitBatch on a clean
    # job directory queued jobs whose .cpp did not exist.
    # --checkBatch and --haddBatch only read what is already there.
    if opt.patch:
        factory.makeScripts(patch_mode=True, force=opt.force)
    elif not any([opt.checkBatch, opt.haddBatch]):
        factory.makeScripts(force=opt.force)

    if opt.localCompile:
        subprocess.run(["make", "-f", "MakefileSnapshot", "-j8"])
        # makeScripts wrote this one too, so build it in the same breath
        if factory.doShapeNuisances():
            subprocess.run(["make", "-f", "MakefileSnapshotNuisances", "-j8"])

    if opt.submitBatch:
        factory.submitBatch()

    if opt.haddBatch:
        tag = getattr(opt, 'tag', 'snapshot')
        haddnano_path = getattr(opt, 'haddnanoPath', None)
        factory.haddBatch(tag, haddnano_path=haddnano_path)

    if opt.checkBatch:
        factory.checkBatch(resubmit=opt.resubmit)

    print("\n\n I'm done ...\n\n")