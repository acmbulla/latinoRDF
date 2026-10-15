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

    # _____________________________________________________________________________
    def setSnapshotValues(self, variables_to_dump, cut_to_dump, folder_where_to_save_trees):
        self._variables_to_dump          = variables_to_dump
        self._cut_to_dump                = cut_to_dump
        self._folder_where_to_save_trees = os.path.abspath(
            os.path.expandvars(folder_where_to_save_trees))

    # _____________________________________________________________________________
    def job_prefix(self):
        return "my_snapshot"

    def output_file_name(self, sampleName, subname, i):
        return f"nanoLatino_snapshot____{sampleName}_{subname}_{i}.root"

    def _output_file_from_suffix(self, sampleName, subname, suffix):
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

        mt_setting = "ROOT::DisableImplicitMT();" if (range_start is not None and range_end is not None) else "ROOT::EnableImplicitMT();"

        snapshot_call = f'node_{self._cut_to_dump}.Snapshot("Events", "{output_root_file_name}", columns_to_save);\n'

        cpp_code = f"""
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RVec.hxx"
#include "TFile.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TChain.h"
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
        with open(f"{output_name}.cpp", "w") as f:
            f.write(cpp_code)

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

    def makeScripts(self):
        print("=====================\n==== makeScripts ====\n=====================")
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
                    out_file = self.output_file_name(sampleName, subname, i)

                    name_code = (self._scripts_run_folder + "/" + sampleName + "/" + subname
                                 + "/" + self.job_prefix() + "_" + sampleName + "_" + subname + "_" + str(i))

                    self.create_cpp_source_list_of_files(
                        name_code, "Events", chunk['files'],
                        out_file,
                        sampleName, weight,
                        {}, {}, aliases_ordered=aliases_ordered,
                        range_start=chunk.get('range_start'), range_end=chunk.get('range_end'),
                        extra_weight_defines=extra_weight_defines)

                    list_of_files_to_compile.append(name_code + ".cpp")

        self.generate_makefile(list_of_files_to_compile)
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

                            # build alternative file list
                            alt_files = []
                            for f_path in chunk['files']:
                                fname = os.path.basename(f_path)
                                alt_files.append(os.path.join(folder, fname))

                            job_name  = (f"my_snapshot_{nuisanceName}{direction}_"
                                         f"{sampleName}_{subname}_{i}")
                            out_file  = (f"nanoLatino_snapshot____{sampleName}_{subname}_{i}"
                                         f"_{nuisanceName}{direction}.root")
                            name_code = os.path.join(self._scripts_run_folder, sampleName,
                                                     subname, job_name)

                            os.makedirs(os.path.dirname(name_code), exist_ok=True)

                            weight = (f"({global_weight}) * ({weights_subname[subname]})"
                                      if subname in weights_subname else global_weight)

                            self.create_cpp_source_list_of_files(
                                name_code, "Events", alt_files, out_file,
                                sampleName, weight,
                                {}, {}, aliases_ordered=aliases_ordered,
                                range_start=chunk.get('range_start'),
                                range_end=chunk.get('range_end'))

                            list_of_files_to_compile.append(name_code + ".cpp")

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

    def _extra_transfer_files_for_submit(self, precompiled, x509_proxy, extra_files_str):
        return f"{x509_proxy}{extra_files_str}"

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

    def _output_file_from_suffix(self, sampleName, subname, suffix):
        return os.path.join(self._folder_where_to_save_trees,
                            f"nanoLatino_snapshot____{suffix}.root")


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
    parser.add_argument("--shapeNuisances",  action='store_true', dest="shapeNuisances",
                        help="Generate alternative TTree files for shape nuisances")

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

    if not any([opt.submitBatch, opt.checkBatch, opt.localCompile,
                opt.shapeNuisances, opt.haddBatch]):
        factory.makeScripts()

    if opt.shapeNuisances:
        factory.makeShapeNuisances()

    if opt.localCompile:
        factory.makeScripts()
        subprocess.run(["make", "-f", "MakefileSnapshot", "-j8"])

    if opt.submitBatch:
        factory.submitBatch()

    if opt.haddBatch:
        tag = getattr(opt, 'tag', 'snapshot')
        haddnano_path = getattr(opt, 'haddnanoPath', None)
        factory.haddBatch(tag, haddnano_path=haddnano_path)

    if opt.checkBatch:
        factory.checkBatch(resubmit=opt.resubmit)

    print("\n\n I'm done ...\n\n")