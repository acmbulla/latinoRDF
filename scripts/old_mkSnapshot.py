#!/usr/bin/env python3

import sys
argv = sys.argv
sys.argv = argv[:1]

import argparse

import os

import ROOT


#
# default parser
#

def defaultParser():
    sys.argv = argv

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        "--pycfg",
        dest="pycfg",
        help="input configuration file",
        default="configuration.py",
    )

    return parser


#
# load configurations into the parser, according to configuration.py
#

def loadDefaultOptions(parser, pycfg=None, quiet=False):

    if pycfg != None and os.path.exists(pycfg):
        # print ("loadDefaultOptions: pycfg = ", pycfg)
        handle = open(pycfg,'r')
        local_variables = {}
        exec(handle.read(),local_variables)
        handle.close()
        # clean the dictionary to remove globals due to "exec" funcionality
        local_variables = {k: v for k, v in local_variables.items() if not (k.startswith('__') and k.endswith('__'))}
        for opt_name, opt_value in local_variables.items():
          parser.add_argument('--' + opt_name, default=opt_value)
        return
    else:
        return




# ----------------------------------------------------- SnapshotFactory --------------------------------------

import subprocess

class SnapshotFactory:

    # _____________________________________________________________________________
    def __init__(self):

        variables = {}   # FIXME needed?
        self._variables = variables

        cuts = {}
        self._cuts = cuts

        self._supercut = ''

        samples = {}
        self._samples = samples

        aliases = {}
        self._aliases = aliases

        self._lumi = 1

        self._treeName = 'latino'

        self._outputRootFile = 'myhistos.root'

        self._script_batch_location = './scripts_snapshot_batch/'

        self._scripts_run_folder = "scripts_snapshot_run"


        self._variables_to_dump           = []
        self._cut_to_dump                 = ''
        self._folder_where_to_save_trees  = ''


        # conditions
        self._silentMode = False



    # _____________________________________________________________________________
    def __del__(self):
        pass


    # _____________________________________________________________________________
    def unroll_cuts(self):

      temp_cuts = {}

      for cutName, cut in self._cuts.items():
        #
        # cuts could be nested, i.e. cut -> categories
        # here is where you exploit the full power of RDataFrame
        #
        # Procedure:
        # - check if categories
        # - if yes, then define the cuts in a nested way
        # - otherwise if "expr" is defined use it
        #   if note use directly "cut"
        #
        # Different possibilities:
        #
        # cut["DY"] = "mll>50 && mll<120"
        # cut["DY"] = {   'expr': 'mll>70 && mll<110' }
        # cut["DY"] = {
        #    'expr': 'mll>70 && mll<110',
        #    'categories' : {
        #       'eleele' : 'ee',   # "ee" is defined in aliases.py
        #       'mumu'   : 'mm',   # "mm" is defined in aliases.py
        #      }
        #    }
        #


        list_cuts = {}
        if isinstance(cut, dict):
          # "cut" is a dictionary
          expression = cut['expr']
          categories = cut.get('categories', {})
          if categories :
            temp_cuts[ cutName ] = "(" + expression + ")"   # add the un-categorized phase space too! it comes "for free"
            for category_name, category in categories.items():
              temp_cuts[ cutName + "_" + category_name] = "(" + expression + ") && (" + category + ")"
          else :
            temp_cuts[ cutName ] = expression
        else :
          temp_cuts[ cutName ] = cut

      self._cuts = temp_cuts



    # _____________________________________________________________________________
    #                                                              "file_path" is a list of root files!
    def create_cpp_source_list_of_files(self, output_name, tree_name, file_path, output_root_file_name, sampleName, weight,
                                        aliases_to_be_defined, functions_to_be_defined,
                                        range_start=None, range_end=None):

        define_input_files_logic = ""
        define_input_files_logic += f"//  Nominal input files\n"
        define_input_files_logic += f'    auto* nominal = new TChain("{tree_name}");\n'

        for file_name in file_path:
          define_input_files_logic += f'    nominal->Add("{file_name}");\n'

        define_input_files_logic += f'    \n'
        define_input_files_logic += f'    ROOT::RDataFrame base_df(*nominal);\n'
        define_input_files_logic += f'    auto varied_df = ROOT::RDF::RNode(base_df);\n'

        #
        # define the weights needed for this specific sample
        #
        define_weights = ""
        define_weights += f"//  weigths needed\n"
        define_weights += f'    varied_df = SafeDefine(varied_df, "my_sample_weight", "{weight}");\n'

        #
        # define the aliases needed for this specific sample
        #
        define_aliases = ""
        define_aliases += f"//  aliases/define needed\n"

        for aliasName, alias in aliases_to_be_defined.items():
          define_aliases += f'    varied_df = SafeDefine(varied_df, "{aliasName}", "{alias}");\n'

        #
        # define the aliases and functions needed for this specific sample
        #
        code_of_function_to_include = ""

        define_aliases_functions = ""
        define_aliases_functions += f"//  aliases of functions\n"

        for aliasName, alias in functions_to_be_defined.items():
          if 'external' in alias.keys():
            with open(alias["external"], "r") as file:
              code_of_function_to_include += file.read()
          code_of_function_to_include += "\n"

          define_aliases_functions += "    varied_df = varied_df"

          if 'xmlfile' not in alias.keys():
            for i, var in enumerate(alias["variables"]):
              define_aliases_functions += f'\n                     .Define("_var_{aliasName}_{i}", "{var}")'
          elif 'xmlfile' in alias.keys():
            for i, var in enumerate(alias["variables"]):
              define_aliases_functions += f'\n                     .Define("_var_{aliasName}_{i}_f", "(float){var}")'

          define_aliases_functions += ";\n"

          if 'xmlfile' not in alias.keys():
            define_aliases_functions += f'\n'
            new_names = [f"_var_{aliasName}_{i}" for i in range(len(alias["variables"]))]
            sintax_variables = '{ "' + '", "'.join(new_names) + '" }'
            define_aliases_functions += f'    varied_df = varied_df.Define("{aliasName}", {alias["function"]}, {sintax_variables} );\n'
            define_aliases_functions += f'\n'
          elif 'xmlfile' in alias.keys():
            define_aliases_functions += f'\n'
            submission_dir = os.getcwd()
            define_aliases_functions += f'''    TMVA::Experimental::RReader my_model_{aliasName}("{submission_dir}/{alias['xmlfile']}");\n'''

            new_names = [f"_var_{aliasName}_{i}_f" for i in range(len(alias["variables"]))]
            sintax_variables = '{ "' + '", "'.join(new_names) + '" }'

            define_aliases_functions += f'''    varied_df = varied_df.Define("{aliasName}",
                     TMVA::Experimental::Compute<{len(alias['variables'])}, float>(my_model_{aliasName}),
                     {sintax_variables});\n'''
            define_aliases_functions += f'\n'


        booking_logic = ""
        booking_logic += f"//  Initial ...\n"
        booking_logic += f'    auto current_node = ROOT::RDF::RNode(varied_df);\n'
        if range_start is not None and range_end is not None:
          booking_logic += f'    current_node = ROOT::RDF::RNode(current_node.Range({range_start}, {range_end}));\n'
          booking_logic += f'    std::cout << "   Range: [{range_start}, {range_end})" << std::endl;\n'


        define_variables_logic = ""
        define_variables_logic += f"//  I need to define the RNode, otherwise SafeDefine will not work\n"
        define_variables_logic += f"    // Define variables \n"
        for variableName, variable in self._variables.items():
          if 'is2d' in variable.keys() and variable['is2d'] == 1:
            pass
          else:
            # only for 1D variables
            name = variable['name']
            define_variables_logic += f'    current_node = SafeDefine(current_node, "{variableName}", "{name}");\n'

        #
        # if "supercut" is defined, use it to speed up
        #
        if self._supercut != '' :
          define_variables_logic += f'    current_node = current_node.Filter("{self._supercut}", "supercut");\n'


        # remove nested definition of cuts
        self.unroll_cuts()

        define_snapshot_logic = ''
        # check if the cut to dump is in the list
        if self._cut_to_dump in self._cuts.keys() :

          this_cutName = self._cut_to_dump
          this_cut = self._cuts[self._cut_to_dump]

          define_variables_logic += f'\n'
          define_variables_logic += f'    auto node_{this_cutName} = current_node.Filter("{this_cut}", "{this_cutName}");\n'
          define_variables_logic += f'\n'

          # save only the columns that exist, e.g. if data there is no "promptgenmatched"
          columns_str = ", ".join([f'"{v}"' for v in self._variables_to_dump])
          define_variables_logic += f'    std::vector<std::string> requested = {{{columns_str}}};\n'
          define_variables_logic += f'    std::vector<std::string> existing = node_{this_cutName}.GetColumnNames();\n'
          define_variables_logic += f'    std::vector<std::string> columns_to_save;\n'
          define_variables_logic += f'    for (const auto& col : requested) {{                                          \n'
          define_variables_logic += f'      if (std::find(existing.begin(), existing.end(), col) != existing.end()) {{  \n'
          define_variables_logic += f'        columns_to_save.push_back(col);                                           \n'
          define_variables_logic += f'      }}                                                                          \n'
          define_variables_logic += f'    }}                                                                            \n'
          define_variables_logic += f'\n'

          define_snapshot_logic += f'node_{this_cutName}.Snapshot("Events", "{output_root_file_name}", columns_to_save);   \n'



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
    if (std::find(colNames.begin(), colNames.end(), name) == colNames.end()) {{
        return df.Define(name, expr);
    }}
    return df;
}}



// --- Automatically generated: code to be added for additional functions ---

{code_of_function_to_include}

// ----------------------------------------




int main() {{
    ROOT::EnableImplicitMT();

    // --- Automatically generated input root files ---

    {define_input_files_logic}

    // ----------------------------------------

    // --- Automatically generated define aliases ---

    {define_aliases}

    // ----------------------------------------


    // --- Automatically generated define aliases for functions, no JIT ---

    {define_aliases_functions}

    // ----------------------------------------


    // --- Automatically generated define weights ---

    {define_weights}

    // ----------------------------------------


    // ----------------------------------------

    // --- Automatically generated ... just the main node ---

    {booking_logic}

    //
    // from now on: current_node
    // ----------------------------------------

    // --- Automatically generated define cuts what to snapshot ---

    {define_variables_logic}

    // ----------------------------------------


    // --- Automatically generated: perform snapshot ---

    {define_snapshot_logic}

    // ----------------------------------------


    return 0;
}}

        """


        with open(f"{output_name}.cpp", "w") as f:
            f.write(cpp_code)


    # _____________________________________________________________________________
    def generate_makefile(self, file_paths, makefile_name="MakefileSnapshot"):
        # Get ROOT configuration via shell calls
        cpp_flags = "$(shell root-config --cflags)"
        ld_flags = "$(shell root-config --libs) -lTMVA -lXMLIO"

        type_of_compilation = "-O2"
        # from "blabla.cpp" to "blabla"
        targets = [os.path.splitext(f)[0] for f in file_paths]

        with open(makefile_name, "w") as f:
            f.write("# Generated Makefile\n")
            f.write("CXX = g++\n")
            f.write(f"CXXFLAGS = {type_of_compilation} -Wall {cpp_flags} \n")
            f.write(f"LDFLAGS = {ld_flags}\n\n")

            # 'all' target: the list of all executables to be created
            all_targets = " ".join(targets)
            f.write(f"all: {all_targets}\n\n")

            # Pattern rule: works across subdirectories automatically
            f.write("%: %.cpp\n")
            f.write("\t$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)\n\n")

            # Clean target — use find to avoid "Argument list too long" with many targets
            f.write("clean:\n")
            f.write(f"\tfind {self._scripts_run_folder} -type f ! -name '*.cpp' ! -name '*.h' -delete\n")


    # _____________________________________________________________________________
    def setValues(self, variables, cuts, supercut, samples, aliases,
                  variables_to_dump, cut_to_dump, folder_where_to_save_trees,
                  batch_script_location='./scripts_snapshot_batch/',
                  scripts_run_folder=None):

        self._variables = variables
        self._samples   = samples
        self._cuts      = cuts
        self._supercut  = supercut
        self._aliases   = aliases

        self._variables_to_dump           = variables_to_dump
        self._cut_to_dump                 = cut_to_dump
        self._folder_where_to_save_trees  = os.path.expandvars(folder_where_to_save_trees)
        self._script_batch_location       = os.path.expandvars(batch_script_location)

        if scripts_run_folder is not None:
            self._scripts_run_folder = os.path.expandvars(scripts_run_folder)


    # _____________________________________________________________________________
    def setConditions(self, silentMode):

        self._silentMode = silentMode


    # _____________________________________________________________________________
    @staticmethod
    def _get_nevents(file_path, tree_name="Events"):
        try:
            f = ROOT.TFile.Open(file_path)
            if not f or f.IsZombie():
                print(f"  WARNING: cannot open {file_path}, assuming 0 events")
                return 0
            t = f.Get(tree_name)
            if not t:
                print(f"  WARNING: tree '{tree_name}' not found in {file_path}, assuming 0 events")
                f.Close()
                return 0
            n = t.GetEntries()
            f.Close()
            return n
        except Exception as e:
            print(f"  WARNING: exception reading {file_path}: {e}")
            return 0


    # _____________________________________________________________________________
    @staticmethod
    def _get_nevents_worker(args):
        """Worker for multiprocessing: returns (file_path, mtime, nevents)."""
        file_path, tree_name = args
        try:
            mtime = os.path.getmtime(file_path)
        except Exception:
            mtime = 0
        n = SnapshotFactory._get_nevents(file_path, tree_name)
        return file_path, mtime, n


    # _____________________________________________________________________________
    def _load_events_cache(self):
        import json
        base = os.path.expandvars("${MY_PATH_JOBS}")
        cache_path = os.path.join(base, "events_cache.json")
        if os.path.exists(cache_path):
            try:
                with open(cache_path, 'r') as f:
                    return json.load(f), cache_path
            except Exception:
                pass
        return {}, cache_path


    # _____________________________________________________________________________
    def _save_events_cache(self, cache, cache_path):
        import json
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        with open(cache_path, 'w') as f:
            json.dump(cache, f, indent=2)


    # _____________________________________________________________________________
    def _count_events(self, all_files, tree_name="Events", n_workers=8):
        """
        Count events for all_files using multiprocessing + JSON cache.
        Cache key: file_path -> {mtime, nevents}. If mtime matches, skip recount.
        """
        import multiprocessing

        cache, cache_path = self._load_events_cache()

        results  = {}
        to_count = []

        for fp in all_files:
            try:
                mtime = os.path.getmtime(fp)
            except Exception:
                mtime = 0
            entry = cache.get(fp)
            if entry and abs(entry.get('mtime', -1) - mtime) < 1.0:
                results[fp] = entry['nevents']
            else:
                to_count.append(fp)

        if to_count:
            print(f"  -> Counting events in {len(to_count)} files ({len(results)} from cache) using {n_workers} workers...")
            with multiprocessing.Pool(processes=n_workers) as pool:
                counted = pool.map(SnapshotFactory._get_nevents_worker,
                                   [(fp, tree_name) for fp in to_count])
            for fp, mtime, n in counted:
                results[fp] = n
                cache[fp] = {'mtime': mtime, 'nevents': n}
                print(f"    {fp}  ->  {n} events")
            self._save_events_cache(cache, cache_path)
            print(f"  -> Cache saved: {cache_path}")
        else:
            print(f"  -> All {len(results)} files found in cache, skipping event counting.")

        return results


    # _____________________________________________________________________________
    @staticmethod
    def _make_event_chunks(list_root_files, events_per_job, nevents_map):
        """Group files into chunks using a precomputed {file: nevents} map."""
        chunks = []
        current_files  = []
        current_events = 0

        for file_path in list_root_files:
            n = nevents_map.get(file_path, 0)
            if n == 0:
                continue
            if n > events_per_job:
                if current_files:
                    chunks.append({'files': current_files, 'range_start': None, 'range_end': None, 'nevents': current_events})
                    current_files  = []
                    current_events = 0
                start = 0
                while start < n:
                    end = min(start + events_per_job, n)
                    chunks.append({'files': [file_path], 'range_start': start, 'range_end': end, 'nevents': end - start})
                    start = end
            else:
                if current_events + n > events_per_job and current_files:
                    chunks.append({'files': current_files, 'range_start': None, 'range_end': None, 'nevents': current_events})
                    current_files  = []
                    current_events = 0
                current_files.append(file_path)
                current_events += n

        if current_files:
            chunks.append({'files': current_files, 'range_start': None, 'range_end': None, 'nevents': current_events})

        return chunks


    # _____________________________________________________________________________
    def makeScripts(self):

        print ("=====================")
        print ("==== makeScripts ====")
        print ("=====================")

        os.makedirs(self._folder_where_to_save_trees, exist_ok=True)

        os.makedirs(self._scripts_run_folder, exist_ok=True)

        list_of_files_to_compile = []

        # collect all files that need event counting (EventsPerJob samples)
        all_files_to_count = []
        for sampleName, sample in self._samples.items():
          if 'EventsPerJob' in sample:
            for subname, list_root_files in sample['name'].items():
              all_files_to_count.extend(list_root_files)

        nevents_map = {}
        if all_files_to_count:
          print(f"\n  -> Pre-counting events for {len(all_files_to_count)} files...")
          nevents_map = self._count_events(all_files_to_count, tree_name="Events")
          print(f"  -> Event counting done.\n")

        for sampleName, sample in self._samples.items():
          os.makedirs(self._scripts_run_folder + "/" + sampleName, exist_ok=True)

          global_weight = sample['weight'] if 'weight' in sample.keys() else "1."
          weights_subname = sample['weights'] if 'weights' in sample.keys() else {}

          if 'isData' in sample.keys() :
            if not (len(sample ['isData']) == 1 and sample ['isData'][0] == 'all') :
              pass
          else :
            global_weight = f"({global_weight}) * {self._lumi}"

          aliases_to_be_defined = {}
          functions_to_be_defined = {}
          for aliasName, alias in self._aliases.items():
            if 'samples' in alias.keys():
              if sampleName in alias['samples'] and 'expr' in alias.keys():
                aliases_to_be_defined[aliasName] = alias['expr']
              if sampleName in alias['samples'] and 'function' in alias.keys():
                functions_to_be_defined[aliasName] = alias

          for subname, list_root_files in sample['name'].items():
            print (f"  {sampleName}/{subname}: {len(list_root_files)} files")
            os.makedirs(self._scripts_run_folder + "/" + sampleName + "/" + subname, exist_ok=True)

            if subname in weights_subname.keys():
              weight = f"({global_weight}) * ({weights_subname[subname]})"
            else :
              weight = global_weight

            tree = "Events"

            if 'EventsPerJob' in sample:
              chunks = self._make_event_chunks(list_root_files, sample['EventsPerJob'], nevents_map)
              print(f"  -> {len(chunks)} jobs created")
            else:
              make_job_every_N = sample.get('FilesPerJob', 1)
              chunks = [{'files': list_root_files[i:i + make_job_every_N], 'range_start': None, 'range_end': None, 'nevents': None}
                        for i in range(0, len(list_root_files), make_job_every_N)]

            for i, chunk in enumerate(chunks):
              name_code = self._scripts_run_folder + "/" + sampleName + "/" + subname + "/" + "my_snapshot_" + sampleName + "_" + subname + "_" + str(i)
              output_root_file_name = "nanoLatino_snapshot____" + sampleName + "_" + subname + "_" + str(i) + ".root"

              self.create_cpp_source_list_of_files(name_code, tree, chunk['files'], output_root_file_name, sampleName, weight,
                                                   aliases_to_be_defined, functions_to_be_defined,
                                                   range_start=chunk.get('range_start'),
                                                   range_end=chunk.get('range_end'))

              list_of_files_to_compile.append(name_code + ".cpp")


        # [MODIFIED] generate Makefile but do NOT compile automatically — same pattern as mkShapes
        self.generate_makefile(list_of_files_to_compile, "MakefileSnapshot")
        subprocess.run(["make", "-f", "MakefileSnapshot", "clean"])
        print("Scripts generated. To compile run:")
        print("    make -f MakefileSnapshot -j8")
        print("Or use --localCompile to compile automatically.")


    # _____________________________________________________________________________
    def checkBatch(self):

        print ("=====================")
        print ("==== check Batch ====")
        print ("=====================")
        print ("\n")

        submission_dir = os.getcwd()

        list_jobs_with_error = []
        list_missing_root_files = []
        #
        # Loop over samples
        #
        for sampleName, sample in self._samples.items():
          for subname, list_root_files in sample['name'].items():
            make_job_every_N = 1
            if "FilesPerJob" in sample.keys():
              make_job_every_N = sample['FilesPerJob']
            chunks_list_root_files = [list_root_files[i:i + make_job_every_N] for i in range(0, len(list_root_files), make_job_every_N)]
            for i, root_files in enumerate(chunks_list_root_files):
              name_err_file = self._script_batch_location + "/" + sampleName + "/" + subname + "/log/" + "my_script_" + str(i) + ".sh.err"
              # print (" name_err_file = ", name_err_file)
              if os.path.exists(name_err_file):
                with open(name_err_file, 'r') as f:
                  line_count = sum(1 for line in f)
                  if line_count > 1: # should I remove the warnings to have 0?
                    # FIXME
                    list_jobs_with_error.append(name_err_file)
              output_root_file_name = "nanoLatino_snapshot____" + sampleName + "_" + subname + "_" + str(i) + ".root"
              root_file_name = f"{self._folder_where_to_save_trees}/{output_root_file_name}"
              if os.path.exists(root_file_name):
                pass
              else :
                list_missing_root_files.append(root_file_name)


        if len(list_jobs_with_error) == 0:
          print (" No jobs with error ")
        else :
          print(" Jobs with errors:")
          for index, job in enumerate(list_jobs_with_error):
            print(f"    Job {index}: {job}")


        if len(list_missing_root_files) == 0:
          print (" All files are ready ")
        else :
          print(" Missing files:")
          for index, job in enumerate(list_missing_root_files):
            print(f"     File: {index}: {job}")




    def submitBatch(self):

        print ("======================")
        print ("==== submit Batch ====")
        print ("======================")

        os.makedirs(self._script_batch_location, exist_ok=True)

        submission_dir = os.getcwd()

        try:
            root_cflags = subprocess.check_output(["root-config", "--cflags"]).decode().strip()
            root_libs   = subprocess.check_output(["root-config", "--libs"]).decode().strip()
        except Exception:
            root_cflags = ""
            root_libs   = ""

        # detect mode: check if first expected executable exists
        first_exe = None
        for sampleName, sample in self._samples.items():
          for subname, list_root_files in sample['name'].items():
            if list_root_files:
              first_exe = os.path.join(self._scripts_run_folder, sampleName, subname,
                                       "my_snapshot_" + sampleName + "_" + subname + "_0")
              break
          if first_exe:
            break

        precompiled = first_exe is not None and os.path.isfile(first_exe)

        if precompiled:
            print("  -> Pre-compiled executables found: transferring binaries to worker nodes.")
        else:
            print("  -> No pre-compiled executables found: worker nodes will compile before running.")

        # X509 proxy check
        _proxy_candidates = [
            os.environ.get('X509_USER_PROXY', ''),
            os.path.expanduser('~/.proxy'),
            f"/tmp/x509up_u{os.getuid()}",
        ]
        x509_proxy = next((p for p in _proxy_candidates if p and os.path.exists(p)), None)
        if x509_proxy:
            proxy_check = subprocess.run(["voms-proxy-info", "--exists", "--valid", "1:00"],
                                         capture_output=True)
            if proxy_check.returncode != 0:
                print(f"ERROR: proxy found at {x509_proxy} but it is expired or invalid.")
                print("       Please run: voms-proxy-init --voms cms --valid 168:00")
                print("       then rerun --submitBatch")
                return
            proxy_basename = os.path.basename(x509_proxy)
            proxy_lines = f"export X509_USER_PROXY=$(pwd)/{proxy_basename}\nexport EOS_MGM_URL=root://eosuser.cern.ch\nvoms-proxy-info\n"
            print(f"  -> X509 proxy found and valid: {x509_proxy}")
        else:
            print("ERROR: no X509 proxy found.")
            print("       Please run: voms-proxy-init --voms cms --valid 168:00")
            print("       then rerun --submitBatch")
            return

        queue_rows = []

        # collect all files that need event counting (EventsPerJob samples)
        all_files_to_count = []
        for sampleName, sample in self._samples.items():
          if 'EventsPerJob' in sample:
            for subname, list_root_files in sample['name'].items():
              all_files_to_count.extend(list_root_files)

        nevents_map = {}
        if all_files_to_count:
          print(f"\n  -> Pre-counting events for {len(all_files_to_count)} files...")
          nevents_map = self._count_events(all_files_to_count, tree_name="Events")
          print(f"  -> Event counting done.\n")

        for sampleName, sample in self._samples.items():
          os.makedirs(os.path.join(self._script_batch_location, sampleName), exist_ok=True)

          for subname, list_root_files in sample['name'].items():
            os.makedirs(os.path.join(self._script_batch_location, sampleName, subname), exist_ok=True)

            if not self._silentMode:
              os.makedirs(os.path.join(self._script_batch_location, sampleName, subname, "log"), exist_ok=True)

            tree = "Events"
            if 'EventsPerJob' in sample:
              chunks = self._make_event_chunks(list_root_files, sample['EventsPerJob'], nevents_map)
            else:
              make_job_every_N = sample.get('FilesPerJob', 1)
              chunks = [{'files': list_root_files[i:i + make_job_every_N], 'range_start': None, 'range_end': None, 'nevents': None}
                        for i in range(0, len(list_root_files), make_job_every_N)]

            for i, chunk in enumerate(chunks):
              job_name              = "my_snapshot_" + sampleName + "_" + subname + "_" + str(i)
              name_bash             = os.path.join(self._script_batch_location, sampleName, subname, job_name + ".sh")
              output_root_file_name = "nanoLatino_snapshot____" + sampleName + "_" + subname + "_" + str(i) + ".root"
              name_folder_code      = os.path.join(self._scripts_run_folder, sampleName, subname)
              abs_executable        = os.path.join(name_folder_code, job_name)
              abs_cpp               = os.path.join(name_folder_code, job_name + ".cpp")

              if precompiled:
                bash_code = f"""#!/bin/bash
set -e
echo "Job started at $(date)"
echo "Running on node $(hostname)"
source /cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/setup.sh
{proxy_lines}mkdir -p {self._folder_where_to_save_trees}
./{job_name}
mv {output_root_file_name} {self._folder_where_to_save_trees}/
echo "Done at $(date)"
ls -lh
"""
                job_exe = abs_executable

              else:
                bash_code = f"""#!/bin/bash
set -e
echo "Job started at $(date)"
echo "Running on node $(hostname)"
source /cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/setup.sh
{proxy_lines}echo "Compiling snapshot executable..."
g++ -O2 -Wall {root_cflags} \\
    {job_name}.cpp -o {job_name} \\
    -lTMVA -lXMLIO \\
    -I/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/include \\
    -L/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/lib \\
    -lcorrectionlib {root_libs}
echo "Compilation done at $(date)"
mkdir -p {self._folder_where_to_save_trees}
./{job_name}
mv {output_root_file_name} {self._folder_where_to_save_trees}/
echo "Done at $(date)"
ls -lh
"""
                job_exe = abs_cpp

              with open(name_bash, "w") as f:
                f.write(bash_code)
              os.system("chmod +x " + name_bash)

              name_folder = os.path.join(self._script_batch_location, sampleName, subname)

              if self._silentMode:
                output_file = "/dev/null"
                error_file  = "/dev/null"
                log_file    = "/dev/null"
              else:
                log_dir     = os.path.join(name_folder, "log")
                output_file = os.path.join(log_dir, job_name + ".sh.out")
                error_file  = os.path.join(log_dir, job_name + ".sh.err")
                log_file    = os.path.join(log_dir, job_name + ".sh.log")

              queue_rows.append((name_bash, job_exe, output_file, error_file, log_file))

        # single submit file
        global_submit_path = os.path.join(self._script_batch_location, "submit_all.sub")

        extra_files = getattr(self, '_extra_transfer_files', [])
        extra_files_str = (", " + ", ".join(extra_files)) if extra_files else ""

        fixed_extra = f"{x509_proxy}{extra_files_str}"

        queue_lines = "\n".join(
            f"  {exe}, {job_exe}, {out}, {err}, {log}"
            for exe, job_exe, out, err, log in queue_rows
        )

        submit_code = (
            "# Auto-generated by mkSnapshot.py\n"
            f"# Mode: {'pre-compiled' if precompiled else 'compile-on-node'}\n"
            "\n"
            "should_transfer_files   = YES\n"
            "when_to_transfer_output = ON_EXIT\n"
            f"transfer_input_files    = $(job_exe), {fixed_extra}\n"
            "getenv                  = True\n"
            "+JobFlavour             = \"longlunch\"\n"
            "\n"
            "queue executable, job_exe, output, error, log from (\n"
            + queue_lines + "\n"
            ")\n"
        )

        with open(global_submit_path, "w") as f:
            f.write(submit_code)

        print(f"Created: {global_submit_path}  ({len(queue_rows)} jobs)")
        print("Submitting all jobs with a single condor_submit ...")
        subprocess.run(["condor_submit", global_submit_path],
                       cwd=self._script_batch_location)




#
#
#  Main
#
#

if __name__ == '__main__':
    sys.argv = argv

    header = """
         --------------------------------------------------------------------------------------------------
         '                                                                                                '
         '      ___|                              |             |         \\  |         |                  '
         '    \\___ \\   __ \\    _` |  __ \\    __|  __ \\    _ \\   __|      |\\/ |   _` |  |  /   _ \\   __|   '
         '          |  |   |  (   |  |   | \\__ \\  | | |  (   |  |        |   |  (   |    <    __/  |      '
         '    _____/  _|  _| \\__,_|  .__/  ____/ _| |_| \\___/  \\__|     _|  _| \\__,_| _|\\_\\ \\___| _|      '
         '                          _|                                                                    '
         '          --------------------------------------------------------------------------------------------------
         """

    print(header)

    parser = defaultParser()

    # [MODIFIED] added --localCompile, aligned with mkShapes
    parser.add_argument("--localCompile",  action='store_true', dest="localCompile",   help="Compile locally with: make -f MakefileSnapshot -j8")
    parser.add_argument("--submitBatch",   action='store_true', dest="submitBatch",    help="Trigger the submission to lxbatch")
    parser.add_argument("--checkBatch",    action='store_true', dest="checkBatch",     help="Check if jobs are done succesfully")
    parser.add_argument("--silentMode",    action='store_true', dest="silentMode",     help="Remove as much as possible the print and the log/err files production")


    opt = parser.parse_args()
    print ("opt.pycfg            = ", opt.pycfg)
    loadDefaultOptions(parser, opt.pycfg)
    opt = parser.parse_args()


    print ("opt.localCompile     = ", opt.localCompile)
    print ("opt.submitBatch      = ", opt.submitBatch)
    print ("opt.checkBatch       = ", opt.checkBatch)
    print ("opt.silentMode       = ", opt.silentMode)


    variables_to_dump = []


    #
    # logic of dependencies:
    #
    #    samples.py
    #      --> aliases.py     since the aliases/defines will be "Defined" only for selected samples
    #    cuts.py
    #    variables.py
    #    snapshot.py
    #


    #
    # read list of samples
    #
    samples = {}
    if os.path.exists(opt.samplesFile) :
      handle = open(opt.samplesFile,'r')
      exec(handle.read())
      handle.close()
      # clean the dictionary to remove globals due to "exec" funcionality
      samples = {k: v for k, v in samples.items() if not (k.startswith('__') and k.endswith('__'))}

    # print ("samples = ", samples)

    #
    # read list of variables
    #
    variables = {}
    if os.path.exists(opt.variablesFile) :
      handle = open(opt.variablesFile,'r')
      exec(handle.read())
      handle.close()
      # clean the dictionary to remove globals due to "exec" funcionality
      variables = {k: v for k, v in variables.items() if not (k.startswith('__') and k.endswith('__'))}

    print ("variables = ", variables)


    #
    # read list of cuts
    #
    cuts = {}
    supercut = ''
    if os.path.exists(opt.cutsFile) :
      handle = open(opt.cutsFile,'r')
      exec(handle.read())
      handle.close()
      # clean the dictionary to remove globals due to "exec" funcionality
      cuts = {k: v for k, v in cuts.items() if not (k.startswith('__') and k.endswith('__'))}

    print ("cuts = ", cuts)
    print ("supercut = ", supercut)



    #
    # read list of aliases
    #
    aliases = {}
    if os.path.exists(opt.aliasesFile) :
      handle = open(opt.aliasesFile,'r')
      exec(handle.read())
      handle.close()
      # clean the dictionary to remove globals due to "exec" funcionality
      aliases = {k: v for k, v in aliases.items() if not (k.startswith('__') and k.endswith('__'))}

    print ("aliases = ", aliases)



    #
    # read configuration for the snapshot
    #
    variables_to_dump = []
    cut_to_dump = ''
    folder_where_to_save_trees = ''

    if os.path.exists(opt.snapshotFile) :
      handle = open(opt.snapshotFile,'r')
      exec(handle.read())
      handle.close()

    print ("variables_to_dump = ", variables_to_dump)
    print ("cut_to_dump = ", cut_to_dump)
    print ("folder_where_to_save_trees = ", folder_where_to_save_trees)


    # [MODIFIED] read batchScriptLocation and snapshotRunFolder from configuration.py if present
    batch_script_location = './scripts_snapshot_batch/'
    if hasattr(opt, 'batchScriptLocation'):
      batch_script_location = opt.batchScriptLocation

    scripts_run_folder = None
    if hasattr(opt, 'snapshotRunFolder'):
      scripts_run_folder = opt.snapshotRunFolder

    print ("batch_script_location = ", batch_script_location)
    print ("scripts_run_folder    = ", scripts_run_folder)


    factory = SnapshotFactory()
    factory.setValues(variables, cuts, supercut, samples, aliases,
                      variables_to_dump, cut_to_dump, folder_where_to_save_trees,
                      batch_script_location,
                      scripts_run_folder)
    factory.setConditions(opt.silentMode)

    # [MODIFIED] workflow aligned with mkShapes:
    #   (default)       -> makeScripts only (generate .cpp + Makefile, no compile)
    #   --localCompile  -> make -f MakefileSnapshot -j8
    #   --submitBatch   -> submit compiled executables to HTCondor (single submit file)
    #   --checkBatch    -> check job outputs and missing files

    if not opt.submitBatch and not opt.checkBatch and not opt.localCompile:
      factory.makeScripts()

    if opt.localCompile:
      factory.makeScripts()
      print("Compiling with: make -f MakefileSnapshot -j8")
      subprocess.run(["make", "-f", "MakefileSnapshot", "-j8"])

    if opt.submitBatch:
      factory.submitBatch()

    if opt.checkBatch:
      print ("Checking if the jobs finished succesfully")
      factory.checkBatch()



    print ("\n\n")
    print (" I'm done ... \n\n")