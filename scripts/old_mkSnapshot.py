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
                                        aliases_to_be_defined, functions_to_be_defined):

        # This is your C++ template as a Python string

        define_input_files_logic = ""

        define_input_files_logic += f"//  Nominal input files\n"

        define_input_files_logic += f'    auto* nominal = new TChain("{tree_name}");\n'

        for file_name in file_path:
          define_input_files_logic += f'    nominal->Add("{file_name}");\n'


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

        # not SafeDefine since these are new variables ... make it so!
        for aliasName, alias in functions_to_be_defined.items():
          if 'external' in alias.keys():
            with open(alias["external"], "r") as file: # this is the file with the c++ code to be included
              code_of_function_to_include += file.read()
          code_of_function_to_include += "\n"

          define_aliases_functions += "    varied_df = varied_df"

          if 'xmlfile' not in alias.keys(): # if "xmlfile" is present, it's a TMVA and it requires casting to float
            for i, var in enumerate(alias["variables"]):
              define_aliases_functions += f'\n                     .Define("_var_{aliasName}_{i}", "{var}")'
          # if "xmlfile" is present, it's a TMVA
          elif 'xmlfile' in alias.keys():
            for i, var in enumerate(alias["variables"]):
              define_aliases_functions += f'\n                     .Define("_var_{aliasName}_{i}_f", "(float){var}")'

          define_aliases_functions += ";\n"

          # if a simple c++ function
          if 'xmlfile' not in alias.keys():
            define_aliases_functions += f'\n'
            new_names = [f"_var_{aliasName}_{i}" for i in range(len(alias["variables"]))]
            sintax_variables = '{ "' + '", "'.join(new_names) + '" }'
            define_aliases_functions += f'    varied_df = varied_df.Define("{aliasName}", {alias["function"]}, {sintax_variables} );\n'
            define_aliases_functions += f'\n'
          # if "xmlfile" is present, it's a TMVA
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

        #
        # In RDataFrame, all variables must be defined before being plotted/used
        # The pro of this is that they can be used also in cuts
        #
        #
        # In case of already defined variables, for example "mll", and "mll" is already defined in the TTree,
        # the code "SafeDefine" should handle this, and "Define" a variable only when needed
        #

        booking_logic += f"//  Initial ...\n"
        booking_logic += f'    auto current_node = ROOT::RDF::RNode(varied_df);\n'


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
    def generate_makefile(self, file_paths, makefile_name):
        # Get ROOT configuration via shell calls
        cpp_flags = "$(shell root-config --cflags)"
        ld_flags = "$(shell root-config --libs) -lTMVA -lXMLIO"  # why TMVA is not by default?!?
        type_of_compilation = "-O2"
        # from "blabla.cpp" to "blabla"
        targets = [os.path.splitext(f)[0] for f in file_paths]

        with open(makefile_name, "w") as f:
            f.write("# Generated Makefile\n")
            f.write("CXX = g++\n")
            # f.write(f"CXXFLAGS = {type_of_compilation} -Wall {cpp_flags} -std=c++20 -Wcpp \n")  # c++20 used for "ends_with"
            f.write(f"CXXFLAGS = {type_of_compilation} -Wall {cpp_flags} \n")
            f.write(f"LDFLAGS = {ld_flags}\n\n")

            # 'all' target: the list of all executables to be created
            all_targets = " ".join(targets)
            f.write(f"all: {all_targets}\n\n")

            # The Pattern Rule
            # This says: To create ANY executable 'X', look for 'X.cpp'
            # It works across subdirectories automatically.
            f.write("%: %.cpp\n")
            f.write("\t$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)\n\n")

            # Clean target
            f.write("clean:\n")
            f.write(f"\trm -f {all_targets}\n")




    # _____________________________________________________________________________
    def setValues(self, variables, cuts, supercut, samples , aliases, variables_to_dump, cut_to_dump, folder_where_to_save_trees):

        self._variables = variables
        self._samples   = samples
        self._cuts      = cuts
        self._supercut  = supercut
        self._aliases   = aliases

        self._variables_to_dump           = variables_to_dump
        self._cut_to_dump                 = cut_to_dump
        self._folder_where_to_save_trees  = folder_where_to_save_trees


    # _____________________________________________________________________________
    def setConditions(self, silentMode):

        self._silentMode = silentMode


    # _____________________________________________________________________________
    def makeScripts(self):

        print ("=====================")
        print ("==== makeScripts ====")
        print ("=====================")

        os.system ("mkdir " + self._folder_where_to_save_trees + "/")

        os.system ("mkdir " + self._scripts_run_folder + "/")

        list_of_files_to_compile = []

        # print ("length of samples = ", len(self._samples))
        #
        # Loop over samples
        #
        for sampleName, sample in self._samples.items():
          os.system ("mkdir " + self._scripts_run_folder + "/" + sampleName + "/")
          #
          # create an executable for each file in each sample
          #
          # samples is a dictionary: {'subname' : [list of root files] }
          #

          global_weight = sample['weight'] if 'weight' in sample.keys() else "1."

          # weight per "sub-dataset"
          weights_subname = sample['weights'] if 'weights' in sample.keys() else {}


          #
          # check if "isData". If NOT isData then multiply by lumi otherwise no
          #
          if 'isData' in sample.keys() :
            if not (len(sample ['isData']) == 1 and sample ['isData'][0] == 'all') :
              # if you put 'all', all the root files are considered "data"
              # ok, it's data ... and so now?
              #
              pass
          else : # default is "scale to luminosity"
            global_weight = f"({global_weight}) * {self._lumi}"


          #
          # check if additional "Define" is needed, from "alises"
          #
          aliases_to_be_defined = {}
          functions_to_be_defined = {}
          for aliasName, alias in self._aliases.items():
            if 'samples' in alias.keys():
              if sampleName in alias['samples'] and 'expr' in alias.keys():
                aliases_to_be_defined[aliasName] = alias['expr']
              if sampleName in alias['samples'] and 'function' in alias.keys():
                functions_to_be_defined[aliasName] = alias


          for subname, list_root_files in sample['name'].items():
            print ("length of list_root_files = ", len(list_root_files))
            os.system ("mkdir " + self._scripts_run_folder + "/" + sampleName + "/" + subname + "/")

            # if weights per subname are listed, use them
            # e.g. per PD weights in data
            #
            if subname in weights_subname.keys():
              weight = f"({global_weight}) * ({weights_subname[subname]})"
            else :
              weight = global_weight

            #
            # Merge together different root files in one single job, not to have billions of jobs :)
            #
            make_job_every_N = 1
            if "FilesPerJob" in sample.keys():
              make_job_every_N = sample['FilesPerJob']
            #
            # create a list of lists of root files, one per job to be submitted
            # [ [root1, root2, root3], [root4, root5]]
            #
            chunks_list_root_files = [list_root_files[i:i + make_job_every_N] for i in range(0, len(list_root_files), make_job_every_N)]

            for i, root_files in enumerate(chunks_list_root_files):
              name_code = self._scripts_run_folder + "/" + sampleName + "/" + subname + "/" + "my_snapshot_" + sampleName + "_" + subname + "_" + str(i)
              tree = "Events"

              output_root_file_name = "nanoLatino_snapshot____" + sampleName + "_" + subname + "_" + str(i) + ".root"

              self.create_cpp_source_list_of_files(name_code, tree, root_files, output_root_file_name, sampleName, weight,
                                                   aliases_to_be_defined, functions_to_be_defined)

              list_of_files_to_compile.append(name_code)



        self.generate_makefile(list_of_files_to_compile, "MakefileSnapshot")

        # Run the compilation in parallel
        print("Start parallel compilation...")
        subprocess.run(["make", "-f", "MakefileSnapshot", "-j8"])
        #
        # subprocess.run() is a blocking call -> it will make the compilation to end before the next step
        #




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
                  if line_count > 1: # should I remove the warnings to have 0? FIXME
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

        os.system ("mkdir " + self._script_batch_location + "/")

        submission_dir = os.getcwd()

        #
        # Loop over samples
        #
        for sampleName, sample in self._samples.items():
          os.system ("mkdir " + self._script_batch_location + "/" + sampleName + "/")

          for subname, list_root_files in sample['name'].items():
            os.system ("mkdir " + self._script_batch_location + "/" + sampleName + "/" + subname + "/")

            make_job_every_N = 1
            if "FilesPerJob" in sample.keys():
              make_job_every_N = sample['FilesPerJob']
            chunks_list_root_files = [list_root_files[i:i + make_job_every_N] for i in range(0, len(list_root_files), make_job_every_N)]

            for i, root_files in enumerate(chunks_list_root_files):
              name_code = self._scripts_run_folder + "/" + sampleName + "/" + subname + "/" + "my_snapshot_" + sampleName + "_" + subname + "_" + str(i)
              name_code_no_folder = "my_snapshot_" + sampleName + "_" + subname + "_" + str(i)
              name_bash = self._script_batch_location + "/" + sampleName + "/" + subname + "/" + "my_script_" + str(i) + ".sh"
              output_root_file_name = "nanoLatino_snapshot____" + sampleName + "_" + subname + "_" + str(i) + ".root"

              bash_code = f"""#!/bin/bash
set -e  # Exit on error
echo "Job started at $(date)"
echo "Running on node $(hostname)"
echo "Current full path before running the snapshot: $(pwd)"
./{name_code_no_folder}
echo "Current full path after running the snapshot: $(pwd)"
mv {output_root_file_name}  {self._folder_where_to_save_trees}/
echo "Current full path after moving the files: $(pwd)"
echo "Current directory content after running:"
ls -lh
echo "Current full path: $(pwd)"

"""

#
# mkdir {self._folder_where_to_save_trees}
#

              with open(f"{name_bash}", "w") as f:
                f.write(bash_code)
                os.system ("chmod +x " + name_bash)
                # print ("name_bash = ", name_bash)

              name_submit = submission_dir + "/" + self._script_batch_location + "/" + sampleName + "/" + subname + "/" + "my_script_" + str(i) + ".sub"
              # name_folder = submission_dir + "/" + self._script_batch_location + "/" + sampleName + "/" + subname + "/"
              name_folder = self._script_batch_location + "/" + sampleName + "/" + subname + "/"
              name_folder_code = self._scripts_run_folder + "/" + sampleName + "/" + subname + "/"
              name_bash_no_folder = "my_script_" + str(i) + ".sh"

              if self._silentMode :
                output_file = "/dev/null"
                error_file  = "/dev/null"
                log_file    = "/dev/null"
              else :
                output_file = f"log/{name_bash_no_folder}.out"
                error_file  = f"log/{name_bash_no_folder}.err"
                log_file    = f"log/{name_bash_no_folder}.log"


              submit_code = f"""
initialdir            = {name_folder}
executable            = {name_bash}
transfer_input_files  = {submission_dir}/{name_folder_code}{name_code_no_folder}
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
output                = {output_file}
error                 = {error_file}
log                   = {log_file}
getenv                = True
+JobFlavour           = "longlunch"
queue
"""

              with open(f"{name_submit}", "w") as f:
                f.write(submit_code)

              print("Submitting job to HTCondor: condor_submit " + name_submit)
              subprocess.run(["condor_submit", f"{name_submit}"])





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
         '                                                                                                '
         --------------------------------------------------------------------------------------------------
         """

    print(header)

    parser = defaultParser()

    parser.add_argument("--submitBatch", action='store_true', dest="submitBatch",   help="Trigger the submission to lxbatch")
    parser.add_argument("--checkBatch",  action='store_true', dest="checkBatch",    help="Check if jobs are done succesfully")


    opt = parser.parse_args()
    print ("opt.pycfg            = ", opt.pycfg)
    loadDefaultOptions(parser, opt.pycfg)
    opt = parser.parse_args()


    print ("opt.submitBatch      = ", opt.submitBatch)
    print ("opt.checkBatch       = ", opt.checkBatch)


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




    factory = SnapshotFactory()
    factory.setValues( variables, cuts, supercut, samples , aliases, variables_to_dump, cut_to_dump, folder_where_to_save_trees)

    # factory._treeName  = opt.treeName
    # factory._energy    = opt.energy
    # factory._lumi      = opt.lumi
    # factory._tag       = opt.tag

    if not opt.submitBatch and not opt.checkBatch:
      factory.makeScripts()

    if opt.submitBatch :
      factory.submitBatch()

    if opt.checkBatch :
      print ("Checking if the jobs finished succesfully")
      factory.checkBatch()



    print ("\n\n")
    print (" I'm done ... \n\n")