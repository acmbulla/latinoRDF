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



# ----------------------------------------------------- ShapeFactory --------------------------------------

import logging
import subprocess


class ShapeFactory:
    _logger = logging.getLogger('ShapeFactory')

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

        self._outputDir = './test/'

        self._outputRootFile = 'myhistos.root'

        self._script_batch_location = './scripts_batch/'

        self._scripts_run_folder = "scripts_run"

        # conditions
        self._silentMode = False



    # _____________________________________________________________________________
    def __del__(self):
        pass


    # _____________________________________________________________________________
    #
    def create_library_cpp(self):


        cpp_code_library_header = f"""

#ifndef LIBRARY_UTILS_H
#define LIBRARY_UTILS_H

#include <string>
#include <vector>
#include "TH1.h"
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"

// Prototypes for your functions
void FoldHistogram(TH1* h, int doFold);
TH1D* UnrollHistogram(TH2D* h2);
ROOT::RDF::RNode SafeDefine(ROOT::RDF::RNode df, std::string var, std::string expr);
std::vector<std::string> getBranchNames(TTree* tree);

#endif // LIBRARY_UTILS_H

        """

        cpp_code_library = f"""

#include "library_utils.h"

//
// doFold
//   0 = no
//   1 = overflow
//   2 = underflow
//   3 = overflow and underflow
//

//                 TH1 works for TH1F, TH2F, ...
void FoldHistogram(TH1* h, int doFold) {{

  // 1D historgram
  if (h->GetDimension() == 1) {{

    // Get the number of visible bins
    int nBins = h->GetNbinsX();

    // overflow
    if (doFold == 1 || doFold == 3) {{
      double content_last     = h->GetBinContent(nBins);
      double error_last       = h->GetBinError(nBins);
      double content_overflow = h->GetBinContent(nBins+1);
      double error_overflow   = h->GetBinError(nBins+1);

      h->SetBinContent(nBins, content_last + content_overflow);
      h->SetBinError(nBins, std::sqrt(error_last*error_last + error_overflow*error_overflow));

      h->SetBinContent(nBins + 1, 0);
      h->SetBinError(nBins + 1, 0);

    }}

    // underflow
    if (doFold == 2 || doFold == 3) {{

      double content_first     = h->GetBinContent(1);
      double error_first       = h->GetBinError(1);
      double content_underflow = h->GetBinContent(0);
      double error_underflow   = h->GetBinError(0);

      h->SetBinContent(1, content_first + content_underflow);
      h->SetBinError(1, std::sqrt(error_first*error_first + error_underflow*error_underflow));

      h->SetBinContent(0, 0);
      h->SetBinError(0, 0);

    }}

  }}
  // 2D histogram (com'on you don't need a 3D version ...)
  else {{

    int nX = h->GetNbinsX();
    int nY = h->GetNbinsY();

    // overflow
    if (doFold == 1 || doFold == 3) {{
      for (int iy = 0; iy <= nY + 1; ++iy) {{
        h->SetBinContent(nX, iy, h->GetBinContent(nX, iy) + h->GetBinContent(nX + 1, iy));
        h->SetBinError(nX, iy, std::hypot(h->GetBinError(nX, iy), h->GetBinError(nX + 1, iy)));
        h->SetBinContent(nX + 1, iy, 0);
        h->SetBinError(nX + 1, iy, 0);
      }}
      for (int ix = 0; ix <= nX + 1; ++ix) {{
        h->SetBinContent(ix, nY, h->GetBinContent(ix, nY) + h->GetBinContent(ix, nY + 1));
        h->SetBinError(ix, nY, std::hypot(h->GetBinError(ix, nY), h->GetBinError(ix, nY + 1)));
        h->SetBinContent(ix, nY + 1, 0);
        h->SetBinError(ix, nY + 1, 0);
      }}
    }}

    // underflow
    if (doFold == 2 || doFold == 3) {{
      for (int iy = 0; iy <= nY + 1; ++iy) {{
        h->SetBinContent(1, iy, h->GetBinContent(1, iy) + h->GetBinContent(0, iy));
        h->SetBinError(1, iy, std::hypot(h->GetBinError(1, iy), h->GetBinError(0, iy)));
        h->SetBinContent(0, iy, 0);
        h->SetBinError(0, iy, 0);
      }}
      for (int ix = 0; ix <= nX + 1; ++ix) {{
        h->SetBinContent(ix, 1, h->GetBinContent(ix, 1) + h->GetBinContent(ix, 0));
        h->SetBinError(ix, 1, std::hypot(h->GetBinError(ix, 1), h->GetBinError(ix, 0)));
        h->SetBinContent(ix, 0, 0);
        h->SetBinError(ix, 0, 0);
      }}
    }}

  }}


}}


// --- transform 2D into 1D: unrolling
//
//      3    6    9
//      2    5    8
//      1    4    7
//


TH1D* UnrollHistogram(TH2D* h2) {{

  int nBinsX = h2->GetNbinsX();
  int nBinsY = h2->GetNbinsY();
  int totalBins = nBinsX * nBinsY;

  // Create a 1D histogram with the same name (plus a suffix)
  TString name = h2->GetName();
  TString title = h2->GetTitle();
  h2->SetName(name + "_old"); // Rename original to avoid conflict

  TH1D* h1 = new TH1D(name, title, totalBins, 0.5, totalBins + 0.5);
  h1->Sumw2();

  int i1D = 1;
  for (int iy = 1; iy <= nBinsY; ++iy) {{
      for (int ix = 1; ix <= nBinsX; ++ix) {{
          h1->SetBinContent(i1D, h2->GetBinContent(ix, iy));
          h1->SetBinError(i1D, h2->GetBinError(ix, iy));
          i1D++;
      }}
  }}

  return h1;
}}




#include "TMVA/RReader.hxx"
#include "TMVA/RInferenceUtils.hxx"



ROOT::RDF::RNode SafeDefine(ROOT::RDF::RNode df, std::string name, std::string expr) {{
    auto colNames = df.GetColumnNames();
    if (std::find(colNames.begin(), colNames.end(), name) == colNames.end()) {{
        return df.Define(name, expr);
    }}
    return df;
}}


std::vector<std::string> getBranchNames(TTree* tree) {{
  std::vector<std::string> names;
  TObjArray* branches = tree->GetListOfBranches();
  for (int i = 0; i < branches->GetEntries(); ++i)
    names.push_back(branches->At(i)->GetName());
  return names;
}}


        """



        with open(f"{self._scripts_run_folder}/library_utils.cpp", "w") as f:
            f.write(cpp_code_library)

        with open(f"{self._scripts_run_folder}/library_utils.h", "w") as f:
            f.write(cpp_code_library_header)



    # _____________________________________________________________________________
    #                                                              "file_path" is a list of root files!
    def create_cpp_source_list_of_files(self, output_name, tree_name, file_path, output_root_file_name, sampleName, weight,
                                        aliases_to_be_defined, functions_to_be_defined,
                                        subsamples, aliases_ordered=None):

        # This is your C++ template as a Python string

        define_input_files_logic = ""

        define_input_files_logic += f"//  Nominal input files\n"

        define_input_files_logic += f'    auto* nominal = new TChain("{tree_name}");\n'

        # define_input_files_logic += f'    nominal->Add("{file_path}");\n'

        for file_name in file_path:
          define_input_files_logic += f'    nominal->Add("{file_name}");\n'

        define_input_files_logic += f'    \n'
        define_input_files_logic += f'    auto nominalBranches = getBranchNames(nominal);\n'
        define_input_files_logic += f'    \n'

        #
        # now define the nuisances based on alternative trees
        #
        # 'kind': 'suffix',
        # 'type': 'shape',
        #
        # 'mapUp': 'ElepTup',
        # 'mapDown': 'ElepTdo',
        #
        #      'folderUp'    :   '/eos/cms/store/group/phys_higgs/cmshww/amassiro/HWWNano/Summer20UL18_106x_nAODv9_Full2018v9/MCl1loose2018v9__MCCorr2018v9NoJERInHorn__l2tightOR2018v9__ElepTup_suffix'
        #      'folderDown'   : ...
        #
        #

        #
        # extract name of root file nominal, without folder path
        #


        define_input_files_logic += f"//  Variations input files (if any)\n"

        for nuisanceName, nuisance in self._nuisances.items():

          if sampleName in nuisance['samples'].keys() and nuisance['type'] == 'shape' and ('kind' in nuisance.keys() and nuisance['kind'] == 'suffix'):

            # One TChain up and one TChain down, then add all the root files
            define_input_files_logic += f'''    auto* friend_{nuisance['mapUp']} = new TChain("{tree_name}");\n'''
            define_input_files_logic += f'''    auto* friend_{nuisance['mapDown']} = new TChain("{tree_name}");\n'''

            for file_name in file_path:

              # Find the index of the first "/" from the right
              index = file_name.rfind("/")
              # Extract everything from that index + 1 to the end
              # (We add 1 so we don't include the "/" itself)
              file_path_only_file_no_folder = file_name[index + 1:]

              define_input_files_logic += f'''    friend_{nuisance['mapUp']}->Add("{nuisance['folderUp']}/{file_path_only_file_no_folder}");\n'''
              define_input_files_logic += f'''    friend_{nuisance['mapDown']}->Add("{nuisance['folderDown']}/{file_path_only_file_no_folder}");\n'''

            # add friend TChain only once
            define_input_files_logic += f'''    nominal->AddFriend(friend_{nuisance['mapUp']}, "{nuisance['mapUp']}");\n'''
            define_input_files_logic += f'''    nominal->AddFriend(friend_{nuisance['mapDown']}, "{nuisance['mapDown']}");\n'''

            # varied branches are the same for all "added" trees
            define_input_files_logic += f'''    auto varBranches_{nuisance['mapUp']} = getBranchNames(friend_{nuisance['mapUp']});\n'''
            define_input_files_logic += f'''    auto varBranches_{nuisance['mapDown']} = getBranchNames(friend_{nuisance['mapDown']});\n'''


        define_input_files_logic += f'    ROOT::RDataFrame base_df(*nominal);\n'
        # define_input_files_logic += f'    ROOT::RDF::RInterface varied_df = base_df;\n'
        define_input_files_logic += f'    auto varied_df = ROOT::RDF::RNode(base_df);\n'


        #
        # define the weights needed for this specific sample
        #
        define_weights = ""
        define_weights += f"//  weigths needed\n"
        define_weights += f'    varied_df = SafeDefine(varied_df, "my_sample_weight", "{weight}");\n'


        #
        #
        # define the aliases and functions in the SAME ORDER as aliases.py
        # so that expr aliases depending on function results are defined after them
        #
        code_of_function_to_include = ""
        _included_external_files = set()  # track already-included files to avoid redefinitions

        lines_to_add     = ""   # extra #include / gSystem->Load() before main()
        lines_to_process = ""   # class instantiations inside main() before dataframe
        lines_to_global  = ""   # global scope declarations after external code, before main()
        _included_lines_to_add = set()  # avoid duplicates
        _included_lines_to_global = set()  # avoid duplicates

        define_aliases = ""
        define_aliases += f"//  aliases and functions (in aliases.py order)\n"

        define_aliases_functions = ""  # kept for backwards compat but will be empty

        _aliases_ordered = aliases_ordered if aliases_ordered is not None else {}
        # fallback: if no ordered dict, process expr first then functions (old behaviour)
        if not _aliases_ordered:
            for aliasName, expr in aliases_to_be_defined.items():
                _aliases_ordered[aliasName] = ('expr', expr)
            for aliasName, alias in functions_to_be_defined.items():
                _aliases_ordered[aliasName] = ('function', alias)

        for aliasName, (kind, payload) in _aliases_ordered.items():
            if kind == 'expr':
                # Replace single quotes with \" so strings inside the JIT expr
                # are correctly delimited in the generated C++ string literal.
                # Users can write 'M' in their expr instead of worrying about escaping.
                _payload = payload.replace("'", '\\"')
                define_aliases += f'    varied_df = SafeDefine(varied_df, "{aliasName}", "{_payload}");\n'
            else:
                alias = payload

                # --- linesToAdd: #include or gSystem->Load() before main() ---
                if 'linesToAdd' in alias.keys():
                    for line in alias['linesToAdd']:
                        if line not in _included_lines_to_add:
                            lines_to_add += line + "\n"
                            _included_lines_to_add.add(line)

                # --- linesToProcess: class instantiation inside main() ---
                if 'linesToProcess' in alias.keys():
                    for line in alias['linesToProcess']:
                        lines_to_process += "    " + line + "\n"

                # --- linesToGlobal: global scope after external code, before main() ---
                if 'linesToGlobal' in alias.keys():
                    for line in alias['linesToGlobal']:
                        if line not in _included_lines_to_global:
                            lines_to_global += line + "\n"
                            _included_lines_to_global.add(line)

                # --- external: include code before continue ---
                if 'external' in alias.keys():
                    if alias["external"] not in _included_external_files:
                        with open(alias["external"], "r") as file:
                            code_of_function_to_include += file.read()
                        code_of_function_to_include += "\n" 
                        _included_external_files.add(alias["external"])

                # --- expr: use directly as Define string (no variables loop) ---
                if ('linesToAdd' in alias.keys() or 'linesToGlobal' in alias.keys()) and 'expr' in alias.keys():
                    _expr = alias['expr'].replace("'", '\\"')
                    define_aliases += f'    varied_df = SafeDefine(varied_df, "{aliasName}", "{_expr}");\n'
                    continue

                define_aliases += "    varied_df = varied_df"
                if 'xmlfile' not in alias.keys():
                    for i, var in enumerate(alias["variables"]):
                        define_aliases += f'\n                     .Define("_var_{aliasName}_{i}", "{var}")'
                else:
                    for i, var in enumerate(alias["variables"]):
                        define_aliases += f'\n                     .Define("_var_{aliasName}_{i}_f", "(float){var}")'
                define_aliases += ";\n"

                if 'xmlfile' not in alias.keys():
                    define_aliases += f'\n'
                    new_names = [f"_var_{aliasName}_{i}" for i in range(len(alias["variables"]))]
                    sintax_variables = '{ "' + '", "'.join(new_names) + '" }'
                    define_aliases += f'    varied_df = varied_df.Define("{aliasName}", {alias["function"]}, {sintax_variables} );\n'
                    define_aliases += f'\n'
                else:
                    define_aliases += f'\n'
                    submission_dir = os.getcwd()
                    define_aliases += f'''    TMVA::Experimental::RReader my_model_{aliasName}("{submission_dir}/{alias['xmlfile']}");\n'''
                    new_names = [f"_var_{aliasName}_{i}_f" for i in range(len(alias["variables"]))]
                    sintax_variables = '{ "' + '", "'.join(new_names) + '" }'
                    define_aliases += f'''    varied_df = varied_df.Define("{aliasName}",
                     TMVA::Experimental::Compute<{len(alias['variables'])}, float>(my_model_{aliasName}),
                     {sintax_variables});\n'''
                    define_aliases += f'\n'


        #
        # register the variations
        #
        #  some aliases do be defined AFTER Vary : e.g. BDT after vary of "scale of lepton"
        #     register_variations_logic_for_alternative_trees --> variations defined before aliases
        #  some aliases do be defined BEFORE Vary: define of "weights"
        #     register_variations_logic --> variations defined after aliases
        #

        register_variations_logic_for_alternative_trees = ""
        register_variations_logic_for_alternative_trees += f"//  Register the variations for alternative trees. They must defined as soon as possible\n"
        register_variations_logic_for_alternative_trees = "\n"

        register_maps_logic_for_alternative_trees_on_the_fly_calculation = ""

        register_variations_for_special_nuisances_on_the_fly = ""

        code_of_functions_and_classes_to_include_for_special_nuisances = ""

        register_variations_logic = ""
        register_variations_logic += f"//  Register the variations\n"

        if len(self._nuisances) != 0 :
          first_time_suffix = 0
          for nuisanceName, nuisance in self._nuisances.items():


            #
            # 'kind': 'suffixFly'
            # 'name': 'CMS_res_j_2018_fly'
            # 'pattern' : 'CMS_my_res_j_XXXX'
            # 'numVariations' : 3,
            # 'external' : "code/jer_on_the_fly.c",
            # 'code_in_main' : "code/jer_on_the_fly_main.c",
            #
            # Nuisance based on alternative tree, calculated on the fly
            # and whose effect is calculated as the "averaged" effect, on numVariations times
            # E.g. JER -> affected by the choice of the seed
            # Purpose: remove unwanted fluctuations
            #
            # 'external' is the code where the class/function is defined
            # 'code_in_main' is the code that will be copied in the main, e.g. the class creation
            #

            if sampleName in nuisance['samples'].keys() and nuisance['type'] == 'shape' and ('kind' in nuisance.keys() and nuisance['kind'] == 'suffixFly'):

              register_maps_logic_for_alternative_trees_on_the_fly_calculation += f'''    map_nuisances_name_for_averaged_histograms["{nuisance['pattern']}"] = std::pair<std::string, int>("{nuisance['name']}", {nuisance['numVariations']});\n'''

              with open(nuisance["external"], "r") as file:
                code_of_functions_and_classes_to_include_for_special_nuisances += file.read()
                code_of_functions_and_classes_to_include_for_special_nuisances += "\n"

              with open(nuisance["code_in_main"], "r") as file:
                register_variations_for_special_nuisances_on_the_fly += file.read()
                register_variations_for_special_nuisances_on_the_fly += "\n"


            #
            # "shape" : "suffix" nuisances
            #  a.k.a. alternative friend trees
            #  also the use of a different weight is possible, to be activated with "weigthsPerSample"    FIXME  Not yet implemented ... was this ever used?
            #         NB: suffix AND weight at the same time is yet NOT implemented (is it possible with RDataFrame?) FIXME
            #
            #  e.g.
            #    'mapUp': 'ElepTup',
            #    'mapDown': 'ElepTdo',
            #    'samples': dict((skey, ['1', '1']) for skey in mcALL),
            #    'folderUp': makeMCDirectory('ElepTup_suffix'),
            #    'folderDown': makeMCDirectory('ElepTdo_suffix'),
            #    'weigthsPerSample': true
            #
            #
            if sampleName in nuisance['samples'].keys() and nuisance['type'] == 'shape' and ('kind' in nuisance.keys() and nuisance['kind'] == 'suffix'):

              if first_time_suffix == 0 :
                register_variations_logic_for_alternative_trees += f'    int suffix_size = 0;\n'
                first_time_suffix = 1

              register_variations_logic_for_alternative_trees += f'''    suffix_size = {len(nuisance['mapUp']) + 1};\n'''
              register_variations_logic_for_alternative_trees += f'''    for (const auto& branch : varBranches_{nuisance['mapUp']}) {{\n'''
              register_variations_logic_for_alternative_trees += f'''      if (int(branch.size()) >= (suffix_size-1) && branch.compare(branch.size() - (suffix_size-1), (suffix_size-1), "{nuisance['mapUp']}") == 0 ) {{\n'''
              register_variations_logic_for_alternative_trees += f'''        std::string nomCol = branch.substr(0, branch.size() - suffix_size);\n'''
              #
              # check if the varied column exists also in the nominal ttree
              # If it does not exist, skip. But ... how is it possible you have the varied value but not the nominal one? How did you produce the varied variable?
              #
              register_variations_logic_for_alternative_trees += f'''        if (std::find(nominalBranches.begin(), nominalBranches.end(), nomCol) != nominalBranches.end()) {{\n'''
              register_variations_logic_for_alternative_trees += f'''          std::string col_type = varied_df.GetColumnType(nomCol);\n'''
              register_variations_logic_for_alternative_trees += f'''          std::string expression = (col_type.find("RVec") != std::string::npos)\n'''
              register_variations_logic_for_alternative_trees += f'''            ? ("ROOT::RVec<" + col_type + ">{{" + branch + "}}")\n'''
              register_variations_logic_for_alternative_trees += f'''            : ("ROOT::RVec<" + col_type + ">{{static_cast<" + col_type + ">(" + branch + ")}}");\n'''
              # why do I have long float vs float between the nominal and the varied variable?
              # FIXME: some fix might be needed in the post-processing step --> removing the (technically not-needed) casting might make this code faster
              #
              register_variations_logic_for_alternative_trees += f'''          varied_df = varied_df.Vary(\n'''
              register_variations_logic_for_alternative_trees += f'''                                    nomCol,\n'''
              register_variations_logic_for_alternative_trees += f'''                                    expression,\n'''
              register_variations_logic_for_alternative_trees += f'''                                    {{"up"}},\n'''
              register_variations_logic_for_alternative_trees += f'''                                    "{nuisance['name']}"\n'''
              register_variations_logic_for_alternative_trees += f'''                                    );\n'''
              register_variations_logic_for_alternative_trees += f'''        }};\n'''
              register_variations_logic_for_alternative_trees += f'''      }};\n'''
              register_variations_logic_for_alternative_trees += f'''    }};\n'''

              register_variations_logic_for_alternative_trees += f'''    \n'''
              register_variations_logic_for_alternative_trees += f'''    suffix_size = {len(nuisance['mapDown']) + 1};\n'''
              register_variations_logic_for_alternative_trees += f'''    for (const auto& branch : varBranches_{nuisance['mapDown']}) {{\n'''
              register_variations_logic_for_alternative_trees += f'''      if (int(branch.size()) >= (suffix_size-1) && branch.compare(branch.size() - (suffix_size-1), (suffix_size-1), "{nuisance['mapDown']}") == 0 ) {{\n'''
              register_variations_logic_for_alternative_trees += f'''        std::string nomCol = branch.substr(0, branch.size() - suffix_size);\n'''
              register_variations_logic_for_alternative_trees += f'''        if (std::find(nominalBranches.begin(), nominalBranches.end(), nomCol) != nominalBranches.end()) {{\n'''
              register_variations_logic_for_alternative_trees += f'''          std::string col_type = varied_df.GetColumnType(nomCol);\n'''
              register_variations_logic_for_alternative_trees += f'''          std::string expression = (col_type.find("RVec") != std::string::npos)\n'''
              register_variations_logic_for_alternative_trees += f'''            ? ("ROOT::RVec<" + col_type + ">{{" + branch + "}}")\n'''
              register_variations_logic_for_alternative_trees += f'''            : ("ROOT::RVec<" + col_type + ">{{static_cast<" + col_type + ">(" + branch + ")}}");\n'''
              register_variations_logic_for_alternative_trees += f'''          varied_df = varied_df.Vary(\n'''
              register_variations_logic_for_alternative_trees += f'''                                    nomCol,\n'''
              register_variations_logic_for_alternative_trees += f'''                                    expression,\n'''
              register_variations_logic_for_alternative_trees += f'''                                    {{"do"}},\n'''
              register_variations_logic_for_alternative_trees += f'''                                    "{nuisance['name']}"\n'''
              register_variations_logic_for_alternative_trees += f'''                                    );\n'''
              register_variations_logic_for_alternative_trees += f'''        }};\n'''
              register_variations_logic_for_alternative_trees += f'''      }};\n'''
              register_variations_logic_for_alternative_trees += f'''    }};\n'''

            #
            # "shape" : "weight" nuisances
            #  a.k.a. same tree but with a different weight
            #
            #  e.g. nominal,         up,                      down
            #    ['SFweightMu', 'SFweightMuUp',         'SFweightMuDown']
            #    ['LepWPSF',    'LepWPSF*SFweightMuUp', 'LepWPSF*SFweightMuDown']
            #  or
            #         up,             down
            #    ['SFweightMuUp', 'SFweightMuDown']
            #     --> in this second definition the system assumes this is a multiplicative effect, namely the nominal is a "1."
            #         and it's always applied, since there is no check that the "nominal" exists
            #         while above the system first checks that the nominal exists before applying the variation
            #
            #
            if sampleName in nuisance['samples'].keys() and nuisance['type'] == 'shape' and ('kind' in nuisance.keys() and nuisance['kind'] == 'weight'):

              if len(nuisance['samples'][sampleName]) == 3:
                #
                # check if '{nuisance['samples'][sampleName][0]}' is in the string 'weight' (that is actually defined in 'my_sample_weight')
                # if yes, propagate the weight variation

                if nuisance['samples'][sampleName][0] in weight:

                  variation_up   =  weight.replace(nuisance['samples'][sampleName][0], nuisance['samples'][sampleName][1])
                  variation_down =  weight.replace(nuisance['samples'][sampleName][0], nuisance['samples'][sampleName][2])

                  register_variations_logic += f'''    varied_df = varied_df.Vary(\n'''
                  register_variations_logic += f'''      "my_sample_weight",\n'''
                  register_variations_logic += f'''      "ROOT::RVecD{{{variation_up},{variation_down}}}",\n'''
                  register_variations_logic += f'''      {{"up", "do"}},\n'''
                  register_variations_logic += f'''      "{nuisance['name']}"\n'''
                  register_variations_logic += f'''      );\n'''

              #
              # otherwise it means that is a multiplicative weight, meaning {original weight} * up/down weight
              # If such, I can exploit "my_sample_weight" and not use the extended definition
              #

              else :

                  # variation_up     =  f"({weight}) * {nuisance['samples'][sampleName][0]}"
                  # variation_down   =  f"({weight}) * {nuisance['samples'][sampleName][1]}"
                  variation_up     =  f"(my_sample_weight) * {nuisance['samples'][sampleName][0]}"
                  variation_down   =  f"(my_sample_weight) * {nuisance['samples'][sampleName][1]}"

                  register_variations_logic += f'''    varied_df = varied_df.Vary(\n'''
                  register_variations_logic += f'''      "my_sample_weight",\n'''
                  register_variations_logic += f'''      "ROOT::RVecD{{{variation_up},{variation_down}}}",\n'''
                  register_variations_logic += f'''      {{"up", "do"}},\n'''
                  register_variations_logic += f'''      "{nuisance['name']}"\n'''
                  register_variations_logic += f'''      );\n'''


        booking_logic = ""


        #
        # In RDataFrame, all variables must be defined before being plotted
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



        #
        # once all variables are defined, they can be used and plotted with "variableName"
        #
        # In root file:    <cut>/<variable>/histo_<sample>
        #
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
              list_cuts[ cutName ] = "(" + expression + ")"   # add the un-categorized phase space too! it comes "for free"
              for category_name, category in categories.items():
                # list_cuts[ cutName + "_" + category_name] = "(" + expression + ") && (" + category + ")"
                # this above is not needed, as handled by the branching of the nodes
                list_cuts[ cutName + "_" + category_name] = "(" + category + ")"
            else :
              list_cuts[ cutName ] = expression
          else :
            list_cuts[ cutName ] = cut


          mother_cut_name = ""
          for icut, (this_cutName, this_cut) in enumerate(list_cuts.items()):
            if len(list_cuts) > 1 :
              if icut == 0:
                define_variables_logic += f'    auto node_{this_cutName} = current_node.Filter("{this_cut}", "{this_cutName}");\n'
                mother_cut_name = this_cutName
              else:
                # with this I'm nesting the node into the "mother node" -> it's RDataFrame power
                define_variables_logic += f'    auto node_{this_cutName} = node_{mother_cut_name}.Filter("{this_cut}", "{this_cutName}");\n'
            else :
              define_variables_logic += f'    auto node_{this_cutName} = current_node.Filter("{this_cut}", "{this_cutName}");\n'

            #
            # if subsamples are defined for this specific sample
            #
            #
            #   'subsamples': {
            #     #            definition          specific weight
            #     'Low' : { 'cut' : 'mll<40',    'weight' :  '1.34'},
            #     'High': { 'cut' : 'mll>40',    'weight' :  '1.00'}
            #   }
            #
            # subsamples can be useful for:
            #  - EFT samples
            #  - unfolding
            #
            # define the "cuts" nodes for the subsamples
            #
            if subsamples :
              for sub_name, sub_cut_name in subsamples.items():
                if 'cut' in sub_cut_name.keys():
                  define_variables_logic += f'''    ROOT::RDF::RNode node_{this_cutName}___{sampleName}_{sub_name} = node_{this_cutName}.Filter("{sub_cut_name['cut']}", "cut_{sampleName}_{sub_name}");\n'''
                  if 'weight' in sub_cut_name.keys():
                    define_variables_logic += f'    node_{this_cutName}___{sampleName}_{sub_name} = SafeDefine(node_{this_cutName}___{sampleName}_{sub_name}, "my_sample_weight_{sampleName}_{sub_name}", "(my_sample_weight)*{sub_cut_name['weight']}");\n'
                elif 'weight' in sub_cut_name.keys():
                    define_variables_logic += f'    node_{this_cutName} = SafeDefine(node_{this_cutName}, "my_sample_weight_{sampleName}_{sub_name}", "(my_sample_weight)*{sub_cut_name['weight']}");\n'




            for variableName, variable in self._variables.items():
                #  Different options for range definition:
                #
                #  1D histograms:
                #     'range' : (200,10,500),
                #     'range' : ([12, 17, 25, 30, 35, 40, 45, 65, 200]),
                #
                #  2D histograms:
                #     'range' : (5, 0.0, 1.0, 10,  0., 10000.),
                #     'range' : ([12, 17, 25, 30, 35, 40, 45, 65, 200],[60, 95, 110, 135, 200],),
                #

                variable_range = variable['range']
                model_str = ""

                if isinstance(variable_range, list):
                  # 1D
                  #     'range' : ([12, 17, 25, 30, 35, 40, 45, 65, 200]),
                  #
                  bins = len(variable_range) - 1
                  # Convert Python list [12, 17...] to C++ string "12, 17, ..."
                  edges_str = ", ".join(map(str, variable_range))
                  model_str = f'{bins}, (const double[]){{{edges_str}}}'
                elif isinstance(variable_range, tuple):
                  if len(variable_range) == 3:
                    # 1D
                    #     'range' : (200,10,500),
                    #
                    bins, v_min, v_max = variable_range
                    model_str = f'{bins}, {v_min}, {v_max}'

                  elif len(variable_range) == 6:
                    # 2D
                    #     'range' : (5, 0.0, 1.0, 10,  0., 10000.),
                    #
                    model_str = ", ".join(map(str, variable_range))

                  elif len(variable_range) == 2 and isinstance(variable_range[0], list):
                    # 2D
                    #     'range' : ([12, 17, 25, 30, 35, 40, 45, 65, 200],[60, 95, 110, 135, 200],),
                    #
                    nx, ny = len(variable_range[0]) - 1, len(variable_range[1]) - 1
                    ex = ", ".join(map(str, variable_range[0]))
                    ey = ", ".join(map(str, variable_range[1]))
                    model_str = f'{nx}, (const double[]){{{ex}}}, {ny}, (const double[]){{{ey}}}'

                #
                # different booking depending if it is a 1D or 2D histogram
                #

                if 'is2d' in variable.keys() and variable['is2d'] == 1:
                  # print ("variable['name'].split(':')", variable['name'].split(':'))
                  # Split "varX:varY" into separate columns
                  v_x, v_y = variable['name'].split(':')
                  histo_call = f'Histo2D({{"h_{variableName}", "{variableName}", {model_str}}}, "{v_x}", "{v_y}", "my_sample_weight")'

                  # We add each RResultPtr to a vector called 'histograms'
                  define_variables_logic += f'    hist_map_2D["{this_cutName}"].push_back(node_{this_cutName}.{histo_call});\n'

                else:
                  # histo_call = f'Histo1D({{"h_{variableName}", "{variableName}", {model_str}}}, "{variableName}", "my_sample_weight")'
                  histo_call = f'Histo1D<float>({{"h_{variableName}", "{variableName}", {model_str}}}, "{variableName}", "my_sample_weight")'
                  # Possible optimization explicitly saying <float>?
                  # FIXME : maybe trigger "int" when needed or it's ok since most of the times it's a float and at most you cast an int into a float?



                  # We add each RResultPtr to a vector called 'histograms'
                  define_variables_logic += f'    hist_map_1D["{this_cutName}"].push_back(node_{this_cutName}.{histo_call});\n'

                #
                # if subsamples then there are new "cuts available"
                #
                if subsamples :
                  for sub_name, sub_cut_name in subsamples.items():

                    name_weight_for_subsample = "my_sample_weight"
                    if 'weight' in sub_cut_name.keys():
                      name_weight_for_subsample = f"my_sample_weight_{sampleName}_{sub_name}"

                    if 'cut' in sub_cut_name.keys():
                      if 'is2d' in variable.keys() and variable['is2d'] == 1:
                        v_x, v_y = variable['name'].split(':')
                        histo_call = f'Histo2D({{"h_{variableName}", "{variableName}", {model_str}}}, "{v_x}", "{v_y}", "{name_weight_for_subsample}")'
                        define_variables_logic += f'    hist_map_2D["{this_cutName}___{sampleName}_{sub_name}"].push_back(node_{this_cutName}___{sampleName}_{sub_name}.{histo_call});\n'
                      else:
                        histo_call = f'Histo1D<float>({{"h_{variableName}", "{variableName}", {model_str}}}, "{variableName}", "{name_weight_for_subsample}")'
                        define_variables_logic += f'    hist_map_1D["{this_cutName}___{sampleName}_{sub_name}"].push_back(node_{this_cutName}___{sampleName}_{sub_name}.{histo_call});\n'




        define_variables_logic += f'    \n'
        define_variables_logic += f'    std::vector<std::string> list_of_variables_1D;\n'
        define_variables_logic += f'    std::vector<std::string> list_of_variables_2D;\n'
        define_variables_logic += f'    std::vector<int> list_of_variables_fold_1D;\n'
        define_variables_logic += f'    std::vector<int> list_of_variables_fold_2D;\n'

        for variableName, variable in self._variables.items():
          if 'is2d' in variable.keys() and variable['is2d'] == 1:
            #  2D::   vary:varx
            define_variables_logic += f'    list_of_variables_2D.push_back("{variableName}");\n'
            if 'fold' in variable.keys():
              define_variables_logic += f'    list_of_variables_fold_2D.push_back({variable["fold"]});\n'
            else :
              define_variables_logic += f'    list_of_variables_fold_2D.push_back(0);\n'
          else :
            #  1D::   var
            define_variables_logic += f'    list_of_variables_1D.push_back("{variableName}");\n'
            if 'fold' in variable.keys():
              define_variables_logic += f'    list_of_variables_fold_1D.push_back({variable["fold"]});\n'
            else :
              define_variables_logic += f'    list_of_variables_fold_1D.push_back(0);\n'


        cpp_code = f"""

#include "library_utils.h"

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RVec.hxx"

#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TH1D.h"
#include "TDirectory.h"
#include "TChain.h"
#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <typeinfo>
#include <regex>


#include "TMVA/RReader.hxx"
#include "TMVA/RInferenceUtils.hxx"


// --- Automatically generated: linesToAdd (includes and library loads) ---

{lines_to_add}

// ----------------------------------------

// --- Automatically generated: code to be added for additional functions ---

{code_of_function_to_include}

// ----------------------------------------

// --- Automatically generated: global scope declarations ---

{lines_to_global}

// ----------------------------------------

// --- Automatically generated: code to be added for special nuisances ---

{code_of_functions_and_classes_to_include_for_special_nuisances}

// ----------------------------------------


int main() {{
    ROOT::DisableImplicitMT();
    std::cout << "   [DEBUG] ImplicitMT disabled" << std::endl;

    // --- Automatically generated: linesToProcess (class instantiations) ---

    {lines_to_process}
    std::cout << "   [DEBUG] linesToProcess done" << std::endl;

    // ----------------------------------------

    // --- Automatically generated input root files ---

    {define_input_files_logic}
    std::cout << "   [DEBUG] TChain created" << std::endl;

    // ----------------------------------------

    std::map<std::string, std::vector<ROOT::RDF::RResultPtr<TH1D>>> hist_map_1D;
    std::map<std::string, std::vector<ROOT::RDF::RResultPtr<TH2D>>> hist_map_2D;

    // --- Automatically generated register nuisances variations for alternative trees ---

    {register_variations_logic_for_alternative_trees}



    // --- Automatically generated register nuisances variations for alternative-like trees, as they are calculated on the fly ---

    {register_variations_for_special_nuisances_on_the_fly}

    // ----------------------------------------

    // --- Automatically generated define aliases for functions, no JIT ---
    // (must come before expr aliases that depend on them, e.g. WZ_vars before mll_Z)

    {define_aliases_functions}
    std::cout << "   [DEBUG] define_aliases_functions done" << std::endl;

    // ----------------------------------------

    // --- Automatically generated define aliases ---

    {define_aliases}
    std::cout << "   [DEBUG] define_aliases done" << std::endl;

    // ----------------------------------------

    // --- Automatically generated define weights ---

    {define_weights}
    std::cout << "   [DEBUG] define_weights done" << std::endl;

    // ----------------------------------------


    // --- Automatically generated register nuisances variations ---

    {register_variations_logic}
    std::cout << "   [DEBUG] register_variations done" << std::endl;

    // ----------------------------------------

    // --- Automatically generated ... just the main node ---

    {booking_logic}
    std::cout << "   [DEBUG] booking done" << std::endl;

    //
    // from now on: current_node
    // ----------------------------------------

    // --- Automatically generated define cuts and histograms ---

    {define_variables_logic}
    std::cout << "   [DEBUG] define_variables done" << std::endl;

    // ----------------------------------------

    //
    // Define variations
    //

    std::vector<ROOT::RDF::Experimental::RResultMap<TH1D>> results_1D_variations;
    std::vector<ROOT::RDF::Experimental::RResultMap<TH2D>> results_2D_variations;

    std::vector<ROOT::RDF::RResultHandle> all_booking_1D;
    std::vector<ROOT::RDF::RResultHandle> all_booking_2D;
    // FIXME: can these be merged?

    for (auto& [cut_label, h_list] : hist_map_1D) {{
      for (auto& h : h_list) {{
        results_1D_variations.push_back( ROOT::RDF::Experimental::VariationsFor(h) );
        all_booking_1D.push_back(h);
      }}
    }}

    for (auto& [cut_label, h_list] : hist_map_2D) {{
      for (auto& h : h_list) {{
        results_2D_variations.push_back( ROOT::RDF::Experimental::VariationsFor(h) );
        all_booking_2D.push_back(h);
      }}
    }}

    //
    // Now and RunGraphs --> the actual loop on the events
    //

    std::cout << "   ROOT::RDF::RunGraphs 1D" << std::endl;
    if (!hist_map_1D.empty()) {{
      ROOT::RDF::RunGraphs(all_booking_1D);
    }}
    std::cout << "   [DEBUG] RunGraphs 1D done" << std::endl;

    std::cout << "   ROOT::RDF::RunGraphs 2D" << std::endl;
    if (!hist_map_2D.empty()) {{
      ROOT::RDF::RunGraphs(all_booking_2D);
    }}
    std::cout << "   [DEBUG] RunGraphs 2D done" << std::endl;


    std::cout << "   Done, now save in the final root file" << std::endl;


    //
    // Prepare to average histograms for specific nuisances
    //
    //         pattern    nuisance-name
    std::map<std::string, std::pair<std::string, int> > map_nuisances_name_for_averaged_histograms;

    {register_maps_logic_for_alternative_trees_on_the_fly_calculation}

    //
    // In root file:    <cut>/<variable>/histo_<sample>
    //

    TFile out_file("{self._outputDir}/{output_root_file_name}", "RECREATE");
    int big_loop = 0;
    for (auto& [cut_label, h_list] : hist_map_1D) {{

      std::string current_cut = cut_label;
      // check if subsamples case: this_cutName___sampleName_sub_name
      std::string sub_name = "";
      std::string delimiter = "___";
      size_t pos_triple = cut_label.find(delimiter);
      if (pos_triple != std::string::npos) {{
        std::string this_cutName = cut_label.substr(0, pos_triple);
        std::string rest = cut_label.substr(pos_triple + 3);
        std::string sampleName = "{sampleName}";
        size_t pos = rest.find(sampleName);
        sub_name = rest.substr(pos + sampleName.length() + 1);
        current_cut = this_cutName;
      }}

      // Create a folder for this cut
      TDirectory *dir = out_file.GetDirectory(current_cut.c_str());
      if (!dir) {{
        // If it doesn't exist, create it
        dir = out_file.mkdir(current_cut.c_str());
      }}

      for (size_t ivar = 0; ivar < h_list.size(); ivar++) {{
        dir->cd();
        TDirectory *subdir = out_file.GetDirectory( (current_cut+"/"+list_of_variables_1D.at(ivar)).c_str() );
        if (!subdir) {{
          // If it doesn't exist, create it
          subdir = out_file.mkdir( (current_cut+"/"+list_of_variables_1D.at(ivar)).c_str() );
        }}
        subdir->cd();

        // get the nominal and the variations too
        auto all_histos = results_1D_variations.at(big_loop);
        big_loop++;

        std::map<std::string, TH1D> map_nuisances_average_histogram;

        for (auto& [name, histo] : all_histos) {{
          std::string temp_name;
          if (name == "nominal") {{
            if (sub_name == "") temp_name = "histo_{sampleName}";
            else                temp_name = "histo_{sampleName}_" + sub_name;
          }}
          else {{
            temp_name = name.c_str();
            // scale_e_2018_UL:up --> scale_e_2018_ULup
            temp_name.erase(std::remove(temp_name.begin(), temp_name.end(), ':'), temp_name.end());
            if (sub_name == "") temp_name = ("histo_{sampleName}_" + temp_name);
            else                temp_name = ("histo_{sampleName}_" + sub_name + "_" + temp_name);
          }}
          gDirectory = subdir;
          histo->SetName(temp_name.c_str());

          if (list_of_variables_fold_1D.at(ivar) != 0) FoldHistogram(histo.get(), list_of_variables_fold_1D.at(ivar));

          bool found_any_nuisances_for_averaged_histograms = false;
          for (auto const& [string_pattern, nuisance_name_numVariations] : map_nuisances_name_for_averaged_histograms) {{

            std::string nuisance_name = nuisance_name_numVariations.first;
            int numVariations_here = nuisance_name_numVariations.second;

            std::string string_pattern_up = "histo_.*__NORM__.*_" + string_pattern + "up";
            std::string string_pattern_do = "histo_.*__NORM__.*_" + string_pattern + "do";

            std::regex pattern_up(string_pattern_up);
            std::regex pattern_do(string_pattern_do);

            if (std::regex_search(temp_name, pattern_up) || std::regex_search(temp_name, pattern_do)) {{
              std::regex pattern_to_remove("__NORM__[0-9]+_");
              std::string histo_temp_name = std::regex_replace(temp_name, pattern_to_remove, "");
              std::regex pattern(string_pattern);

              histo_temp_name = std::regex_replace(histo_temp_name, pattern, nuisance_name);

              histo->SetName(histo_temp_name.c_str());

              if (std::regex_search(temp_name, pattern_up)) {{
                if (auto it = map_nuisances_average_histogram.find(string_pattern + "up"); it != map_nuisances_average_histogram.end()) {{
                  it->second.Add(histo.get(), 1. / numVariations_here);
                }}
                else {{
                  map_nuisances_average_histogram[nuisance_name + "up"] = (*histo);
                  found_any_nuisances_for_averaged_histograms = true;
                  break;
                }}
              }}
              else {{
                if (auto it = map_nuisances_average_histogram.find(string_pattern + "do"); it != map_nuisances_average_histogram.end()) {{
                  it->second.Add(histo.get(), 1. / numVariations_here);
                }}
                else {{
                  map_nuisances_average_histogram[nuisance_name + "do"] = (*histo);
                  found_any_nuisances_for_averaged_histograms = true;
                  break;
                }}
              }}
            }}
          }}
          if (not found_any_nuisances_for_averaged_histograms) {{
            histo->Write();
          }}

          for (auto const& [name, temp_histo] : map_nuisances_average_histogram) {{
            temp_histo.Write();
          }}

        }}
      }}
      out_file.cd();
    }}

    big_loop = 0;
    for (auto& [cut_label, h_list] : hist_map_2D) {{

      std::string current_cut = cut_label;
      std::string sub_name = "";
      std::string delimiter = "___";
      size_t pos_triple = cut_label.find(delimiter);
      if (pos_triple != std::string::npos) {{
        std::string this_cutName = cut_label.substr(0, pos_triple);
        std::string rest = cut_label.substr(pos_triple + 3);
        std::string sampleName = "{sampleName}";
        size_t pos = cut_label.find(sampleName);
        sub_name = cut_label.substr(pos + sampleName.length() + 1);
        current_cut = this_cutName;
      }}

      TDirectory *dir = out_file.GetDirectory(current_cut.c_str());
      if (!dir) {{
        dir = out_file.mkdir(current_cut.c_str());
      }}
      for (size_t ivar = 0; ivar < h_list.size(); ivar++) {{
        dir->cd();
        TDirectory *subdir = out_file.GetDirectory( (current_cut+"/"+list_of_variables_2D.at(ivar)).c_str() );
        if (!subdir) {{
          subdir = out_file.mkdir( (current_cut+"/"+list_of_variables_2D.at(ivar)).c_str() );
        }}
        subdir->cd();

        auto all_histos = results_2D_variations.at(big_loop);
        big_loop++;

        std::map<std::string, TH1D> map_nuisances_average_histogram;

        for (auto& [name, histo] : all_histos) {{
          std::string temp_name;
          if (name == "nominal") {{
            if (sub_name == "") temp_name = "histo_{sampleName}";
            else                temp_name = "histo_{sampleName}_" + sub_name;
          }}
          else {{
            temp_name = name.c_str();
            temp_name.erase(std::remove(temp_name.begin(), temp_name.end(), ':'), temp_name.end());
            if (sub_name == "") temp_name = ("histo_{sampleName}_" + temp_name);
            else                temp_name = ("histo_{sampleName}_" + sub_name + "_" + temp_name);
          }}
          gDirectory = subdir;
          histo->SetName(temp_name.c_str());

          if (list_of_variables_fold_2D.at(ivar) != 0) FoldHistogram(histo.get(), list_of_variables_fold_2D.at(ivar));

          auto new_histo = UnrollHistogram(dynamic_cast<TH2D*>(histo.get()));

          bool found_any_nuisances_for_averaged_histograms = false;
          for (auto const& [string_pattern, nuisance_name_numVariations] : map_nuisances_name_for_averaged_histograms) {{

            std::string nuisance_name = nuisance_name_numVariations.first;
            int numVariations_here = nuisance_name_numVariations.second;

            std::string string_pattern_up = "histo_.*__NORM__.*_" + string_pattern + "up";
            std::string string_pattern_do = "histo_.*__NORM__.*_" + string_pattern + "do";

            std::regex pattern_up(string_pattern_up);
            std::regex pattern_do(string_pattern_do);

            if (std::regex_search(temp_name, pattern_up) || std::regex_search(temp_name, pattern_do)) {{
              std::regex pattern_to_remove("__NORM__[0-9]+_");
              std::string histo_temp_name = std::regex_replace(temp_name, pattern_to_remove, "");
              std::regex pattern(string_pattern);

              histo_temp_name = std::regex_replace(histo_temp_name, pattern, nuisance_name);

              new_histo->SetName(histo_temp_name.c_str());

              if (std::regex_search(temp_name, pattern_up)) {{
                if (auto it = map_nuisances_average_histogram.find(string_pattern + "up"); it != map_nuisances_average_histogram.end()) {{
                  it->second.Add(new_histo, 1. / numVariations_here);
                }}
                else {{
                  map_nuisances_average_histogram[nuisance_name + "up"] = (*new_histo);
                  found_any_nuisances_for_averaged_histograms = true;
                  break;
                }}
              }}
              else {{
                if (auto it = map_nuisances_average_histogram.find(string_pattern + "do"); it != map_nuisances_average_histogram.end()) {{
                  it->second.Add(new_histo, 1. / numVariations_here);
                }}
                else {{
                  map_nuisances_average_histogram[nuisance_name + "do"] = (*new_histo);
                  found_any_nuisances_for_averaged_histograms = true;
                  break;
                }}
              }}
            }}
          }}
          if (not found_any_nuisances_for_averaged_histograms) {{
            new_histo->Write();
          }}

          for (auto const& [name, temp_histo] : map_nuisances_average_histogram) {{
            temp_histo.Write();
          }}

        }}
      }}
      out_file.cd();
    }}



    out_file.Close();

    return 0;
}}

        """


        with open(f"{output_name}.cpp", "w") as f:
            f.write(cpp_code)


    # _____________________________________________________________________________
    def compile_cpp(self, source_name):
        root_flags = subprocess.check_output(["root-config", "--cflags", "--libs"]).decode().strip()
        compile_cmd = f"g++ -O2 {source_name}.cpp -o {source_name} {root_flags}"
        print(f"Compiling...")
        result = os.system(compile_cmd)
        if result == 0:
            print(f"Compilation successful: ./{source_name}")
        else:
            print("Compilation failed!")


    # _____________________________________________________________________________
    def generate_makefile(self, file_paths, dictionary_of_files_divided_by_sample, makefile_name="Makefile"):
        cpp_flags = "$(shell root-config --cflags)"
        ld_flags = f"$(shell root-config --libs) -lTMVA -lXMLIO -I/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/include -L/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/lib -lcorrectionlib -Wl,-rpath,'$$ORIGIN:$$ORIGIN/../../'"

        type_of_compilation = "-O2"
        targets = [os.path.splitext(f)[0] for f in file_paths]

        with open(makefile_name, "w") as f:
            f.write("# Generated Makefile\n")
            f.write("CXX = g++\n")
            f.write(f"CXXFLAGS = {type_of_compilation} -Wall -fPIC {cpp_flags} -I{self._scripts_run_folder}/ \n")
            f.write(f"LDFLAGS = {ld_flags}\n\n")

            f.write(f"LIB_NAME = liblibrary_utils.so\n\n")

            all_targets = " ".join(targets)
            f.write(f"all: $(LIB_NAME) {all_targets}\n\n")

            for sample_key in dictionary_of_files_divided_by_sample.keys():
              target_sample = [os.path.splitext(f)[0] for f in dictionary_of_files_divided_by_sample[sample_key]]
              all_target_sample = " ".join(target_sample)
              f.write(f"{sample_key}: $(LIB_NAME) {all_target_sample}\n\n")

            f.write(f"$(LIB_NAME): {self._scripts_run_folder}/library_utils.cpp {self._scripts_run_folder}/library_utils.h\n")
            f.write(f"\t$(CXX) $(CXXFLAGS) -shared -o {self._scripts_run_folder}/$@ {self._scripts_run_folder}/library_utils.cpp\n")

            f.write("%: %.cpp $(LIB_NAME)\n")
            f.write(f"\t$(CXX) $(CXXFLAGS) $< -o $@ -L{self._scripts_run_folder}/ -llibrary_utils $(LDFLAGS)\n\n")

            f.write("clean:\n")
            f.write(f"\trm -f $(LIB_NAME) {all_targets}\n")


    # _____________________________________________________________________________
    def setValues(self, outputDir, variables, cuts, supercut, samples, nuisances, aliases, lumi,
                  scripts_run_folder=None, script_batch_location=None, extra_transfer_files=None):
        """Set analysis values. Optionally override folder locations for generated
        C++ sources (scripts_run_folder) and batch submission scripts
        (script_batch_location). If not provided, the instance defaults are kept."""

        self._variables = variables
        self._samples   = samples
        self._cuts      = cuts
        self._supercut  = supercut
        self._nuisances = nuisances
        self._aliases   = aliases
        self._outputDir = outputDir
        self._lumi      = lumi

        # Override folder locations only when explicitly provided.
        # os.path.expandvars handles ${MY_PATH_JOBS} and $MY_PATH_JOBS style
        # environment variables so users can set paths in configuration.py
        # using env-var placeholders (e.g. "${MY_PATH_JOBS}/2024/scripts_run").
        if scripts_run_folder is not None:
            self._scripts_run_folder = os.path.expandvars(scripts_run_folder)
        if script_batch_location is not None:
            self._script_batch_location = os.path.expandvars(script_batch_location)
        if extra_transfer_files is not None:
            self._extra_transfer_files = [os.path.expandvars(f) for f in extra_transfer_files]

        # Also expand env vars in outputDir for consistency
        self._outputDir = os.path.expandvars(outputDir)

    # _____________________________________________________________________________
    def setConditions(self, silentMode):

        self._silentMode = silentMode


    # _____________________________________________________________________________
    @staticmethod
    def _normalize_sample_name(name_field):
        """Accept both dict {'subname': [files]} and list of tuples [(subname, [files])]."""
        if isinstance(name_field, dict):
            return name_field
        # list of tuples: [(subname, [files]), ...]
        result = {}
        for item in name_field:
            subname, files = item
            if subname in result:
                result[subname].extend(files)
            else:
                result[subname] = list(files)
        return result

    # _____________________________________________________________________________
    def makeNominals(self):

        print ("======================")
        print ("==== makeNominals ====")
        print ("======================")

        os.makedirs(self._outputDir, exist_ok=True)

        ROOT.TH1.SetDefaultSumw2(True)

        os.makedirs(self._scripts_run_folder, exist_ok=True)

        self.create_library_cpp()

        list_of_files_to_compile = []
        list_of_files_to_compile_divided_by_sample = {}

        for sampleName, sample in self._samples.items():
          os.makedirs(self._scripts_run_folder + "/" + sampleName, exist_ok=True)

          list_of_files_to_compile_divided_by_sample[sampleName] = []

          global_weight = sample['weight'] if 'weight' in sample.keys() else "1."

          weights_subname = sample['weights'] if 'weights' in sample.keys() else {}

          if 'isData' in sample.keys() :
            if not (len(sample ['isData']) == 1 and sample ['isData'][0] == 'all') :
              pass
          else :
            global_weight = f"({global_weight}) * {self._lumi}"

          aliases_to_be_defined = {}
          functions_to_be_defined = {}
          # aliases_ordered preserves the original order from aliases.py
          # so that expr aliases depending on function aliases are defined after them
          aliases_ordered = {}
          for aliasName, alias in self._aliases.items():
            # if 'samples' is not specified, apply to all samples (default ALL)
            if 'samples' not in alias.keys() or sampleName in alias['samples']:
              if 'expr' in alias.keys() and 'linesToAdd' not in alias.keys():
                aliases_to_be_defined[aliasName] = alias['expr']
                aliases_ordered[aliasName] = ('expr', alias['expr'])
              if 'function' in alias.keys() or 'linesToAdd' in alias.keys() or 'linesToGlobal' in alias.keys():
                  functions_to_be_defined[aliasName] = alias
                  aliases_ordered[aliasName] = ('function', alias)


          for subname, list_root_files in self._normalize_sample_name(sample['name']).items():
            print ("length of list_root_files = ", len(list_root_files))
            os.makedirs(self._scripts_run_folder + "/" + sampleName + "/" + subname, exist_ok=True)

            if subname in weights_subname.keys():
              weight = f"({global_weight}) * ({weights_subname[subname]})"
            else :
              weight = global_weight

            make_job_every_N = 1
            if "FilesPerJob" in sample.keys():
              make_job_every_N = sample['FilesPerJob']

            chunks_list_root_files = [list_root_files[i:i + make_job_every_N] for i in range(0, len(list_root_files), make_job_every_N)]

            for i, root_files in enumerate(chunks_list_root_files):
              name_code = self._scripts_run_folder + "/" + sampleName + "/" + subname + "/" + "my_run_analysis_" + sampleName + "_" + subname + "_" + str(i)
              tree = "Events"

              output_root_file_name = "root_file___" + sampleName + "_" + subname + "_" + str(i) + ".root"

              subsamples_dict = {}
              if 'subsamples' in sample.keys():
                subsamples_dict = sample['subsamples']
              self.create_cpp_source_list_of_files(name_code, tree, root_files, output_root_file_name, sampleName, weight,
                                                   aliases_to_be_defined, functions_to_be_defined,
                                                   subsamples_dict, aliases_ordered)

              list_of_files_to_compile.append(name_code)
              list_of_files_to_compile_divided_by_sample[sampleName].append(name_code)


        self.generate_makefile(list_of_files_to_compile, list_of_files_to_compile_divided_by_sample)

        print("Do not compile ... If you want to compile run the compilation by hand!")

        subprocess.run(["make", "clean"])


    # _____________________________________________________________________________
    def checkBatch(self):

        print ("=====================")
        print ("==== check Batch ====")
        print ("=====================")
        print ("\n")

        submission_dir = os.getcwd()

        list_jobs_with_error = []
        list_missing_root_files = []
        list_running_jobs = []

        stats = {}

        for sampleName, sample in self._samples.items():
          stats[sampleName] = {'done': 0, 'running': 0, 'error': 0, 'missing': 0}
          for subname, list_root_files in self._normalize_sample_name(sample['name']).items():
            make_job_every_N = 1
            if "FilesPerJob" in sample.keys():
              make_job_every_N = sample['FilesPerJob']
            chunks_list_root_files = [list_root_files[i:i + make_job_every_N] for i in range(0, len(list_root_files), make_job_every_N)]
            for i, root_files in enumerate(chunks_list_root_files):
              name_err_file = self._script_batch_location + "/" + sampleName + "/" + subname + "/log/" + "my_script_" + str(i) + ".sh.err"
              name_log_file = self._script_batch_location + "/" + sampleName + "/" + subname + "/log/" + "my_script_" + str(i) + ".log"
              output_root_file_name = "root_file___" + sampleName + "_" + subname + "_" + str(i) + ".root"
              root_file_name = f"{submission_dir}/{self._outputDir}/{output_root_file_name}"

              if os.path.exists(root_file_name):
                stats[sampleName]['done'] += 1
              elif not os.path.exists(name_err_file):
                if os.path.exists(name_log_file):
                  with open(name_log_file, 'r') as f:
                    content = f.read()
                  if 'Job was aborted' in content:
                    list_jobs_with_error.append(name_log_file)
                    stats[sampleName]['error'] += 1
                  elif 'Job was evicted' in content or 'Finished transferring output files' in content:
                    list_missing_root_files.append(root_file_name)
                    stats[sampleName]['missing'] += 1
                  else:
                    list_running_jobs.append(root_file_name)
                    stats[sampleName]['running'] += 1
                else:
                  list_running_jobs.append(root_file_name)
                  stats[sampleName]['running'] += 1
              else:
                with open(name_err_file, 'r') as f:
                  lines = f.readlines()
                if len(lines) > 6:
                  list_jobs_with_error.append(name_err_file)
                  stats[sampleName]['error'] += 1
                else:
                  list_missing_root_files.append(root_file_name)
                  stats[sampleName]['missing'] += 1

        # print tabella
        col_w = max(len(s) for s in stats.keys()) + 2
        print(f"{'Sample':<{col_w}} {'Done':>8} {'Pending':>8} {'Error':>8}")
        print("-" * (col_w + 28))
        total = {'done': 0, 'running': 0, 'error': 0, 'missing': 0}
        for sampleName, s in stats.items():
          pending = s['running'] + s['missing']
          print(f"{sampleName:<{col_w}} {s['done']:>8} {pending:>8} {s['error']:>8}")
          for k in total:
            total[k] += s[k]
        print("-" * (col_w + 28))
        total_pending = total['running'] + total['missing']
        print(f"{'TOTAL':<{col_w}} {total['done']:>8} {total_pending:>8} {total['error']:>8}")


    def parallelCompile(self):

        print ("==========================")
        print ("==== compile on batch ====")
        print ("==========================")

        where_compile = "parallel_compile"
        os.makedirs(where_compile, exist_ok=True)

        submission_dir = os.getcwd()

        for sampleName, sample in self._samples.items():
          name_bash = where_compile + "/my_script_compile_" + sampleName + ".sh"

          bash_code = f"""#!/bin/bash
set -e  # Exit on error
echo "Job started at $(date)"
echo "Running on node $(hostname)"
cd {submission_dir}
make {sampleName}
"""

          with open(f"{name_bash}", "w") as f:
            f.write(bash_code)
            os.system ("chmod +x " + name_bash)


          name_submit = where_compile + "/my_script_compile_" + sampleName + ".sub"

          output_file = "/dev/null"
          error_file  = "/dev/null"
          log_file    = "/dev/null"

          submit_code = f"""
initialdir            = {submission_dir}
executable            = {name_bash}
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




    def submitBatch(self):

        print ("======================")
        print ("==== submit Batch ====")
        print ("======================")

        os.makedirs(self._script_batch_location, exist_ok=True)

        submission_dir = os.getcwd()
        abs_library  = os.path.join(self._scripts_run_folder, "liblibrary_utils.so")
        abs_lib_cpp  = os.path.join(self._scripts_run_folder, "library_utils.cpp")
        abs_lib_h    = os.path.join(self._scripts_run_folder, "library_utils.h")

        try:
            root_cflags = subprocess.check_output(["root-config", "--cflags"]).decode().strip()
            root_libs   = subprocess.check_output(["root-config", "--libs"]).decode().strip()
        except Exception:
            root_cflags = ""
            root_libs   = ""

        # Build liblibrary_utils.so locally if it does not exist yet.
        # This is fast (~5s) and needed in compile-on-node mode so the .so
        # can be transferred to the worker node.
        if not os.path.exists(abs_library):
            print("  -> liblibrary_utils.so not found, compiling it locally...")
            compile_lib_cmd = (
                f"g++ -O2 -Wall -fPIC {root_cflags} -shared "
                f"-o {abs_library} {abs_lib_cpp} {root_libs}"
            )
            ret = os.system(compile_lib_cmd)
            if ret == 0:
                print(f"  -> liblibrary_utils.so compiled: {abs_library}")
            else:
                print("  -> ERROR: failed to compile liblibrary_utils.so")

        # -------------------------------------------------------------------
        # Detect mode once by checking the first expected executable.
        # Either ALL are compiled or NONE are.
        # -------------------------------------------------------------------
        first_exe = None
        for sampleName, sample in self._samples.items():
          for subname, list_root_files in self._normalize_sample_name(sample['name']).items():
            if list_root_files:
              first_exe = os.path.join(self._scripts_run_folder, sampleName, subname,
                                       "my_run_analysis_" + sampleName + "_" + subname + "_0")
              break
          if first_exe:
            break

        precompiled = first_exe is not None and os.path.isfile(first_exe)

        if precompiled:
            print("  -> Pre-compiled executables found: transferring binaries to worker nodes.")
        else:
            print("  -> No pre-compiled executables found: worker nodes will compile before running.")

        # queue rows: (executable, job_exe, output, error, log)
        # "executable" is the .sh script (set per-row in the queue)
        # "job_exe"    is the binary or .cpp (set per-row, used in transfer_input_files via header)
        # The fixed extra transfer files (library/.so or lib sources) go in the header.
        # Detect the X509 proxy path at submission time so the worker node can use it
        _proxy_candidates = [
            os.environ.get('X509_USER_PROXY', ''),
            os.path.expanduser('~/.proxy'),
            f"/tmp/x509up_u{os.getuid()}",
        ]
        x509_proxy = next((p for p in _proxy_candidates if p and os.path.exists(p)), None)
        if x509_proxy:
            # Also check the proxy is actually valid (not expired)
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

        for sampleName, sample in self._samples.items():
          os.makedirs(os.path.join(self._script_batch_location, sampleName), exist_ok=True)

          for subname, list_root_files in self._normalize_sample_name(sample['name']).items():
            os.makedirs(os.path.join(self._script_batch_location, sampleName, subname), exist_ok=True)

            make_job_every_N = sample.get('FilesPerJob', 1)
            chunks_list_root_files = [list_root_files[i:i + make_job_every_N]
                                      for i in range(0, len(list_root_files), make_job_every_N)]

            for i, root_files in enumerate(chunks_list_root_files):
              name_code_no_folder   = "my_run_analysis_" + sampleName + "_" + subname + "_" + str(i)
              name_bash             = os.path.join(self._script_batch_location, sampleName, subname,
                                                   "my_script_" + str(i) + ".sh")
              output_root_file_name = "root_file___" + sampleName + "_" + subname + "_" + str(i) + ".root"
              name_folder      = os.path.join(self._script_batch_location, sampleName, subname)
              name_folder_code = os.path.join(self._scripts_run_folder, sampleName, subname)
              abs_executable   = os.path.join(name_folder_code, name_code_no_folder)
              abs_cpp          = os.path.join(name_folder_code, name_code_no_folder + ".cpp")

              # Build the copy command.
              # If destination is on EOS, avoid hardcoding the /eos path in the script
              # (HTCondor scans executable content and rejects /eos paths).
              # Instead, build the destination at runtime inside the script using shell variables.
              dest_dir = f"{submission_dir}/{self._outputDir}"
              if submission_dir.startswith('/eos'):
                # Normalize /eos/home-X/username to /eos/user/X/username
                # which is the correct path for root://eosuser.cern.ch
                import re as _re
                dest_dir_normalized = _re.sub(
                    r'^/eos/home-([a-z])/',
                    lambda m: f"/eos/user/{m.group(1)}/",
                    dest_dir
                )
                cp_cmd = f"xrdcp {self._outputDir}/{output_root_file_name} root://eosuser.cern.ch//{dest_dir_normalized}/{output_root_file_name}"
              else:
                cp_cmd = f"cp {self._outputDir}/{output_root_file_name}  {dest_dir}/"

              if precompiled:
                bash_code = f"""#!/bin/bash
set -e
echo "Job started at $(date)"
echo "Running on node $(hostname)"
source /cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/setup.sh
{proxy_lines}mkdir -p {self._outputDir}
mkdir -p {sampleName}/{subname}/
mv {name_code_no_folder} {sampleName}/{subname}/
./{sampleName}/{subname}/{name_code_no_folder}
{cp_cmd}
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
{proxy_lines}echo "Compiling analysis executable..."
g++ -O2 -Wall -fPIC {root_cflags} -I. \\
    {name_code_no_folder}.cpp -o {name_code_no_folder} \\
    -L. -llibrary_utils \\
    -lTMVA -lXMLIO \\
    -I/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/include \\
    -L/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/lib \\
    -lcorrectionlib {root_libs} \\
    -Wl,-rpath,'$ORIGIN'
echo "Compilation done at $(date)"
mkdir -p {self._outputDir}
mkdir -p {sampleName}/{subname}/
mv liblibrary_utils.so {sampleName}/{subname}/
mv {name_code_no_folder} {sampleName}/{subname}/
export LD_LIBRARY_PATH=$(pwd)/{sampleName}/{subname}/:$LD_LIBRARY_PATH
./{sampleName}/{subname}/{name_code_no_folder}
{cp_cmd}
echo "Done at $(date)"
ls -lh
"""
                job_exe = abs_cpp

              with open(name_bash, "w") as f:
                f.write(bash_code)
                os.system("chmod +x " + name_bash)

              if self._silentMode:
                output_file = "/dev/null"
                error_file  = "/dev/null"
                log_file    = "/dev/null"
              else:
                log_dir = os.path.join(name_folder, "log")
                os.makedirs(log_dir, exist_ok=True)
                output_file = os.path.join(log_dir, f"my_script_{i}.out")
                error_file  = os.path.join(log_dir, f"my_script_{i}.err")
                log_file    = os.path.join(log_dir, f"my_script_{i}.log")

              queue_rows.append((name_bash, job_exe, output_file, error_file, log_file))

        # ------------------------------------------------------------------
        # Single .sub file using "queue executable, job_exe, output, error, log from (...)"
        # From the HTCondor docs (queue var from (...) syntax), submit commands
        # like "executable", "output", "error", "log" can be set as per-row variables.
        # The fixed extra transfer files (same for all jobs) go in the header.
        # This avoids any macro-expansion issues.
        # ------------------------------------------------------------------
        global_submit_path = os.path.join(self._script_batch_location, "submit_all.sub")

        # Extra files to transfer (e.g. btag efficiency maps, macros)
        # Set via extra_transfer_files in configuration.py
        extra_files = getattr(self, '_extra_transfer_files', [])
        extra_files_str = (", " + ", ".join(extra_files)) if extra_files else ""

        if precompiled:
            fixed_extra = f"{abs_library}, {x509_proxy}{extra_files_str}"
        else:
            fixed_extra = f"{abs_lib_cpp}, {abs_lib_h}, {abs_library}, {x509_proxy}{extra_files_str}"

        queue_lines = "\n".join(
            f"  {exe}, {job_exe}, {out}, {err}, {log}"
            for exe, job_exe, out, err, log in queue_rows
        )

        submit_code = (
            "# Auto-generated by mkShapes.py\n"
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
        # Run condor_submit from the directory containing the .sub file,
        # so that HTCondor sees an AFS working directory regardless of where
        # mkShapes.py was launched from.
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
         '                ___|   |                               \\  |         |                           '
         '              \\___ \\   __ \\    _` |  __ \\    _ \\      |\\/ |   _` |  |  /   _ \\   __|            '
         '                    |  | | |  (   |  |   |   __/      |   |  (   |    <    __/  |               '
         '              _____/  _| |_| \\__,_|  .__/  \\___|     _|  _| \\__,_| _|\\_\\ \\___| _|               '
         '                                    _|                                                          '
         '                                                                                                '
         --------------------------------------------------------------------------------------------------
         """

    print(header)

    parser = defaultParser()

    parser.add_argument("--parallelCompile", action='store_true', dest="parallelCompile",  help="Trigger the parallel compilation of the code, one per sample")
    parser.add_argument("--localCompile",    action='store_true', dest="localCompile",     help="Trigger the local compilation of the code")
    parser.add_argument("--submitBatch",     action='store_true', dest="submitBatch",      help="Trigger the submission to lxbatch")
    parser.add_argument("--hadd",            action='store_true', dest="haddRootFiles",    help="Trigger the merging of the root files")
    parser.add_argument("--checkBatch",      action='store_true', dest="checkBatch",       help="Check if jobs are done succesfully")
    parser.add_argument("--silentMode",      action='store_true', dest="silentMode",       help="Remove as much as possible the print and the log/err files production")

    opt = parser.parse_args()
    print ("opt.pycfg            = ", opt.pycfg)
    loadDefaultOptions(parser, opt.pycfg)
    opt = parser.parse_args()


    print ("opt.tag              = ", opt.tag)
    print ("opt.variablesFile    = ", opt.variablesFile)
    print ("opt.cutsFile         = ", opt.cutsFile)
    print ("opt.samplesFile      = ", opt.samplesFile)
    print ("opt.nuisancesFile    = ", opt.nuisancesFile)
    print ("opt.lumi             = ", opt.lumi)
    print ("opt.outputDir        = ", opt.outputDir)


    print ("opt.parallelCompile  = ", opt.parallelCompile)
    print ("opt.localCompile     = ", opt.localCompile)
    print ("opt.submitBatch      = ", opt.submitBatch)
    print ("opt.hadd             = ", opt.haddRootFiles)
    print ("opt.checkBatch       = ", opt.checkBatch)
    print ("opt.silentMode       = ", opt.silentMode)

    #
    # read list of samples
    #
    samples = {}
    if os.path.exists(opt.samplesFile) :
      handle = open(opt.samplesFile,'r')
      exec(handle.read())
      handle.close()
      samples = {k: v for k, v in samples.items() if not (k.startswith('__') and k.endswith('__'))}

    #
    # read list of variables
    #
    variables = {}
    if os.path.exists(opt.variablesFile) :
      handle = open(opt.variablesFile,'r')
      exec(handle.read())
      handle.close()
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
      cuts = {k: v for k, v in cuts.items() if not (k.startswith('__') and k.endswith('__'))}

    print ("cuts = ", cuts)
    print ("supercut = ", supercut)

    #
    # read list of nuisances
    #
    nuisances = {}
    if os.path.exists(opt.nuisancesFile) :
      handle = open(opt.nuisancesFile,'r')
      exec(handle.read())
      handle.close()
      nuisances = {k: v for k, v in nuisances.items() if not (k.startswith('__') and k.endswith('__'))}

    print ("nuisances = ", nuisances)

    #
    # read list of aliases
    #
    aliases = {}
    if os.path.exists(opt.aliasesFile) :
      handle = open(opt.aliasesFile,'r')
      exec(handle.read())
      handle.close()
      aliases = {k: v for k, v in aliases.items() if not (k.startswith('__') and k.endswith('__'))}

    print ("aliases = ", aliases)


    factory = ShapeFactory()

    # -----------------------------------------------------------------------
    # Read folder locations from configuration.py (via opt), falling back to
    # the ShapeFactory defaults when not set there.
    # -----------------------------------------------------------------------
    scripts_run_folder    = getattr(opt, 'scripts_run_folder',    None)
    script_batch_location = getattr(opt, 'script_batch_location', None)

    factory.setValues(
        opt.outputDir, variables, cuts, supercut, samples, nuisances, aliases, opt.lumi,
        scripts_run_folder    = scripts_run_folder,
        script_batch_location = script_batch_location,
        extra_transfer_files  = getattr(opt, 'extra_transfer_files', None),
    )
    factory.setConditions(opt.silentMode)

    print ("factory._scripts_run_folder    = ", factory._scripts_run_folder)
    print ("factory._script_batch_location = ", factory._script_batch_location)

    if not opt.submitBatch and not opt.haddRootFiles and not opt.checkBatch and not opt.localCompile and not opt.parallelCompile:
      factory.makeNominals()

    if opt.localCompile :
      subprocess.run(["make", "-j4"])

    if opt.parallelCompile :
      factory.parallelCompile()

    if opt.submitBatch :
      factory.submitBatch()

    if opt.haddRootFiles :
      print ("Hadd the different root files ...")

      hadd_cmd = f"hadd -j 8 -f {opt.outputDir}/root_file_joined.root    {opt.outputDir}/root_file___*.root"

      print(f"hadd: {hadd_cmd}")
      result = os.system(hadd_cmd)

      print ("I have hadded all the root files but I have not removed the original ones")
      print ("I have performed an hadd of all the suitable root files in the folder ... ")


    if opt.checkBatch :
      print ("Checking if the jobs finished succesfully")
      factory.checkBatch()



    print ("\n\n")
    print (" I'm done ... \n\n")