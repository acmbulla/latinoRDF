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
                                        subsamples):

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
          else:  # if "xmlfile" is present, it's a TMVA and it requires casting to float
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
              register_variations_logic_for_alternative_trees += f'''          std::string expression = "ROOT::RVec<" + varied_df.GetColumnType(nomCol) + ">{{static_cast<" + varied_df.GetColumnType(nomCol) + ">(" + branch + ")}}";\n'''
              #
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
              register_variations_logic_for_alternative_trees += f'''          std::string expression = "ROOT::RVec<" + varied_df.GetColumnType(nomCol) + ">{{static_cast<" + varied_df.GetColumnType(nomCol) + ">(" + branch + ")}}";\n'''
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




        # # Write C++ code to create a specific node for this sub_name
        # cpp_file.write(f'auto df_{sub_name} = df_base.Filter("{sub_cut}");\n')
        #
        # # If the subsample has a specific weight:
        # sub_weight = subsample_weights.get(sub_name, "1.0")
        # cpp_file.write(f'auto df_{sub_name}_final = df_{sub_name}.Define("total_weight", "{sub_weight}");\n')



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




                # # We add each RResultPtr to a vector called 'histograms'
                # define_variables_logic += f'    hist_map_1D["{this_cutName}"].push_back(node_{this_cutName}.{histo_call});\n'
                # define_variables_logic += f'    hist_map["{this_cutName}"].push_back(node_{this_cutName}.{histo_call});\n'
                # define_variables_logic += f'    hist_map["{this_cutName}"].push_back(ROOT::RDF::RResultPtr<TH1D> (node_{this_cutName}.{histo_call}) );\n'

                # We add each RResultPtr to a vector called 'histograms'
                # define_variables_logic += f'    hist_map["{this_cutName}"].push_back(node_{this_cutName}.Histo1D({{"h_{variableName}", "{variableName}", {model_str}}}, "{variableName}", "my_sample_weight"));\n'


###                 (bins, v_min, v_max) = variable['range']
###                 # We add each RResultPtr to a vector called 'histograms'
###                 define_variables_logic += f'    hist_map["{this_cutName}"].push_back(node_{this_cutName}.Histo1D({{"h_{variableName}", "{variableName}", {bins}, {v_min}, {v_max}}}, "{variableName}", "my_sample_weight"));\n'


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


// --- Automatically generated: code to be added for additional functions ---

{code_of_function_to_include}

// ----------------------------------------



// --- Automatically generated: code to be added for special nuisances ---

{code_of_functions_and_classes_to_include_for_special_nuisances}

// ----------------------------------------


int main() {{
    ROOT::EnableImplicitMT();

    // --- Automatically generated input root files ---

    {define_input_files_logic}

    // ----------------------------------------

    std::map<std::string, std::vector<ROOT::RDF::RResultPtr<TH1D>>> hist_map_1D;
    std::map<std::string, std::vector<ROOT::RDF::RResultPtr<TH2D>>> hist_map_2D;

    // --- Automatically generated register nuisances variations for alternative trees ---
    //        Alternative trees variation nuisances to be defined at the beginning because they might change the
    //        variables that are defined downstream, e.g. njet, BDT(ptjet, ptlepton, ...)
    //

    {register_variations_logic_for_alternative_trees}



    // --- Automatically generated register nuisances variations for alternative-like trees, as they are calculated on the fly ---

    {register_variations_for_special_nuisances_on_the_fly}

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


    // --- Automatically generated register nuisances variations ---

    {register_variations_logic}

    // ----------------------------------------

    // --- Automatically generated ... just the main node ---

    {booking_logic}

    //
    // from now on: current_node
    // ----------------------------------------

    // --- Automatically generated define cuts and histograms ---

    {define_variables_logic}

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

    std::cout << "   ROOT::RDF::RunGraphs 2D" << std::endl;
    if (!hist_map_2D.empty()) {{
      ROOT::RDF::RunGraphs(all_booking_2D);
    }}


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
            //size_t pos = temp_name.find(':');
            //if (pos != std::string::npos) {{
            //  temp_name = temp_name.substr(0, pos);
            //}}
            if (sub_name == "") temp_name = ("histo_{sampleName}_" + temp_name);
            else                temp_name = ("histo_{sampleName}_" + sub_name + "_" + temp_name);
          }}
          gDirectory = subdir;
          histo->SetName(temp_name.c_str());

          if (list_of_variables_fold_1D.at(ivar) != 0) FoldHistogram(histo.get(), list_of_variables_fold_1D.at(ivar));

          //
          // special treatment for nuisances based on averaged histograms
          //
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

              histo_temp_name = std::regex_replace(histo_temp_name, pattern, nuisance_name); // replace with real nuisance name

              histo->SetName(histo_temp_name.c_str());

              if (std::regex_search(temp_name, pattern_up)) {{
                if (auto it = map_nuisances_average_histogram.find(string_pattern + "up"); it != map_nuisances_average_histogram.end()) {{
                  it->second.Add(histo.get(), 1. / numVariations_here);
                }}
                else {{
                  map_nuisances_average_histogram[nuisance_name + "up"] = (*histo);
                  found_any_nuisances_for_averaged_histograms = true;
                  break; // end loop over map to save time
                }}
              }}
              else {{
                if (auto it = map_nuisances_average_histogram.find(string_pattern + "do"); it != map_nuisances_average_histogram.end()) {{
                  it->second.Add(histo.get(), 1. / numVariations_here);
                }}
                else {{
                  map_nuisances_average_histogram[nuisance_name + "do"] = (*histo);
                  found_any_nuisances_for_averaged_histograms = true;
                  break; // end loop over map to save time
                }}
              }}
            }}
          }}
          // if it's NOT a histogram to be averaged, then write it directly
          if (not found_any_nuisances_for_averaged_histograms) {{
            histo->Write();
          }}

          // now write all the averaged histograms, up and down variations, if there are any!
          for (auto const& [name, temp_histo] : map_nuisances_average_histogram) {{
            temp_histo.Write();
          }}


        }}
      }}
      out_file.cd(); // Go back to root for the next directory
    }}

    big_loop = 0;
    for (auto& [cut_label, h_list] : hist_map_2D) {{

      std::string current_cut = cut_label;
      // check if subsamples case: this_cutName___sampleName_sub_name
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

      // Create a folder for this cut
      TDirectory *dir = out_file.GetDirectory(current_cut.c_str());
      if (!dir) {{
        // If it doesn't exist, create it
        dir = out_file.mkdir(current_cut.c_str());
      }}
      for (size_t ivar = 0; ivar < h_list.size(); ivar++) {{
        dir->cd();
        TDirectory *subdir = out_file.GetDirectory( (current_cut+"/"+list_of_variables_2D.at(ivar)).c_str() );
        if (!subdir) {{
          // If it doesn't exist, create it
          subdir = out_file.mkdir( (current_cut+"/"+list_of_variables_2D.at(ivar)).c_str() );
        }}
        subdir->cd();

        // get the nominal and the variations too
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
            // scale_e_2018_UL:up --> scale_e_2018_ULup
            temp_name.erase(std::remove(temp_name.begin(), temp_name.end(), ':'), temp_name.end());
            //size_t pos = temp_name.find(':');
            //if (pos != std::string::npos) {{
            //  temp_name = temp_name.substr(0, pos);
            //}}
            if (sub_name == "") temp_name = ("histo_{sampleName}_" + temp_name);
            else                temp_name = ("histo_{sampleName}_" + sub_name + "_" + temp_name);
          }}
          gDirectory = subdir;
          histo->SetName(temp_name.c_str());


          if (list_of_variables_fold_2D.at(ivar) != 0) FoldHistogram(histo.get(), list_of_variables_fold_2D.at(ivar));

          //UnrollHistogram(dynamic_cast<TH2D*>(histo.get()))->Write();

          auto new_histo = UnrollHistogram(dynamic_cast<TH2D*>(histo.get()));


          //
          // special treatment for nuisances based on averaged histograms
          //
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

              histo_temp_name = std::regex_replace(histo_temp_name, pattern, nuisance_name); // replace with real nuisance name

              new_histo->SetName(histo_temp_name.c_str());

              if (std::regex_search(temp_name, pattern_up)) {{
                if (auto it = map_nuisances_average_histogram.find(string_pattern + "up"); it != map_nuisances_average_histogram.end()) {{
                  it->second.Add(new_histo, 1. / numVariations_here);
                }}
                else {{
                  map_nuisances_average_histogram[nuisance_name + "up"] = (*new_histo);
                  found_any_nuisances_for_averaged_histograms = true;
                  break; // end loop over map to save time
                }}
              }}
              else {{
                if (auto it = map_nuisances_average_histogram.find(string_pattern + "do"); it != map_nuisances_average_histogram.end()) {{
                  it->second.Add(new_histo, 1. / numVariations_here);
                }}
                else {{
                  map_nuisances_average_histogram[nuisance_name + "do"] = (*new_histo);
                  found_any_nuisances_for_averaged_histograms = true;
                  break; // end loop over map to save time
                }}
              }}
            }}
          }}
          // if it's NOT a histogram to be averaged, then write it directly
          if (not found_any_nuisances_for_averaged_histograms) {{
            new_histo->Write();
          }}

          // now write all the averaged histograms, up and down variations, if there are any!
          for (auto const& [name, temp_histo] : map_nuisances_average_histogram) {{
            temp_histo.Write();
          }}


        }}
      }}
      out_file.cd(); // Go back to root for the next directory
    }}



    out_file.Close();

    return 0;
}}

        """


        with open(f"{output_name}.cpp", "w") as f:
            f.write(cpp_code)
        # print(f"Created {output_name}.cpp")




    # _____________________________________________________________________________
    def compile_cpp(self, source_name):
        # Get ROOT flags automatically using 'root-config'
        root_flags = subprocess.check_output(["root-config", "--cflags", "--libs"]).decode().strip()

        # Construct the compilation command
        compile_cmd = f"g++ -O2 {source_name}.cpp -o {source_name} {root_flags}"
        # compile_cmd = f"g++ -O3 {source_name}.cpp -o {source_name} {root_flags}"

        print(f"Compiling...")
        result = os.system(compile_cmd)

        if result == 0:
            print(f"Compilation successful: ./{source_name}")
        else:
            print("Compilation failed!")


    # _____________________________________________________________________________
    def generate_makefile(self, file_paths, dictionary_of_files_divided_by_sample, makefile_name="Makefile"):
        # Get ROOT configuration via shell calls
        cpp_flags = "$(shell root-config --cflags)"
        # ld_flags = "$(shell root-config --libs) -lTMVA -lXMLIO"  # why TMVA is not by default?!?
        # ld_flags = f"$(shell root-config --libs) -lTMVA -lXMLIO -Wl,-rpath,'$$ORIGIN:$$ORIGIN{self._scripts_run_folder}/'"

        # ld_flags = f"$(shell root-config --libs) -lTMVA -lXMLIO -I/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/include -L/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/lib -lcorrectionlib -Wl,-rpath,'$$ORIGIN:$$ORIGIN{self._scripts_run_folder}/'"

        # $ORIGIN is where the executable is
        # the libraries are in ../../
        # NB: FIXME it might be useful to move them into a "lib" folder later ... for cleaner setup, but not mandatory
        #
        ld_flags = f"$(shell root-config --libs) -lTMVA -lXMLIO -I/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/include -L/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/lib -lcorrectionlib -Wl,-rpath,'$$ORIGIN:$$ORIGIN/../../'"

        # ld_flags = f"$(shell root-config --libs) -lTMVA -lXMLIO -I/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/include -L/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/lib -lcorrectionlib -Wl,-rpath,'$$ORIGIN:$$ORIGIN{self._scripts_run_folder}/'"

        # -I/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/include -L/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/lib -lcorrectionlib


        # ld_flags = f"$(shell root-config --libs) -lTMVA -lXMLIO -I/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/include -L/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/lib -lcorrectionlib -Wl,-rpath,'$$ORIGIN:$$ORIGIN{self._scripts_run_folder}/'"

        # ld_flags = f"$(shell root-config --libs) -lTMVA -lXMLIO -I/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/include -L/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/lib -lcorrectionlib -Wl,-rpath,'$$ORIGIN:$$ORIGIN{self._scripts_run_folder}/'"



        type_of_compilation = "-O2"
        # from "blabla.cpp" to "blabla"
        targets = [os.path.splitext(f)[0] for f in file_paths]


        with open(makefile_name, "w") as f:
            f.write("# Generated Makefile\n")
            f.write("CXX = g++\n")
            # f.write(f"CXXFLAGS = {type_of_compilation} -Wall {cpp_flags} -std=c++20 -Wcpp \n")  # c++20 used for "ends_with"
            # f.write(f"CXXFLAGS = {type_of_compilation} -Wall {cpp_flags} \n")
            f.write(f"CXXFLAGS = {type_of_compilation} -Wall -fPIC {cpp_flags} -I{self._scripts_run_folder}/ \n")
            f.write(f"LDFLAGS = {ld_flags}\n\n")

            f.write(f"LIB_NAME = liblibrary_utils.so\n\n")


            # 'all' target: the list of all executables to be created
            all_targets = " ".join(targets)
            f.write(f"all: $(LIB_NAME) {all_targets}\n\n")


            # from "blabla.cpp" to "blabla"
            for sample_key in dictionary_of_files_divided_by_sample.keys():
              target_sample = [os.path.splitext(f)[0] for f in dictionary_of_files_divided_by_sample[sample_key]]
              all_target_sample = " ".join(target_sample)
              f.write(f"{sample_key}: $(LIB_NAME) {all_target_sample}\n\n")


            f.write(f"$(LIB_NAME): {self._scripts_run_folder}/library_utils.cpp {self._scripts_run_folder}/library_utils.h\n")
            f.write(f"\t$(CXX) $(CXXFLAGS) -shared -o {self._scripts_run_folder}/$@ {self._scripts_run_folder}/library_utils.cpp\n")

            # The Pattern Rule
            # This says: To create ANY executable 'X', look for 'X.cpp'
            # It works across subdirectories automatically.
            f.write("%: %.cpp $(LIB_NAME)\n")
            # f.write("\t$(CXX) $(CXXFLAGS) $< -o $@ -L. -llibrary_utils $(LDFLAGS)\n\n")
            f.write(f"\t$(CXX) $(CXXFLAGS) $< -o $@ -L{self._scripts_run_folder}/ -llibrary_utils $(LDFLAGS)\n\n")

            # Clean target
            f.write("clean:\n")
            f.write(f"\trm -f $(LIB_NAME) {all_targets}\n")


    # _____________________________________________________________________________
    def setValues(self, outputDir, variables, cuts, supercut, samples, nuisances, aliases, lumi):

        self._variables = variables
        self._samples   = samples
        self._cuts      = cuts
        self._supercut  = supercut
        self._nuisances = nuisances
        self._aliases   = aliases
        self._outputDir = outputDir
        self._lumi      = lumi

    # _____________________________________________________________________________
    def setConditions(self, silentMode):

        self._silentMode = silentMode


    # _____________________________________________________________________________
    def makeNominals(self):

        print ("======================")
        print ("==== makeNominals ====")
        print ("======================")

        os.system ("mkdir " + self._outputDir + "/")

        ROOT.TH1.SetDefaultSumw2(True)

        os.system ("mkdir " + self._scripts_run_folder + "/")

        #
        # create the header and cpp code for the library with useful functions
        #
        self.create_library_cpp()

        list_of_files_to_compile = []
        list_of_files_to_compile_divided_by_sample = {}

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
          # print ("length of sample['name'] = ", len(sample['name']))

          list_of_files_to_compile_divided_by_sample[sampleName] = []

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
              # name_code = self._scripts_run_folder + "/" + sampleName + "/" + subname + "/" + "my_run_analysis_" + str(i)
              name_code = self._scripts_run_folder + "/" + sampleName + "/" + subname + "/" + "my_run_analysis_" + sampleName + "_" + subname + "_" + str(i)
              tree = "Events"

              output_root_file_name = "root_file___" + sampleName + "_" + subname + "_" + str(i) + ".root"

              subsamples_dict = {}
              if 'subsamples' in sample.keys():
                subsamples_dict = sample['subsamples']
              self.create_cpp_source_list_of_files(name_code, tree, root_files, output_root_file_name, sampleName, weight,
                                                   aliases_to_be_defined, functions_to_be_defined,
                                                   subsamples_dict)

              list_of_files_to_compile.append(name_code)
              list_of_files_to_compile_divided_by_sample[sampleName].append(name_code)

              # self.compile_cpp(name_code)


            # for i, root_file in enumerate(list_root_files):
            #   # name_code = self._scripts_run_folder + "/" + sampleName + "/" + subname + "/" + "my_run_analysis_" + str(i)
            #   name_code = self._scripts_run_folder + "/" + sampleName + "/" + subname + "/" + "my_run_analysis_" + sampleName + "_" + subname + "_" + str(i)
            #   tree = "Events"
            #
            #   output_root_file_name = "root_file___" + sampleName + "_" + subname + "_" + str(i) + ".root"
            #
            #   self.create_cpp_source_single_file(name_code, tree, root_file, output_root_file_name, sampleName, weight, aliases_to_be_defined)
            #
            #   list_of_files_to_compile.append(name_code)
            #   # self.compile_cpp(name_code)


        self.generate_makefile(list_of_files_to_compile, list_of_files_to_compile_divided_by_sample)

        # Run the compilation in parallel
        # print("Start parallel compilation...")
        print("Do not compile ... If you want to compile run the compilation by hand!")

        subprocess.run(["make", "clean"])
        # subprocess.run(["make", "-j4"])
        # subprocess.run(["make", "-j8"])
        #
        # subprocess.run() is a blocking call -> it will make the compilation to end before the next step
        #

        # Now you can run it:
        # subprocess.run([f"./{name_code}"])


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
              if os.path.exists(name_err_file):
                with open(name_err_file, 'r') as f:
                  line_count = sum(1 for line in f)
                  if line_count > 6: # should I remove the warnings to have 0? FIXME
                    list_jobs_with_error.append(name_err_file)
              output_root_file_name = "root_file___" + sampleName + "_" + subname + "_" + str(i) + ".root"
              root_file_name = f"{submission_dir}/{self._outputDir}/{output_root_file_name}"
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



    def parallelCompile(self):

        print ("==========================")
        print ("==== compile on batch ====")
        print ("==========================")

        where_compile = "parallel_compile"
        os.system ("mkdir " + where_compile + "/")

        submission_dir = os.getcwd()

        #
        # Loop over samples
        #
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
              name_code = self._scripts_run_folder + "/" + sampleName + "/" + subname + "/" + "my_run_analysis_" + sampleName + "_" + subname + "_" + str(i)
              name_code_no_folder = "my_run_analysis_" + sampleName + "_" + subname + "_" + str(i)
              name_bash = self._script_batch_location + "/" + sampleName + "/" + subname + "/" + "my_script_" + str(i) + ".sh"
              output_root_file_name = "root_file___" + sampleName + "_" + subname + "_" + str(i) + ".root"

#
# moving into a folder structure is needed to handle libraries properly
#

              bash_code = f"""#!/bin/bash
set -e  # Exit on error
echo "Job started at $(date)"
echo "Running on node $(hostname)"
mkdir {self._outputDir}
mkdir {sampleName}/
mkdir {sampleName}/{subname}/
mv {name_code_no_folder} {sampleName}/{subname}/
./{sampleName}/{subname}/{name_code_no_folder}
cp {self._outputDir}/{output_root_file_name}  {submission_dir}/{self._outputDir}/

echo "Current directory content after running:"
ls -lh
echo "Current full path: $(pwd)"

"""

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
transfer_input_files  = {submission_dir}/{name_folder_code}{name_code_no_folder}, {submission_dir}/{self._scripts_run_folder}/liblibrary_utils.so
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




    # not used by mkShapes
    # print ("opt.structureFile    = ", opt.structureFile)

    #
    # logic of dependencies:
    #
    #    samples.py
    #      --> nuisances.py   since some nuisances are defined only for some samples
    #      --> aliases.py     since the aliases/defines will be "Defined" only for selected samples
    #    cuts.py
    #      --> nuisances.py   since some nuisances are defined only for some cuts (e.g. lnN of migration)
    #    nuisances.py
    #    variables.py
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
    # read list of nuisances
    #
    nuisances = {}
    if os.path.exists(opt.nuisancesFile) :
      handle = open(opt.nuisancesFile,'r')
      exec(handle.read())
      handle.close()
      # clean the dictionary to remove globals due to "exec" funcionality
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
      # clean the dictionary to remove globals due to "exec" funcionality
      aliases = {k: v for k, v in aliases.items() if not (k.startswith('__') and k.endswith('__'))}

    print ("aliases = ", aliases)






    factory = ShapeFactory()
    factory.setValues( opt.outputDir, variables, cuts, supercut, samples, nuisances , aliases, opt.lumi)
    factory.setConditions (opt.silentMode)

    # factory._treeName  = opt.treeName
    # factory._energy    = opt.energy
    # factory._lumi      = opt.lumi
    # factory._tag       = opt.tag

    if not opt.submitBatch and not opt.haddRootFiles and not opt.checkBatch and not opt.localCompile and not opt.parallelCompile:
      factory.makeNominals()

    if opt.localCompile :
      # no need to call the factory, just compile
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