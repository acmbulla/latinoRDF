#!/usr/bin/env python3

import sys
argv = sys.argv
sys.argv = argv[:1]

import os
import subprocess
import logging
import ROOT


sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysisFactory import (AnalysisFactory, defaultParser, loadDefaultOptions,
                             read_exec, load_config_files, _SCRIPTS_DIR)

# Variation tags used by RDataFrame Vary(): they end up as suffixes of the histogram names
# and must match what combine expects ($SYSTEMATIC -> <nuisance>Up / <nuisance>Down).
TAG_UP, TAG_DOWN = 'Up', 'Down'


# -----------------------------------------------------------------------
#  ShapeFactory  —  mkShapes-specific implementation
# -----------------------------------------------------------------------

class ShapeFactory(AnalysisFactory):
    """Produces ROOT files with TH1D/TH2D histograms per cut/variable/sample."""

    _logger = logging.getLogger('ShapeFactory')

    def __init__(self):
        super().__init__()
        self._outputRootFile = 'myhistos.root'

    # _____________________________________________________________________________
    def job_prefix(self):
        return "my_run_analysis"

    def output_file_name(self, sampleName, subname, i):
        return f"root_file___{sampleName}_{subname}_{i}.root"

    def _output_file_from_suffix(self, sampleName, subname, suffix):
        return f"root_file___{suffix}.root"

    # _____________________________________________________________________________
    #  library_utils.cpp/h  (mkShapes-specific — not needed by mkSnapshot)
    # _____________________________________________________________________________

    def create_library_cpp(self):
        cpp_code_library_header = """
#ifndef LIBRARY_UTILS_H
#define LIBRARY_UTILS_H

#include <string>
#include <vector>
#include "TH1.h"
#include "TTree.h"
#include "ROOT/RDataFrame.hxx"

void FoldHistogram(TH1* h, int doFold);
TH1D* UnrollHistogram(TH2D* h2);
ROOT::RDF::RNode SafeDefine(ROOT::RDF::RNode df, std::string var, std::string expr);
std::vector<std::string> getBranchNames(TTree* tree);

#endif // LIBRARY_UTILS_H
"""
        cpp_code_library = r"""
#include "library_utils.h"

void FoldHistogram(TH1* h, int doFold) {
  if (h->GetDimension() == 1) {
    int nBins = h->GetNbinsX();
    if (doFold == 1 || doFold == 3) {
      double content_last     = h->GetBinContent(nBins);
      double error_last       = h->GetBinError(nBins);
      double content_overflow = h->GetBinContent(nBins+1);
      double error_overflow   = h->GetBinError(nBins+1);
      h->SetBinContent(nBins, content_last + content_overflow);
      h->SetBinError(nBins, std::sqrt(error_last*error_last + error_overflow*error_overflow));
      h->SetBinContent(nBins + 1, 0); h->SetBinError(nBins + 1, 0);
    }
    if (doFold == 2 || doFold == 3) {
      double content_first     = h->GetBinContent(1);
      double error_first       = h->GetBinError(1);
      double content_underflow = h->GetBinContent(0);
      double error_underflow   = h->GetBinError(0);
      h->SetBinContent(1, content_first + content_underflow);
      h->SetBinError(1, std::sqrt(error_first*error_first + error_underflow*error_underflow));
      h->SetBinContent(0, 0); h->SetBinError(0, 0);
    }
  } else {
    int nX = h->GetNbinsX(), nY = h->GetNbinsY();
    if (doFold == 1 || doFold == 3) {
      for (int iy = 0; iy <= nY + 1; ++iy) {
        h->SetBinContent(nX, iy, h->GetBinContent(nX, iy) + h->GetBinContent(nX+1, iy));
        h->SetBinError(nX, iy, std::hypot(h->GetBinError(nX, iy), h->GetBinError(nX+1, iy)));
        h->SetBinContent(nX+1, iy, 0); h->SetBinError(nX+1, iy, 0);
      }
      for (int ix = 0; ix <= nX + 1; ++ix) {
        h->SetBinContent(ix, nY, h->GetBinContent(ix, nY) + h->GetBinContent(ix, nY+1));
        h->SetBinError(ix, nY, std::hypot(h->GetBinError(ix, nY), h->GetBinError(ix, nY+1)));
        h->SetBinContent(ix, nY+1, 0); h->SetBinError(ix, nY+1, 0);
      }
    }
    if (doFold == 2 || doFold == 3) {
      for (int iy = 0; iy <= nY + 1; ++iy) {
        h->SetBinContent(1, iy, h->GetBinContent(1, iy) + h->GetBinContent(0, iy));
        h->SetBinError(1, iy, std::hypot(h->GetBinError(1, iy), h->GetBinError(0, iy)));
        h->SetBinContent(0, iy, 0); h->SetBinError(0, iy, 0);
      }
      for (int ix = 0; ix <= nX + 1; ++ix) {
        h->SetBinContent(ix, 1, h->GetBinContent(ix, 1) + h->GetBinContent(ix, 0));
        h->SetBinError(ix, 1, std::hypot(h->GetBinError(ix, 1), h->GetBinError(ix, 0)));
        h->SetBinContent(ix, 0, 0); h->SetBinError(ix, 0, 0);
      }
    }
  }
}

TH1D* UnrollHistogram(TH2D* h2) {
  int nBinsX = h2->GetNbinsX(), nBinsY = h2->GetNbinsY();
  int totalBins = nBinsX * nBinsY;
  TString name = h2->GetName(), title = h2->GetTitle();
  h2->SetName(name + "_old");
  TH1D* h1 = new TH1D(name, title, totalBins, 0.5, totalBins + 0.5);
  h1->Sumw2();
  int i1D = 1;
  for (int iy = 1; iy <= nBinsY; ++iy)
    for (int ix = 1; ix <= nBinsX; ++ix) {
      h1->SetBinContent(i1D, h2->GetBinContent(ix, iy));
      h1->SetBinError(i1D, h2->GetBinError(ix, iy));
      i1D++;
    }
  return h1;
}

#include "TMVA/RReader.hxx"
#include "TMVA/RInferenceUtils.hxx"

ROOT::RDF::RNode SafeDefine(ROOT::RDF::RNode df, std::string name, std::string expr) {
    auto colNames = df.GetColumnNames();
    if (std::find(colNames.begin(), colNames.end(), name) == colNames.end())
        return df.Define(name, expr);
    return df;
}

std::vector<std::string> getBranchNames(TTree* tree) {
  std::vector<std::string> names;
  TObjArray* branches = tree->GetListOfBranches();
  for (int i = 0; i < branches->GetEntries(); ++i)
    names.push_back(branches->At(i)->GetName());
  return names;
}
"""
        with open(f"{self._scripts_run_folder}/library_utils.cpp", "w") as f:
            f.write(cpp_code_library)
        with open(f"{self._scripts_run_folder}/library_utils.h", "w") as f:
            f.write(cpp_code_library_header)

    # _____________________________________________________________________________
    #  nuisances — mkShapes-specific overrides
    # _____________________________________________________________________________

    def _build_cpp_input_files_extra(self, file_path, sampleName, tree_name="Events"):
        """Add friend trees for suffix nuisances + getBranchNames calls."""
        out = "    auto nominalBranches = getBranchNames(nominal);\n\n"
        out += "//  Variations input files (if any)\n"
        for nuisanceName, nuisance in self._nuisances.items():
            #print(nuisanceName) ### mcortino debug
            if (sampleName in nuisance['samples'] and nuisance['type'] == 'shape'
                    and nuisance.get('kind') == 'suffix'):
                mUp, mDo = nuisance['mapUp'], nuisance['mapDown']
                #print('---', mUp, mDo) ### mcortino debug
                out += f'    auto* friend_{mUp} = new TChain("{tree_name}");\n' ### mcortino: {self._treeName} --> {tree_name}
                out += f'    auto* friend_{mDo} = new TChain("{tree_name}");\n' ### mcortino: {self._treeName} --> {tree_name}
                for file_name in file_path:
                    fname = file_name[file_name.rfind("/")+1:]
                    out += f'    friend_{mUp}->Add("{nuisance["folderUp"]}/{fname}");\n'
                    out += f'    friend_{mDo}->Add("{nuisance["folderDown"]}/{fname}");\n'
                out += f'    nominal->AddFriend(friend_{mUp}, "{mUp}");\n'
                out += f'    nominal->AddFriend(friend_{mDo}, "{mDo}");\n'
                out += f'    auto varBranches_{mUp} = getBranchNames(friend_{mUp});\n'
                out += f'    auto varBranches_{mDo} = getBranchNames(friend_{mDo});\n'
        return out

    def _build_cpp_variations_alt_trees(self, sampleName, file_path, weight):
        out = "\n"
        first_time_suffix = 0
        for nuisanceName, nuisance in self._nuisances.items():
            if (sampleName in nuisance['samples'] and nuisance['type'] == 'shape'
                    and nuisance.get('kind') == 'suffix'):
                if first_time_suffix == 0:
                    out += '    int suffix_size = 0;\n'
                    first_time_suffix = 1
                for (mapName, direction) in [(nuisance['mapUp'], f'"{TAG_UP}"'), (nuisance['mapDown'], f'"{TAG_DOWN}"')]:
                    out += f'    suffix_size = {len(mapName) + 1};\n'
                    out += f'    for (const auto& branch : varBranches_{mapName}) {{\n'
                    out += f'      if (int(branch.size()) >= (suffix_size-1) && branch.compare(branch.size() - (suffix_size-1), (suffix_size-1), "{mapName}") == 0 ) {{\n'
                    out += f'        std::string nomCol = branch.substr(0, branch.size() - suffix_size);\n'
                    out += f'        if (std::find(nominalBranches.begin(), nominalBranches.end(), nomCol) != nominalBranches.end()) {{\n'
                    out += f'          std::string col_type = varied_df.GetColumnType(nomCol);\n'
                    out += f'          std::string expression = (col_type.find("RVec") != std::string::npos)\n'
                    out += f'            ? ("ROOT::RVec<" + col_type + ">{{" + branch + "}}")\n'
                    out += f'            : ("ROOT::RVec<" + col_type + ">{{static_cast<" + col_type + ">(" + branch + ")}}");\n'
                    out += f'          varied_df = varied_df.Vary(nomCol, expression, {{{direction}}}, "{nuisance["name"]}");\n'
                    out += f'        }};\n      }};\n    }};\n\n'
        return out


    def _build_cpp_variations_weights(self, sampleName, weight):
        out = "//  Register the variations\n"
        for nuisanceName, nuisance in self._nuisances.items():
            if (sampleName in nuisance['samples'] and nuisance['type'] == 'shape'
                    and nuisance.get('kind') == 'weight'):
                samples_entry = nuisance['samples'][sampleName]
                if len(samples_entry) == 3:
                    if samples_entry[0] in weight:
                        var_up = weight.replace(samples_entry[0], samples_entry[1])
                        var_do = weight.replace(samples_entry[0], samples_entry[2])
                        out += f'    varied_df = varied_df.Vary("my_sample_weight", "ROOT::RVecD{{{var_up},{var_do}}}", {{"{TAG_UP}", "{TAG_DOWN}"}}, "{nuisance["name"]}");\n'
                else:
                    var_up = f"(my_sample_weight) * {samples_entry[0]}"
                    var_do = f"(my_sample_weight) * {samples_entry[1]}"
                    out += f'    varied_df = varied_df.Vary("my_sample_weight", "ROOT::RVecD{{{var_up},{var_do}}}", {{"{TAG_UP}", "{TAG_DOWN}"}}, "{nuisance["name"]}");\n'
        return out

    def _build_cpp_variations_special(self, sampleName):
        out = ""
        for nuisanceName, nuisance in self._nuisances.items():
            if (sampleName in nuisance['samples'] and nuisance['type'] == 'shape'
                    and nuisance.get('kind') == 'suffixFly'):
                with open(nuisance["code_in_main"], "r") as fh:
                    out += fh.read() + "\n"
        return out

    def _build_cpp_code_special_nuisances(self, sampleName):
        out = ""
        seen = set()
        for nuisanceName, nuisance in self._nuisances.items():
            if (sampleName in nuisance['samples'] and nuisance['type'] == 'shape'
                    and nuisance.get('kind') == 'suffixFly'):
                if nuisance["external"] not in seen:
                    with open(nuisance["external"], "r") as fh:
                        out += fh.read() + "\n"
                    seen.add(nuisance["external"])
        return out

    def _build_cpp_maps_special(self, sampleName):
        out = ""
        for nuisanceName, nuisance in self._nuisances.items():
            if (sampleName in nuisance['samples'] and nuisance['type'] == 'shape'
                    and nuisance.get('kind') == 'suffixFly'):
                out += f'    map_nuisances_name_for_averaged_histograms["{nuisance["pattern"]}"] = std::pair<std::string, int>("{nuisance["name"]}", {nuisance["numVariations"]});\n'
        return out

    # _____________________________________________________________________________
    #  Cuts + histograms — the mkShapes-specific output fragment
    # _____________________________________________________________________________

    def _build_cpp_cuts_and_output(self, sampleName):
        """Returns C++ code for filter nodes + Histo1D/Histo2D bookings."""
        out = ""
        subsamples = getattr(self, '_current_subsamples', {})

        for cutName, cut in self._cuts.items():
            list_cuts = {}
            if isinstance(cut, dict):
                expression = cut['expr']
                categories = cut.get('categories', {})
                if categories:
                    list_cuts[cutName] = "(" + expression + ")"
                    for cat_name, cat in categories.items():
                        list_cuts[cutName + "_" + cat_name] = "(" + cat + ")"
                else:
                    list_cuts[cutName] = expression
            else:
                list_cuts[cutName] = cut

            mother_cut_name = ""
            for icut, (this_cutName, this_cut) in enumerate(list_cuts.items()):
                if len(list_cuts) > 1:
                    if icut == 0:
                        out += f'    auto node_{this_cutName} = current_node.Filter("{this_cut}", "{this_cutName}");\n'
                        mother_cut_name = this_cutName
                    else:
                        out += f'    auto node_{this_cutName} = node_{mother_cut_name}.Filter("{this_cut}", "{this_cutName}");\n'
                else:
                    out += f'    auto node_{this_cutName} = current_node.Filter("{this_cut}", "{this_cutName}");\n'

                if subsamples:
                    for sub_name, sub_cut_name in subsamples.items():
                        if 'cut' in sub_cut_name:
                            out += f'    ROOT::RDF::RNode node_{this_cutName}___{sampleName}_{sub_name} = node_{this_cutName}.Filter("{sub_cut_name["cut"]}", "cut_{sampleName}_{sub_name}");\n'
                            if 'weight' in sub_cut_name:
                                out += f'    node_{this_cutName}___{sampleName}_{sub_name} = SafeDefine(node_{this_cutName}___{sampleName}_{sub_name}, "my_sample_weight_{sampleName}_{sub_name}", "(my_sample_weight)*{sub_cut_name["weight"]}");\n'
                        elif 'weight' in sub_cut_name:
                            out += f'    node_{this_cutName} = SafeDefine(node_{this_cutName}, "my_sample_weight_{sampleName}_{sub_name}", "(my_sample_weight)*{sub_cut_name["weight"]}");\n'

                for variableName, variable in self._variables.items():
                    variable_range = variable['range']
                    model_str = ""
                    if isinstance(variable_range, list):
                        bins = len(variable_range) - 1
                        edges_str = ", ".join(map(str, variable_range))
                        model_str = f'{bins}, (const double[]){{{edges_str}}}'
                    elif isinstance(variable_range, tuple):
                        if len(variable_range) == 3:
                            bins, v_min, v_max = variable_range
                            model_str = f'{bins}, {v_min}, {v_max}'
                        elif len(variable_range) == 6:
                            model_str = ", ".join(map(str, variable_range))
                        elif len(variable_range) == 2 and isinstance(variable_range[0], list):
                            nx = len(variable_range[0]) - 1
                            ny = len(variable_range[1]) - 1
                            ex = ", ".join(map(str, variable_range[0]))
                            ey = ", ".join(map(str, variable_range[1]))
                            model_str = f'{nx}, (const double[]){{{ex}}}, {ny}, (const double[]){{{ey}}}'

                    if 'is2d' in variable and variable['is2d'] == 1:
                        v_x, v_y = variable['name'].split(':')
                        histo_call = f'Histo2D({{"h_{variableName}", "{variableName}", {model_str}}}, "{v_x}", "{v_y}", "my_sample_weight")'
                        out += f'    hist_map_2D["{this_cutName}"].push_back(node_{this_cutName}.{histo_call});\n'
                    else:
                        histo_call = f'Histo1D<float>({{"h_{variableName}", "{variableName}", {model_str}}}, "{variableName}", "my_sample_weight")'
                        out += f'    hist_map_1D["{this_cutName}"].push_back(node_{this_cutName}.{histo_call});\n'

                    if subsamples:
                        for sub_name, sub_cut_name in subsamples.items():
                            name_weight = "my_sample_weight"
                            if 'weight' in sub_cut_name:
                                name_weight = f"my_sample_weight_{sampleName}_{sub_name}"
                            if 'cut' in sub_cut_name:
                                if 'is2d' in variable and variable['is2d'] == 1:
                                    v_x, v_y = variable['name'].split(':')
                                    histo_call = f'Histo2D({{"h_{variableName}", "{variableName}", {model_str}}}, "{v_x}", "{v_y}", "{name_weight}")'
                                    out += f'    hist_map_2D["{this_cutName}___{sampleName}_{sub_name}"].push_back(node_{this_cutName}___{sampleName}_{sub_name}.{histo_call});\n'
                                else:
                                    histo_call = f'Histo1D<float>({{"h_{variableName}", "{variableName}", {model_str}}}, "{variableName}", "{name_weight}")'
                                    out += f'    hist_map_1D["{this_cutName}___{sampleName}_{sub_name}"].push_back(node_{this_cutName}___{sampleName}_{sub_name}.{histo_call});\n'

        # variable lists for fold
        out += '\n    std::vector<std::string> list_of_variables_1D;\n'
        out += '    std::vector<std::string> list_of_variables_2D;\n'
        out += '    std::vector<int> list_of_variables_fold_1D;\n'
        out += '    std::vector<int> list_of_variables_fold_2D;\n'
        for variableName, variable in self._variables.items():
            if 'is2d' in variable and variable['is2d'] == 1:
                out += f'    list_of_variables_2D.push_back("{variableName}");\n'
                out += f'    list_of_variables_fold_2D.push_back({variable.get("fold", 0)});\n'
            else:
                out += f'    list_of_variables_1D.push_back("{variableName}");\n'
                out += f'    list_of_variables_fold_1D.push_back({variable.get("fold", 0)});\n'
        return out

    # _____________________________________________________________________________
    #  create_cpp_source_list_of_files
    # _____________________________________________________________________________

    def create_cpp_source_list_of_files(self, output_name, tree_name, file_path,
                                         output_root_file_name, sampleName, weight,
                                         aliases_to_be_defined, functions_to_be_defined,
                                         subsamples, aliases_ordered=None,
                                         range_start=None, range_end=None):
        self._current_subsamples = subsamples
        frags = self._build_cpp_common(tree_name, file_path, sampleName, weight,
                                       aliases_ordered or {}, range_start, range_end)

        # The big C++ template — same as the original mkShapes
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
#include <string>
#include <typeinfo>
#include <regex>
#include "TMVA/RReader.hxx"
#include "TMVA/RInferenceUtils.hxx"

{frags['lines_to_add']}

{frags['code_of_function_to_include']}

{frags['lines_to_global']}

{frags['code_for_special_nuisances']}

int main() {{
    ROOT::DisableImplicitMT();

    {frags['lines_to_process']}

    {frags['define_input_files_logic']}

    std::map<std::string, std::vector<ROOT::RDF::RResultPtr<TH1D>>> hist_map_1D;
    std::map<std::string, std::vector<ROOT::RDF::RResultPtr<TH2D>>> hist_map_2D;

    {frags['register_variations_alt_trees']}

    {frags['register_variations_special']}

    {frags['define_aliases']}

    {frags['define_weights']}

    {frags['register_variations_weights']}

    {frags['booking_logic']}

    {frags['define_variables_logic']}

    // Variations + RunGraphs
    std::vector<ROOT::RDF::Experimental::RResultMap<TH1D>> results_1D_variations;
    std::vector<ROOT::RDF::Experimental::RResultMap<TH2D>> results_2D_variations;
    std::vector<ROOT::RDF::RResultHandle> all_booking_1D;
    std::vector<ROOT::RDF::RResultHandle> all_booking_2D;

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

    if (!hist_map_1D.empty()) ROOT::RDF::RunGraphs(all_booking_1D);
    if (!hist_map_2D.empty()) ROOT::RDF::RunGraphs(all_booking_2D);

    std::map<std::string, std::pair<std::string, int>> map_nuisances_name_for_averaged_histograms;
    {frags['register_maps_special']}

    TFile out_file("{self._outputDir}/{output_root_file_name}", "RECREATE");
    int big_loop = 0;
    for (auto& [cut_label, h_list] : hist_map_1D) {{
      std::string current_cut = cut_label;
      std::string sub_name = "";
      size_t pos_triple = cut_label.find("___");
      if (pos_triple != std::string::npos) {{
        std::string rest = cut_label.substr(pos_triple + 3);
        std::string sN = "{sampleName}";
        size_t pos = rest.find(sN);
        sub_name = rest.substr(pos + sN.length() + 1);
        current_cut = cut_label.substr(0, pos_triple);
      }}
      TDirectory *dir = out_file.GetDirectory(current_cut.c_str());
      if (!dir) dir = out_file.mkdir(current_cut.c_str());
      for (size_t ivar = 0; ivar < h_list.size(); ivar++) {{
        dir->cd();
        TDirectory *subdir = out_file.GetDirectory((current_cut+"/"+list_of_variables_1D.at(ivar)).c_str());
        if (!subdir) subdir = out_file.mkdir((current_cut+"/"+list_of_variables_1D.at(ivar)).c_str());
        subdir->cd();
        auto all_histos = results_1D_variations.at(big_loop++);
        std::map<std::string, TH1D> map_nuisances_average_histogram;
        for (auto& [name, histo] : all_histos) {{
          std::string temp_name;
          if (name == "nominal") {{
            temp_name = (sub_name == "") ? "histo_{sampleName}" : "histo_{sampleName}_" + sub_name;
          }} else {{
            temp_name = name;
            temp_name.erase(std::remove(temp_name.begin(), temp_name.end(), ':'), temp_name.end());
            temp_name = (sub_name == "") ? "histo_{sampleName}_" + temp_name : "histo_{sampleName}_" + sub_name + "_" + temp_name;
          }}
          gDirectory = subdir;
          histo->SetName(temp_name.c_str());
          if (list_of_variables_fold_1D.at(ivar) != 0) FoldHistogram(histo.get(), list_of_variables_fold_1D.at(ivar));
          bool found_avg = false;
          for (auto const& [sp, nn] : map_nuisances_name_for_averaged_histograms) {{
            std::string nu = nn.first; int nv = nn.second;
            std::regex rup("histo_.*__NORM__.*_" + sp + "{TAG_UP}");
            std::regex rdo("histo_.*__NORM__.*_" + sp + "{TAG_DOWN}");
            if (std::regex_search(temp_name, rup) || std::regex_search(temp_name, rdo)) {{
              std::string hn = std::regex_replace(std::regex_replace(temp_name, std::regex("__NORM__[0-9]+_"), ""), std::regex(sp), nu);
              histo->SetName(hn.c_str());
              bool isUp = std::regex_search(temp_name, rup);
              std::string key = isUp ? sp+"{TAG_UP}" : sp+"{TAG_DOWN}"; std::string nukey = isUp ? nu+"{TAG_UP}" : nu+"{TAG_DOWN}";
              if (auto it = map_nuisances_average_histogram.find(key); it != map_nuisances_average_histogram.end())
                it->second.Add(histo.get(), 1./nv);
              else {{ map_nuisances_average_histogram[nukey] = (*histo); found_avg = true; break; }}
            }}
          }}
          if (!found_avg) histo->Write();
          for (auto const& [n2, th] : map_nuisances_average_histogram) th.Write();
        }}
      }}
      out_file.cd();
    }}
    big_loop = 0;
    for (auto& [cut_label, h_list] : hist_map_2D) {{
      std::string current_cut = cut_label;
      std::string sub_name = "";
      size_t pos_triple = cut_label.find("___");
      if (pos_triple != std::string::npos) {{
        std::string rest = cut_label.substr(pos_triple + 3);
        std::string sN = "{sampleName}";
        size_t pos = cut_label.find(sN);
        sub_name = cut_label.substr(pos + sN.length() + 1);
        current_cut = cut_label.substr(0, pos_triple);
      }}
      TDirectory *dir = out_file.GetDirectory(current_cut.c_str());
      if (!dir) dir = out_file.mkdir(current_cut.c_str());
      for (size_t ivar = 0; ivar < h_list.size(); ivar++) {{
        dir->cd();
        TDirectory *subdir = out_file.GetDirectory((current_cut+"/"+list_of_variables_2D.at(ivar)).c_str());
        if (!subdir) subdir = out_file.mkdir((current_cut+"/"+list_of_variables_2D.at(ivar)).c_str());
        subdir->cd();
        auto all_histos = results_2D_variations.at(big_loop++);
        for (auto& [name, histo] : all_histos) {{
          std::string temp_name;
          if (name == "nominal") {{
            temp_name = (sub_name == "") ? "histo_{sampleName}" : "histo_{sampleName}_" + sub_name;
          }} else {{
            temp_name = name;
            temp_name.erase(std::remove(temp_name.begin(), temp_name.end(), ':'), temp_name.end());
            temp_name = (sub_name == "") ? "histo_{sampleName}_" + temp_name : "histo_{sampleName}_" + sub_name + "_" + temp_name;
          }}
          gDirectory = subdir;
          histo->SetName(temp_name.c_str());
          if (list_of_variables_fold_2D.at(ivar) != 0) FoldHistogram(histo.get(), list_of_variables_fold_2D.at(ivar));
          auto new_histo = UnrollHistogram(dynamic_cast<TH2D*>(histo.get()));
          new_histo->Write();
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
    #  generate_makefile
    # _____________________________________________________________________________

    def generate_makefile(self, file_paths, dictionary_of_files_divided_by_sample, makefile_name="Makefile"):
        cpp_flags = "$(shell root-config --cflags)"
        ld_flags  = (f"$(shell root-config --libs) -lTMVA -lXMLIO "
                     f"-I/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/include "
                     f"-L/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/lib "
                     f"-lcorrectionlib -Wl,-rpath,'$$ORIGIN:$$ORIGIN/../../'")
        targets = [os.path.splitext(f)[0] for f in file_paths]
        with open(makefile_name, "w") as f:
            f.write("# Generated Makefile\nCXX = g++\n")
            f.write(f"CXXFLAGS = -O2 -Wall -fPIC {cpp_flags} -I{self._scripts_run_folder}/\n")
            f.write(f"LDFLAGS = {ld_flags}\n\nLIB_NAME = liblibrary_utils.so\n\n")
            all_targets = " ".join(targets)
            f.write(f"all: $(LIB_NAME) {all_targets}\n\n")
            for sample_key, sample_files in dictionary_of_files_divided_by_sample.items():
                st = " ".join(os.path.splitext(x)[0] for x in sample_files)
                f.write(f"{sample_key}: $(LIB_NAME) {st}\n\n")
            f.write(f"$(LIB_NAME): {self._scripts_run_folder}/library_utils.cpp {self._scripts_run_folder}/library_utils.h\n")
            f.write(f"\t$(CXX) $(CXXFLAGS) -shared -o {self._scripts_run_folder}/$@ {self._scripts_run_folder}/library_utils.cpp\n")
            f.write("%: %.cpp $(LIB_NAME)\n")
            f.write(f"\t$(CXX) $(CXXFLAGS) $< -o $@ -L{self._scripts_run_folder}/ -llibrary_utils $(LDFLAGS)\n\n")
            f.write("clean:\n")
            f.write(f"\trm -f $(LIB_NAME) {all_targets}\n")

    # _____________________________________________________________________________
    #  makeNominals
    # _____________________________________________________________________________

    def makeNominals(self):
        print("======================\n==== makeNominals ====\n======================")

        os.makedirs(self._outputDir, exist_ok=True)
        ROOT.TH1.SetDefaultSumw2(True)
        os.makedirs(self._scripts_run_folder, exist_ok=True)
        self.create_library_cpp()

        self._pre_count_events()

        list_of_files_to_compile = []
        list_of_files_divided    = {}

        for sampleName, sample in self._samples.items():
            os.makedirs(self._scripts_run_folder + "/" + sampleName, exist_ok=True)
            list_of_files_divided[sampleName] = []

            global_weight   = sample.get('weight', "1.")
            weights_subname = sample.get('weights', {})
            if 'isData' not in sample:
                global_weight = f"({global_weight}) * {self._lumi}"

            _, _, aliases_ordered = self._build_aliases_ordered(sampleName)

            for subname, list_root_files in self._normalize_sample_name(sample['name']).items():
                os.makedirs(self._scripts_run_folder + "/" + sampleName + "/" + subname, exist_ok=True)

                weight = (f"({global_weight}) * ({weights_subname[subname]})"
                          if subname in weights_subname else global_weight)

                chunks = self._get_chunks(sample, list_root_files)

                subsamples_dict = sample.get('subsamples', {})

                for i, chunk in enumerate(chunks):
                    name_code = (self._scripts_run_folder + "/" + sampleName + "/" + subname
                                 + "/" + self.job_prefix() + "_" + sampleName + "_" + subname + "_" + str(i))
                    out_file  = self.output_file_name(sampleName, subname, i)

                    self.create_cpp_source_list_of_files(
                        name_code, "Events", chunk['files'], out_file, sampleName, weight,
                        {}, {}, subsamples_dict, aliases_ordered,
                        range_start=chunk['range_start'], range_end=chunk['range_end'])

                    list_of_files_to_compile.append(name_code)
                    list_of_files_divided[sampleName].append(name_code)

        self.generate_makefile(list_of_files_to_compile, list_of_files_divided)
        print("Do not compile ... If you want to compile run the compilation by hand!")
        subprocess.run(["make", "clean"])

    # _____________________________________________________________________________
    #  haddBatch — submit a single condor job to hadd all root files
    # _____________________________________________________________________________

    def haddBatch(self, tag, haddfast_path=None):
        if haddfast_path is None:
            haddfast_path = os.path.join(_SCRIPTS_DIR, "haddfast.py")
        super().haddBatch(
            tag          = tag,
            output_dir   = self._outputDir,
            name_prefix  = "root_file___",
            # same name as --hadd: mkPlot and mkDatacards look for this file
            output_name  = "root_file_joined.root",
            hadd_tool_path = haddfast_path,
            # -C purges the duplicate key cycles haddfast leaves behind. Its
            # merge step calls TH1::Write() without kOverwrite, so every
            # histogram that is *added* to an existing one gains a cycle:
            # ;1 holds the first chunk alone and the highest cycle the running
            # sum. TFile::Get() resolves to the highest cycle, so mkPlot and
            # mkDatacards read the right one, but any consumer that iterates
            # GetListOfKeys() sees all of them, in descending cycle order.
            hadd_tool_args = "-f -j 8 -E -C",
        )

    # _____________________________________________________________________________
    #  submitBatch hooks
    # _____________________________________________________________________________

    def _prepare_submit_batch(self, root_cflags, root_libs):
        """Compile liblibrary_utils.so locally if missing."""
        abs_library = os.path.join(self._scripts_run_folder, "liblibrary_utils.so")
        abs_lib_cpp = os.path.join(self._scripts_run_folder, "library_utils.cpp")
        if not os.path.exists(abs_library):
            print("  -> liblibrary_utils.so not found, compiling it locally...")
            ret = os.system(f"g++ -O2 -Wall -fPIC {root_cflags} -shared -o {abs_library} {abs_lib_cpp} {root_libs}")
            if ret == 0:
                print(f"  -> liblibrary_utils.so compiled: {abs_library}")
            else:
                print("  -> ERROR: failed to compile liblibrary_utils.so")
        self._abs_library = abs_library
        self._abs_lib_cpp = abs_lib_cpp
        self._abs_lib_h   = os.path.join(self._scripts_run_folder, "library_utils.h")

    def _extra_transfer_files_for_submit(self, precompiled, x509_proxy, extra_files_str):
        if precompiled:
            return f"{self._abs_library}, {x509_proxy}{extra_files_str}"
        return f"{self._abs_lib_cpp}, {self._abs_lib_h}, {self._abs_library}, {x509_proxy}{extra_files_str}"

    def _build_bash_precompiled(self, job_name, out_file, sampleName, subname,
                                proxy_lines, submission_dir, root_cflags, root_libs):
        cp_cmd = self._build_cp_cmd(out_file, self._outputDir, src_prefix=self._outputDir)
        return f"""#!/bin/bash
set -e
echo "Job started at $(date)"; echo "Running on node $(hostname)"
source /cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/setup.sh
{proxy_lines}mkdir -p {sampleName}/{subname}/
mv {job_name} {sampleName}/{subname}/
mkdir -p {self._outputDir}
./{sampleName}/{subname}/{job_name}
{cp_cmd}
echo "Done at $(date)"; ls -lh
"""

    def _build_bash_compile_on_node(self, job_name, out_file, sampleName, subname,
                                    proxy_lines, submission_dir, root_cflags, root_libs):
        cp_cmd = self._build_cp_cmd(out_file, self._outputDir, src_prefix=self._outputDir)
        return f"""#!/bin/bash
set -e
echo "Job started at $(date)"; echo "Running on node $(hostname)"
source /cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/setup.sh
{proxy_lines}echo "Compiling..."
g++ -O2 -Wall -fPIC {root_cflags} -I. \\
    {job_name}.cpp -o {job_name} -L. -llibrary_utils \\
    -lTMVA -lXMLIO \\
    -I/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/include \\
    -L/cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/lib \\
    -lcorrectionlib {root_libs} -Wl,-rpath,'$ORIGIN'
echo "Compilation done at $(date)"
mkdir -p {sampleName}/{subname}/
mv liblibrary_utils.so {sampleName}/{subname}/
mv {job_name} {sampleName}/{subname}/
export LD_LIBRARY_PATH=$(pwd)/{sampleName}/{subname}/:$LD_LIBRARY_PATH
mkdir -p {self._outputDir}
./{sampleName}/{subname}/{job_name}
{cp_cmd}
echo "Done at $(date)"; ls -lh
"""


# -----------------------------------------------------------------------
#  Main
# -----------------------------------------------------------------------

if __name__ == '__main__':
    sys.argv = argv

    header = """
         --------------------------------------------------------------------------------------------------
         '                ___|   |                               \\  |         |                           '
         '              \\___ \\   __ \\    _` |  __ \\    _ \\      |\\/ |   _` |  |  /   _ \\   __|            '
         '                    |  | | |  (   |  |   |   __/      |   |  (   |    <    __/  |               '
         '              _____/  _| |_| \\__,_|  .__/  \\___|     _|  _| \\__,_| _|\\_\\ \\___| _|               '
         --------------------------------------------------------------------------------------------------
    """
    print(header)

    parser = defaultParser()
    parser.add_argument("--parallelCompile", action='store_true', dest="parallelCompile")
    parser.add_argument("--localCompile",    action='store_true', dest="localCompile")
    parser.add_argument("--submitBatch",     action='store_true', dest="submitBatch")
    parser.add_argument("--hadd",            action='store_true', dest="haddRootFiles")
    parser.add_argument("--haddBatch",       action='store_true', dest="haddBatch",
                        help="Submit a single condor job to hadd all root files")
    parser.add_argument("--checkBatch",      action='store_true', dest="checkBatch")
    parser.add_argument("--resubmit",        action='store_true', dest="resubmit")
    parser.add_argument("--silentMode",      action='store_true', dest="silentMode")

    opt = parser.parse_args()
    loadDefaultOptions(parser, opt.pycfg)
    opt = parser.parse_args()

    cfg       = load_config_files(opt)
    samples   = cfg['samples']
    variables = cfg['variables']
    cuts      = cfg['cuts']
    supercut  = cfg['supercut']
    nuisances = cfg['nuisances']
    aliases   = cfg['aliases']

    factory = ShapeFactory()
    factory.setValues(
        opt.outputDir, variables, cuts, supercut, samples, nuisances, aliases, opt.lumi,
        scripts_run_folder    = getattr(opt, 'scripts_run_folder',    None),
        script_batch_location = getattr(opt, 'script_batch_location', None),
        extra_transfer_files  = getattr(opt, 'extra_transfer_files',  None),
        jobs_config           = getattr(opt, 'jobs_config',           None)
    )
    factory.setConditions(opt.silentMode)

    if not any([opt.submitBatch, opt.haddRootFiles, opt.haddBatch, opt.checkBatch,
                opt.localCompile, opt.parallelCompile]):
        factory.makeNominals()

    if opt.localCompile:
        subprocess.run(["make", "-j4"])
    if opt.parallelCompile:
        factory.parallelCompile()
    if opt.submitBatch:
        factory.submitBatch()
    if opt.haddBatch:
        tag = getattr(opt, 'tag', 'merged')
        haddfast_path = getattr(opt, 'haddfastPath', None)
        factory.haddBatch(tag, haddfast_path=haddfast_path)

    if opt.haddRootFiles:
        all_done, missing = factory._check_all_done()
        if not all_done:
            print(f"ERROR: {len(missing)} output files missing — cannot hadd incomplete results.")
            for f in missing:
                print(f"  {f}")
            print("Run --checkBatch to get the full picture, --resubmit to resubmit missing jobs.")
        else:
            # the shell glob would also pick up the per-job files of a sample
            # that samples.py no longer declares, and hadd would add them to
            # the histograms of the sample that replaced it
            inputs, orphans = factory.selectOutputFiles(factory._outputDir, "root_file___")
            factory.reportOrphans(orphans)
            if not inputs:
                print("ERROR: no file belongs to a sample declared in samples.py")
            else:
                joined = os.path.join(factory._outputDir, "root_file_joined.root")
                print(f"  -> hadd of {len(inputs)} files into {joined}")
                subprocess.run(["hadd", "-j", "8", "-f", joined] + inputs)
    if opt.checkBatch:
        factory.checkBatch(resubmit=opt.resubmit)

    print("\n\n I'm done ...\n\n")
