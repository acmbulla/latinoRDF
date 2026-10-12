#!/usr/bin/env python3

import sys
argv = sys.argv
sys.argv = argv[:1]

import argparse

import os
import re

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
        default="merge_configuration.py",
    )

    return parser


#
# load configurations into the parser, according to configuration.py
#

def loadDefaultOptions(parser, pycfg=None, quiet=False):

    if pycfg != None and os.path.exists(pycfg):
        handle = open(pycfg,'r')
        local_variables = {}
        exec(handle.read(),local_variables)
        handle.close()
        local_variables = {k: v for k, v in local_variables.items() if not (k.startswith('__') and k.endswith('__'))}
        for opt_name, opt_value in local_variables.items():
          parser.add_argument('--' + opt_name, default=opt_value)
        return
    else:
        return




# ----------------------------------------------------- Merge Years: e.g. Run 2, Run 3, ... --------------------------------------


class MergerFactory:

    # _____________________________________________________________________________
    def __init__(self):
        self._fileIn = None

    # _____________________________________________________________________________
    def merge(
        self,
        tag,
        variables,
        cuts,
        samples,
        nuisances,
        foldersToMerge,
        foldersToMergeNuisancesFiles,
        outputDir=None,
        plotDir=None,
    ):


       #
       # expand cuts, unrolling the categories
       #

       list_cuts = {}
       for cutName, cut in cuts.items():

         if isinstance(cut, dict):
           expression = cut['expr']
           categories = cut.get('categories', {})
           if categories :
             list_cuts[ cutName ] = "(" + expression + ")"
             for category_name, category in categories.items():
               list_cuts[ cutName + "_" + category_name] = "(" + expression + ") && (" + category + ")"
           else :
             list_cuts[ cutName ] = expression
         else :
           list_cuts[ cutName ] = cut

       self._cuts = list_cuts


       new_list_of_nuisances = nuisances


       #
       # Get all the samples from all the years. This is needed for the definition of the nuisances
       #
       foldersToMergeSamplesFiles = {
           k : os.path.join(p["folder"], 'samples.py') for k, p in foldersToMerge.items()
           }


       #
       # get all the nuisances from all the years
       #
       all_nuisances = {}

       for folderHR, year_nuisances in foldersToMergeNuisancesFiles.items() :
         year_samples = foldersToMergeSamplesFiles[folderHR]

         print ("year_nuisances = ", year_nuisances)
         print ("year_samples   = ", year_samples)

         with open(year_samples, 'r') as f_samples:
           handle_samples = f_samples.read()
         with open(year_nuisances, 'r') as f_nuisances:
           handle_nuisances = f_nuisances.read()

         handle_unique = handle_samples + "\n# --- End of Samples / Start of Nuisances ---\n" + handle_nuisances

         exec_context = {'os': os,
                         'samples': {},
                         'nuisances': {}
                        }
         exec(handle_unique, exec_context)
         nuisances = exec_context.get('nuisances', {})
         nuisances = {k: v for k, v in nuisances.items() if not (k.startswith('__') and k.endswith('__'))}

         all_nuisances[folderHR] = nuisances

         print (" all_nuisances[", folderHR, "] = ", all_nuisances[folderHR])


       #
       # get all the root files with the histograms from the different years
       #
       _out_base = "histos_" + tag + ".root"
       if outputDir:
           os.makedirs(outputDir, exist_ok=True)
           new_root_file_name = os.path.join(outputDir, _out_base)
       else:
           new_root_file_name = _out_base
       print(f"Output root file: {new_root_file_name}")
       if plotDir:
           os.makedirs(plotDir, exist_ok=True)
           print(f"Plot directory:   {plotDir}")
       rootFileNew = ROOT.TFile.Open( new_root_file_name, "RECREATE")

       histograms = {}

       for folderHR, folder in foldersToMerge.items():

         # [FEATURE] sampleRename: optional dict {oldName: newName} per year
         # Used to rename sample histograms when sample names differ between years
         # Example: {'DY_2024': 'DY', 'top_v2': 'top'}
         sample_rename = folder.get('sampleRename', {})
         # build histo-level rename map: histo_oldName -> histo_newName
         # also handles nuisance variations: histo_oldName_<nuis>up -> histo_newName_<nuis>up
         histo_rename = {f"histo_{old}": f"histo_{new}" for old, new in sample_rename.items()}

         old_root_file_name = os.path.join(folder["folder"], "rootFile", folder["tag"], "root_file_joined.root")
         import tempfile
         tmp_dir = tempfile.gettempdir()
         new_root_file_name = os.path.join(tmp_dir, "year_" + folderHR + "_histos_" + folder["tag"] + ".root")
         # use xrdcp for EOS paths, plain cp otherwise
         if "/eos/" in old_root_file_name:
             import re as _re
             eos_normalized = _re.sub(r'^/eos/home-([a-z])/', lambda m: f"/eos/user/{m.group(1)}/", old_root_file_name)
             ret = os.system(f"xrdcp -f root://eosuser.cern.ch//{eos_normalized} {new_root_file_name}")
         else:
             ret = os.system(f"cp {old_root_file_name} {new_root_file_name}")
         if ret != 0:
             raise RuntimeError(f"Failed to copy {old_root_file_name} -> {new_root_file_name}")
         if not os.path.exists(new_root_file_name) or os.path.getsize(new_root_file_name) == 0:
             raise RuntimeError(f"Local copy is missing or empty: {new_root_file_name} (source: {old_root_file_name})")
         rootFile    = ROOT.TFile.Open( new_root_file_name, "READ")

         rootFileNew.cd()
         histograms[folderHR] = {}

         for cutName in cuts :
           for variableName, variable in variables.items():
             folder_name = cutName + "/" + variableName
             histograms[folderHR][folder_name] = {}
             folder_in_file = rootFile.Get(folder_name)

             if not folder_in_file:
               print(f"  WARNING: {folder_name} not found in {new_root_file_name}, skipping")
               continue

             for k in folder_in_file.GetListOfKeys():
               h = k.ReadObj()
               histoName = h.GetName()
               match = re.search("histo_", histoName)
               if not match:
                 continue

               # [FEATURE] apply sample rename: check if histoName starts with any old histo prefix
               # Handles both nominal (histo_DY) and variations (histo_DY_puWeightup)
               renamed_histoName = histoName
               for old_histo, new_histo in histo_rename.items():
                 if histoName == old_histo or histoName.startswith(old_histo + '_'):
                   renamed_histoName = new_histo + histoName[len(old_histo):]
                   break

               cloned = h.Clone()
               cloned.SetName(renamed_histoName)
               histograms[folderHR][folder_name][renamed_histoName] = cloned


       #
       # add the nominals for all the years
       #

       rootFileNew.cd()

       for cutName in cuts :
         for variableName, variable in variables.items():
           folder_name = cutName + "/" + variableName
           rootFileNew.cd()
           rootFileNew.mkdir(folder_name)
           rootFileNew.cd(folder_name)

           for sampleName, sample in samples.items():

             histos_to_be_summed = []

             for folderHR, folder in foldersToMerge.items():
               histoName = "histo_" + sampleName

               if folder_name not in histograms[folderHR]:
                 print(f"  WARNING: {folder_name} missing for {folderHR}, skipping")
                 continue
               if histoName not in histograms[folderHR][folder_name]:
                 print(f"  WARNING: {histoName} missing in {folderHR}/{folder_name}, skipping")
                 continue
               histos_to_be_summed.append(histograms[folderHR][folder_name][histoName])

             if not histos_to_be_summed:
               continue

             summed_histo = histos_to_be_summed[0].Clone()

             for hh in histos_to_be_summed[1:] :
               summed_histo.Add(hh)

             summed_histo.Write()


       #
       # Now let's handle the nuisances
       #

       rootFileNew.cd()

       print (" new_list_of_nuisances = ", new_list_of_nuisances)

       for cutName in cuts :
         for variableName, variable in variables.items():
           folder_name = cutName + "/" + variableName
           rootFileNew.cd()
           rootFileNew.cd(folder_name)

           for nuisanceName, nuisance in new_list_of_nuisances.items():

             if nuisance['type'] == 'lnN' :
               continue

             elif 'name' in nuisance.keys() :

               for sampleName, sample in samples.items():

                 nameTempUp   = 'histo_' + str(sampleName) + '_' + (nuisance['name']) + 'up'
                 nameTempDown = 'histo_' + str(sampleName) + '_' + (nuisance['name']) + 'do'
                 nameTemp     = 'histo_' + str(sampleName)

                 histos_up_to_be_summed = []
                 histos_down_to_be_summed = []

                 histos_up_to_be_summed_weights = []
                 histos_down_to_be_summed_weights = []

                 for folderHR, folder in foldersToMerge.items():

                   if folder_name not in histograms[folderHR]:
                     continue

                   if nameTempUp in histograms[folderHR][folder_name].keys() and nameTempDown in histograms[folderHR][folder_name].keys() :
                     histos_up_to_be_summed.append  ( histograms[folderHR][folder_name][nameTempUp] )
                     histos_down_to_be_summed.append( histograms[folderHR][folder_name][nameTempDown] )
                     histos_up_to_be_summed_weights.append ( 1.0 )
                     histos_down_to_be_summed_weights.append ( 1.0 )

                   else :

                     if nuisanceName in all_nuisances[folderHR].keys() :

                       if all_nuisances[folderHR][nuisanceName]['type'] == 'lnN' and sampleName in all_nuisances[folderHR][nuisanceName]['samples'].keys() :

                         if nameTemp not in histograms[folderHR][folder_name]:
                           continue

                         histos_up_to_be_summed.append  ( histograms[folderHR][folder_name][nameTemp] )
                         histos_down_to_be_summed.append( histograms[folderHR][folder_name][nameTemp] )
                         if "/" not in all_nuisances[folderHR][nuisanceName]['samples'][sampleName]:
                           histos_up_to_be_summed_weights.append ( float(all_nuisances[folderHR][nuisanceName]['samples'][sampleName]) )
                           histos_down_to_be_summed_weights.append ( 1. / float(all_nuisances[folderHR][nuisanceName]['samples'][sampleName]) )
                         else :
                           val_up, val_down = all_nuisances[folderHR][nuisanceName]['samples'][sampleName].split("/")
                           histos_up_to_be_summed_weights.append ( val_up )
                           histos_down_to_be_summed_weights.append ( val_down )

                     else:
                       if nameTemp not in histograms[folderHR][folder_name]:
                         continue
                       histos_up_to_be_summed.append  ( histograms[folderHR][folder_name][nameTemp] )
                       histos_down_to_be_summed.append( histograms[folderHR][folder_name][nameTemp] )
                       histos_up_to_be_summed_weights.append   ( 1.0 )
                       histos_down_to_be_summed_weights.append ( 1.0 )


                 if len(histos_up_to_be_summed) >= 1 and sampleName != "DATA":
                   print ("nuisanceName to be written:", nuisanceName, "for sample: ", sampleName, " cut: ", cutName, " variable: ", variableName)
                   summed_up_histo = histos_up_to_be_summed[0].Clone()
                   summed_up_histo.Scale(histos_up_to_be_summed_weights[0])
                   ihh = 0
                   for hh in histos_up_to_be_summed[1:] :
                     ihh += 1
                     summed_up_histo.Add(hh, histos_up_to_be_summed_weights[ihh])
                   summed_up_histo.SetName (nameTempUp)
                   summed_up_histo.Write()

                   summed_down_histo = histos_down_to_be_summed[0].Clone()
                   summed_down_histo.Scale(histos_down_to_be_summed_weights[0])
                   ihh = 0
                   for hh in histos_down_to_be_summed[1:] :
                     ihh += 1
                     summed_down_histo.Add(hh, histos_down_to_be_summed_weights[ihh])
                   summed_down_histo.SetName (nameTempDown)
                   summed_down_histo.Write()


       return True




if __name__ == '__main__':
    sys.argv = argv

    ROOT.gROOT.SetBatch(True)

    header = """
         ------------------------------------------------------------------------------
         '          __  __                      __   __                               '
         '         |  \\/  | ___ _ __ __ _  ___  \\ \\ / /__  __ _ _ __ ___              '
         '         | |\\/| |/ _ \\ '__/ _` |/ _ \\  \\ V / _ \\/ _` | '__/ __|             '
         '         | |  | |  __/ | | (_| |  __/   | |  __/ (_| | |  \\__ \\             '
         '         |_|  |_|\\___|_|  \\__, |\\___|   |_|\\___|\\__,_|_|  |___/             '
         '                          |___/                                             '
         '                                                                            '
         ------------------------------------------------------------------------------
         """

    print(header)

    parser = defaultParser()

    parser.add_argument("--silentMode",  action='store_true', dest="silentMode",    help="Remove as much as possible the print and the log/err files production")

    opt = parser.parse_args()
    print ("opt.pycfg            = ", opt.pycfg)
    loadDefaultOptions(parser, opt.pycfg)
    opt = parser.parse_args()

    print ("opt.tag              = ", opt.tag)
    print ("opt.variablesFile    = ", opt.variablesFile)
    print ("opt.cutsFile         = ", opt.cutsFile)
    print ("opt.samplesFile      = ", opt.samplesFile)
    print ("opt.structureFile    = ", opt.structureFile)
    print ("opt.nuisancesFile    = ", opt.nuisancesFile)
    print ("opt.foldersToMerge   = ", opt.foldersToMerge)

    foldersToMergeNuisancesFiles = {
        k : os.path.join(p["folder"], 'nuisances.py') for k, p in opt.foldersToMerge.items()
        }
    print ("foldersToMergeNuisancesFiles   = ", foldersToMergeNuisancesFiles)

    factory = MergerFactory()

    samples = {}
    if os.path.exists(opt.samplesFile) :
      handle = open(opt.samplesFile,'r')
      exec(handle.read())
      handle.close()
      samples = {k: v for k, v in samples.items() if not (k.startswith('__') and k.endswith('__'))}

    variables = {}
    if os.path.exists(opt.variablesFile) :
      handle = open(opt.variablesFile,'r')
      exec(handle.read())
      handle.close()
      variables = {k: v for k, v in variables.items() if not (k.startswith('__') and k.endswith('__'))}

    print ("variables = ", variables)

    cuts = {}
    if os.path.exists(opt.cutsFile) :
      handle = open(opt.cutsFile,'r')
      exec(handle.read())
      handle.close()
      cuts = {k: v for k, v in cuts.items() if not (k.startswith('__') and k.endswith('__'))}

    print ("cuts = ", cuts)

    nuisances = {}
    if os.path.exists(opt.nuisancesFile) :
      handle = open(opt.nuisancesFile,'r')
      exec(handle.read())
      handle.close()
      nuisances = {k: v for k, v in nuisances.items() if not (k.startswith('__') and k.endswith('__'))}

    print ("nuisances = ", nuisances)

    outputDir = getattr(opt, 'outputDir', None)
    plotDir   = getattr(opt, 'plotDir',   None)

    factory.merge(
       opt.tag,
       variables,
       cuts,
       samples,
       nuisances,
       opt.foldersToMerge,
       foldersToMergeNuisancesFiles,
       outputDir=outputDir,
       plotDir=plotDir,
       )

    print ("\n\n")
    print (" I'm done ... \n\n")