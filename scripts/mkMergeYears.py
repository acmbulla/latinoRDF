#!/usr/bin/env python3

import sys
argv = sys.argv
sys.argv = argv[:1]

import argparse

import os
import re
import tempfile

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




#
# output file naming: never overwrite an existing merge
#

def versionedOutputPath(directory, base="root_file_joined.root"):
    """
    First merge writes <base>. Every later one is suffixed _v1, _v2, ... so an
    existing file is never clobbered: the suffix goes on the *new* file and the
    ones already there are left untouched.

    NB mkPlot defaults --inputFile to root_file_joined.root and mkDatacards
    hardcodes it (mkDatacards.py:909), so from the second merge on they still
    read the first file. Pass the versioned path explicitly.
    """
    stem, ext = os.path.splitext(base)
    candidate = os.path.join(directory, base) if directory else base
    if not os.path.exists(candidate):
        return candidate

    i = 1
    while True:
        name = "%s_v%d%s" % (stem, i, ext)
        candidate = os.path.join(directory, name) if directory else name
        if not os.path.exists(candidate):
            return candidate
        i += 1


def readFolder(rootFile, folder_name, histo_rename):
    """
    Read one cut/variable folder of one year into {histoName: TH1}.

    Only the latest version of each histogram is taken (highest key cycle
    number). haddfast.py writes a new cycle every time it adds into an existing
    histogram, so a file merged with --haddBatch carries ;1 -- the first chunk
    alone -- next to the highest cycle, which holds the running sum.
    GetListOfKeys() returns them in *descending* cycle order, so keying a dict
    by name and letting the last iteration win would keep the partial histogram
    and silently under-count every sample. Same dedup as haddfast.merge().

    The histograms are cloned, so the caller owns them and has to Delete() them;
    the originals belong to the directory and go away with the file.
    """
    out = {}

    folder_in_file = rootFile.Get(folder_name)
    if not folder_in_file:
        print(f"  WARNING: {folder_name} not found in {rootFile.GetName()}, skipping")
        return out

    inkeys = {}
    for k in folder_in_file.GetListOfKeys():
        keyName = k.GetName()
        if keyName in inkeys and inkeys[keyName].GetCycle() > k.GetCycle():
            continue
        inkeys[keyName] = k

    for histoName, k in inkeys.items():
        if not histoName.startswith("histo_"):
            continue

        # [FEATURE] apply sample rename: check if histoName starts with any old
        # histo prefix. Handles both nominal (histo_DY) and variations
        # (histo_DY_puWeightup)
        renamed_histoName = histoName
        for old_histo, new_histo in histo_rename.items():
            if histoName == old_histo or histoName.startswith(old_histo + '_'):
                renamed_histoName = new_histo + histoName[len(old_histo):]
                break

        cloned = k.ReadObj().Clone()
        cloned.SetName(renamed_histoName)
        cloned.SetDirectory(0)
        out[renamed_histoName] = cloned

    return out


#
# batch splitting: one job per cut
#
# The merge touches one output folder, cut/variable, at a time and a folder of
# one cut never depends on another cut, so the work splits cleanly by cut. Each
# job writes its own partial file and the partials hold disjoint folders, so
# hadding them together is a union and cannot double count.
#

PARTIAL_PREFIX = "root_file_joined___"


def partialName(cutName):
    return PARTIAL_PREFIX + cutName + ".root"


def expandedCutNames(cuts):
    """The cut names the merge will produce, categories unrolled."""
    names = []
    for cutName, cut in cuts.items():
        names.append(cutName)
        if isinstance(cut, dict):
            for category_name in cut.get('categories', {}):
                names.append(cutName + "_" + category_name)
    return names


def checkPartials(outputDir, cutNames):
    """(missing, present) partial files, by cut name."""
    missing, present = [], []
    for c in cutNames:
        path = os.path.join(outputDir, partialName(c))
        if os.path.exists(path) and os.path.getsize(path) > 0:
            present.append(c)
        else:
            missing.append(c)
    return missing, present


def submitBatch(pycfg, outputDir, cutNames, batchDir="jobsMerge", jobFlavour="workday"):
    """Write one condor job per cut and submit them as a single cluster."""
    import subprocess

    batchDir = os.path.abspath(batchDir)
    os.makedirs(batchDir, exist_ok=True)
    outputDir = os.path.abspath(outputDir)
    os.makedirs(outputDir, exist_ok=True)

    script_path = os.path.abspath(__file__)
    pycfg_abs = os.path.abspath(pycfg)
    workdir = os.getcwd()

    sh_path = os.path.join(batchDir, "merge_one_cut.sh")
    with open(sh_path, "w") as f:
        f.write(
            '#!/bin/bash\n'
            'set -e\n'
            'CUT=$1\n'
            'echo "merge job for cut $CUT started at $(date) on $(hostname)"\n'
            'cd ' + workdir + '\n'
            'python ' + script_path + ' \\\n'
            '    --pycfg ' + pycfg_abs + ' \\\n'
            '    --cut "$CUT" \\\n'
            '    --output ' + PARTIAL_PREFIX + '"$CUT".root \\\n'
            '    --silentMode\n'
            'echo "done at $(date)"\n'
            'ls -lh ' + outputDir + '/' + PARTIAL_PREFIX + '"$CUT".root\n'
        )
    os.system("chmod +x " + sh_path)

    sub_path = os.path.join(batchDir, "submit_merge.sub")
    with open(sub_path, "w") as f:
        f.write(
            '# one merge job per cut\n'
            'executable              = ' + sh_path + '\n'
            'arguments               = $(cutname)\n'
            'output                  = ' + batchDir + '/merge_$(cutname).out\n'
            'error                   = ' + batchDir + '/merge_$(cutname).err\n'
            'log                     = ' + batchDir + '/merge_$(cutname).log\n'
            'getenv                  = True\n'
            '+JobFlavour             = "' + jobFlavour + '"\n'
            '\n'
            'queue cutname in (\n'
        )
        for c in cutNames:
            f.write("  " + c + "\n")
        f.write(")\n")

    print("  -> job script  : %s" % sh_path)
    print("  -> submit file : %s" % sub_path)
    print("  -> %d jobs, one per cut" % len(cutNames))
    print("  -> partials    : %s/%s<cut>.root" % (outputDir, PARTIAL_PREFIX))
    print("Submitting...")
    subprocess.run(["condor_submit", sub_path], cwd=batchDir)


def haddPartials(outputDir, cutNames, haddfastPath=None):
    """Merge the per-cut partials into the versioned root_file_joined.root."""
    import subprocess

    missing, present = checkPartials(outputDir, cutNames)
    if missing:
        print("ERROR: %d partial file(s) missing, refusing to hadd:" % len(missing))
        for c in missing:
            print("  %s" % os.path.join(outputDir, partialName(c)))
        print("Run --checkBatch for the full picture.")
        return False

    if haddfastPath is None:
        haddfastPath = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "haddfast.py")

    target = versionedOutputPath(outputDir)
    sources = [os.path.join(outputDir, partialName(c)) for c in present]

    # -C purges the duplicate key cycles haddfast leaves behind, so the merged
    # file carries one cycle per histogram like a local merge does.
    cmd = [sys.executable, haddfastPath, "-f", "-C", target] + sources
    print("  -> %d partials -> %s" % (len(sources), target))
    rc = subprocess.run(cmd).returncode
    if rc != 0:
        print("ERROR: hadd failed with code %d" % rc)
        return False
    print("  -> done: %s" % target)
    return True


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
        silentMode=False,
        cutSubset=None,
        outputName=None,
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

       # The merge loop runs over list_cuts, not over the raw cuts dict, so a cut
       # carrying 'categories' contributes its unrolled sub-regions too. cutSubset
       # restricts the merge to part of that list: this is how a batch job gets its
       # own slice, and because the slices are disjoint folders in the output the
       # partial files can simply be hadded together afterwards.
       if cutSubset:
         unknown = [c for c in cutSubset if c not in list_cuts]
         if unknown:
           raise ValueError("unknown cut(s) %s; available: %s"
                            % (", ".join(unknown), ", ".join(list_cuts)))
         list_cuts = {c: list_cuts[c] for c in cutSubset}
         print("Merging only: %s" % ", ".join(list_cuts))


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
       if outputDir:
           os.makedirs(outputDir, exist_ok=True)
       if outputName:
           out_root_file_name = (os.path.join(outputDir, outputName)
                                 if outputDir else outputName)
       else:
           out_root_file_name = versionedOutputPath(outputDir or "")
       print(f"Output root file: {out_root_file_name}")
       if plotDir:
           os.makedirs(plotDir, exist_ok=True)
           print(f"Plot directory:   {plotDir}")
       rootFileNew = ROOT.TFile.Open( out_root_file_name, "RECREATE")

       #
       # Open every year once. The files stay open for the whole merge: the old
       # code reassigned `rootFile` inside the year loop without closing it, so
       # the handles and the local copies leaked until the process exited.
       #
       year_files  = {}
       year_rename = {}
       local_copies = []

       for folderHR, folder in foldersToMerge.items():

         # [FEATURE] sampleRename: optional dict {oldName: newName} per year
         # Used to rename sample histograms when sample names differ between years
         # Example: {'DY_2024': 'DY', 'top_v2': 'top'}
         sample_rename = folder.get('sampleRename', {})
         # build histo-level rename map: histo_oldName -> histo_newName
         # also handles nuisance variations: histo_oldName_<nuis>up -> histo_newName_<nuis>up
         year_rename[folderHR] = {f"histo_{old}": f"histo_{new}" for old, new in sample_rename.items()}

         old_root_file_name = os.path.join(folder["folder"], "rootFile", folder["tag"], "root_file_joined.root")
         tmp_dir = tempfile.gettempdir()
         local_year_file = os.path.join(tmp_dir, "year_" + folderHR + "_histos_" + folder["tag"] + ".root")
         # use xrdcp for EOS paths, plain cp otherwise
         if "/eos/" in old_root_file_name:
             eos_normalized = re.sub(r'^/eos/home-([a-z])/', lambda m: f"/eos/user/{m.group(1)}/", old_root_file_name)
             ret = os.system(f"xrdcp -f root://eosuser.cern.ch//{eos_normalized} {local_year_file}")
         else:
             ret = os.system(f"cp {old_root_file_name} {local_year_file}")
         if ret != 0:
             raise RuntimeError(f"Failed to copy {old_root_file_name} -> {local_year_file}")
         if not os.path.exists(local_year_file) or os.path.getsize(local_year_file) == 0:
             raise RuntimeError(f"Local copy is missing or empty: {local_year_file} (source: {old_root_file_name})")

         year_files[folderHR]  = ROOT.TFile.Open( local_year_file, "READ")
         local_copies.append(local_year_file)
         print(f"  opened {folderHR}: {local_year_file}")

       #
       # One output folder at a time: read cut/variable from every year, sum the
       # nominals, sum the nuisances, write, free. Only the folder being written
       # is held, instead of every year's whole file at once.
       #
       # What this buys, measured on 7 years x 10 cuts x 8 variables x 25 samples
       # x 2 shape nuisances (9680 output histograms), against the load-everything
       # version it replaces: 352 s -> 157 s, with bit-identical output.
       #
       # The speed comes from the clones no longer piling up in the output file's
       # in-memory object list. The old code called rootFileNew.cd() and then
       # Clone()d every histogram of every year, so all of them were appended to
       # one directory, which makes the cost quadratic: cloning 40k histograms
       # into a directory takes 4.7 s with distinct names and 70 s when the names
       # repeat, as they do here (~500 names reused in every folder of every
       # year), against 0.17 s detached with SetDirectory(0).
       #
       # Peak RSS is NOT the win: it measured 392 MB before and 401 MB after, so
       # at this scale the peak is ROOT plus the open input files, not the
       # histograms. What streaming removes is the growth with the number of
       # variables and shape nuisances, which the old shape could not survive.
       #
       # The remaining floor is per-histogram object lifetime: creating and
       # freeing a histogram read from a file costs ~0.17 ms whichever way it is
       # freed (explicit Delete() and dropping the reference measured the same),
       # and the merge must do that once per histogram per year. That is why this
       # splits by cut for the batch rather than being optimised further here.
       #

       for cutName in list_cuts :
         for variableName, variable in variables.items():
           folder_name = cutName + "/" + variableName

           per_year = {}
           for folderHR in foldersToMerge :
             per_year[folderHR] = readFolder(year_files[folderHR], folder_name,
                                             year_rename[folderHR])

           rootFileNew.cd()
           rootFileNew.mkdir(folder_name)
           rootFileNew.cd(folder_name)

           #
           # add the nominals for all the years
           #

           for sampleName, sample in samples.items():

             histos_to_be_summed = []

             for folderHR, folder in foldersToMerge.items():
               histoName = "histo_" + sampleName

               if not per_year[folderHR]:
                 print(f"  WARNING: {folder_name} missing for {folderHR}, skipping")
                 continue
               if histoName not in per_year[folderHR]:
                 print(f"  WARNING: {histoName} missing in {folderHR}/{folder_name}, skipping")
                 continue
               histos_to_be_summed.append(per_year[folderHR][histoName])

             if not histos_to_be_summed:
               continue

             summed_histo = histos_to_be_summed[0].Clone()

             for hh in histos_to_be_summed[1:] :
               summed_histo.Add(hh)

             summed_histo.Write()
             summed_histo.Delete()

           #
           # Now let's handle the nuisances
           #

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

                   if not per_year[folderHR]:
                     continue

                   if nameTempUp in per_year[folderHR].keys() and nameTempDown in per_year[folderHR].keys() :
                     histos_up_to_be_summed.append  ( per_year[folderHR][nameTempUp] )
                     histos_down_to_be_summed.append( per_year[folderHR][nameTempDown] )
                     histos_up_to_be_summed_weights.append ( 1.0 )
                     histos_down_to_be_summed_weights.append ( 1.0 )

                   else :

                     if nuisanceName in all_nuisances[folderHR].keys() :

                       if all_nuisances[folderHR][nuisanceName]['type'] == 'lnN' and sampleName in all_nuisances[folderHR][nuisanceName]['samples'].keys() :

                         if nameTemp not in per_year[folderHR]:
                           continue

                         histos_up_to_be_summed.append  ( per_year[folderHR][nameTemp] )
                         histos_down_to_be_summed.append( per_year[folderHR][nameTemp] )
                         if "/" not in all_nuisances[folderHR][nuisanceName]['samples'][sampleName]:
                           histos_up_to_be_summed_weights.append ( float(all_nuisances[folderHR][nuisanceName]['samples'][sampleName]) )
                           histos_down_to_be_summed_weights.append ( 1. / float(all_nuisances[folderHR][nuisanceName]['samples'][sampleName]) )
                         else :
                           val_up, val_down = all_nuisances[folderHR][nuisanceName]['samples'][sampleName].split("/")
                           histos_up_to_be_summed_weights.append ( val_up )
                           histos_down_to_be_summed_weights.append ( val_down )

                     else:
                       if nameTemp not in per_year[folderHR]:
                         continue
                       histos_up_to_be_summed.append  ( per_year[folderHR][nameTemp] )
                       histos_down_to_be_summed.append( per_year[folderHR][nameTemp] )
                       histos_up_to_be_summed_weights.append   ( 1.0 )
                       histos_down_to_be_summed_weights.append ( 1.0 )


                 if len(histos_up_to_be_summed) >= 1 and sampleName != "DATA":
                   if not silentMode:
                     print ("nuisanceName to be written:", nuisanceName, "for sample: ", sampleName, " cut: ", cutName, " variable: ", variableName)
                   summed_up_histo = histos_up_to_be_summed[0].Clone()
                   summed_up_histo.Scale(histos_up_to_be_summed_weights[0])
                   ihh = 0
                   for hh in histos_up_to_be_summed[1:] :
                     ihh += 1
                     summed_up_histo.Add(hh, histos_up_to_be_summed_weights[ihh])
                   summed_up_histo.SetName (nameTempUp)
                   summed_up_histo.Write()
                   summed_up_histo.Delete()

                   summed_down_histo = histos_down_to_be_summed[0].Clone()
                   summed_down_histo.Scale(histos_down_to_be_summed_weights[0])
                   ihh = 0
                   for hh in histos_down_to_be_summed[1:] :
                     ihh += 1
                     summed_down_histo.Add(hh, histos_down_to_be_summed_weights[ihh])
                   summed_down_histo.SetName (nameTempDown)
                   summed_down_histo.Write()
                   summed_down_histo.Delete()

           #
           # the folder is on disk: drop everything it needed
           #
           for _byname in per_year.values():
             for _h in _byname.values():
               _h.Delete()
           per_year.clear()

       rootFileNew.Close()

       for folderHR, tf in year_files.items():
         tf.Close()
       for path in local_copies:
         try:
           os.unlink(path)
         except OSError:
           pass


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
    parser.add_argument("--cut",         dest="cut",         default=None,
                        help="Merge only these cuts, comma separated. One batch job gets one cut.")
    parser.add_argument("--output",      dest="output",      default=None,
                        help="Explicit output file name, relative to outputDir. Bypasses the _vN versioning.")
    parser.add_argument("--submitBatch", action='store_true', dest="submitBatch",
                        help="Submit one condor job per cut instead of merging locally")
    parser.add_argument("--checkBatch",  action='store_true', dest="checkBatch",
                        help="Report which per-cut partial files are missing")
    parser.add_argument("--hadd",        action='store_true', dest="hadd",
                        help="Merge the per-cut partials into root_file_joined.root")
    parser.add_argument("--batchDir",    dest="batchDir",    default="jobsMerge",
                        help="Where the condor scripts and logs go")

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

    cutNames = expandedCutNames(cuts)

    if opt.submitBatch:
        submitBatch(opt.pycfg, outputDir or ".", cutNames, batchDir=opt.batchDir)

    elif opt.checkBatch:
        missing, present = checkPartials(outputDir or ".", cutNames)
        print("  -> %d/%d partials present" % (len(present), len(cutNames)))
        for c in missing:
            print("  MISSING: %s" % os.path.join(outputDir or ".", partialName(c)))
        if not missing:
            print("  -> all done, run --hadd")

    elif opt.hadd:
        haddPartials(outputDir or ".", cutNames)

    else:
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
           silentMode=opt.silentMode,
           cutSubset=[c.strip() for c in opt.cut.split(",")] if opt.cut else None,
           outputName=opt.output,
           )

    print ("\n\n")
    print (" I'm done ... \n\n")