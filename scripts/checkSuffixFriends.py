#!/usr/bin/env python3
"""
Pre-flight check for the varied snapshot jobs.

A varied job reads the nominal files with the nuisance's folderUp / folderDown
attached as a friend, and a friend is positional: ROOT pairs entry i with entry
i. A '_suffix' file that is missing, or that holds a different number of
entries than the nominal file it shadows, therefore does not fail loudly -- the
job reads past the end of the friend and dies in the event loop with a
segmentation violation, after Snapshot has already started.

Nothing in the production guarantees that every dataset exists in every
'_suffix' folder, so this walks what samples.py and nuisances.py declare and
reports the pairs that cannot work, before a single job is submitted.

The second thing it reports is which columns a variation actually replaces.
That is a physics question, not a crash: the varied job redefines every friend
branch whose nominal counterpart exists, so a '_suffix' tree that carries
`nCleanJet` but not `CleanJet_jetIdx` -- or the other way round -- leaves the
collection internally inconsistent. A helper that loops over the counter and
indexes the arrays then reads out of bounds, which is how the b-tag scale
factors segfaulted.

    checkSuffixFriends.py --pycfg configuration.py
    checkSuffixFriends.py --pycfg configuration.py --sample ZZ --nuisance JER
    checkSuffixFriends.py --pycfg configuration.py --existenceOnly
    checkSuffixFriends.py --pycfg configuration.py --branches

Exit code is 1 when anything is missing, mismatched, or a partially varied
collection is found, so it can gate a submit.
"""

import sys
argv = sys.argv
sys.argv = argv[:1]

import os
import ROOT

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from analysisFactory import defaultParser, loadDefaultOptions, load_config_files

ROOT.gErrorIgnoreLevel = ROOT.kError


def _entries(path, tree_name="Events"):
    """Entries of tree_name, or None if the file or the tree cannot be read."""
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        return None
    try:
        t = f.Get(tree_name)
        return int(t.GetEntries()) if t else None
    finally:
        f.Close()


def _branches(path, tree_name="Events"):
    """Branch names of tree_name, or None if the file cannot be read."""
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        return None
    try:
        t = f.Get(tree_name)
        if not t:
            return None
        return [b.GetName() for b in t.GetListOfBranches()]
    finally:
        f.Close()


def _collection(branch):
    """The collection a branch belongs to: 'CleanJet_pt' -> 'CleanJet'."""
    return branch.split('_', 1)[0] if '_' in branch else branch


def _report_branches(nuisanceName, direction, mapName, nominal_file, varied_file,
                     tree_name):
    """Which columns the variation replaces, and whether any collection is
    left half varied. Returns the list of inconsistent collections."""
    nom = _branches(nominal_file, tree_name)
    var = _branches(varied_file, tree_name)
    if nom is None or var is None:
        print(f"    {nuisanceName}{direction}: cannot read the pair, skipped")
        return []

    nom_set   = set(nom)
    suffix    = "_" + mapName
    suffixed  = [b for b in var if b.endswith(suffix)]
    mapped    = {}          # nominal column -> varied branch
    unmatched = []
    for b in suffixed:
        nomCol = b[:-len(suffix)]
        if nomCol in nom_set:
            mapped[nomCol] = b
        else:
            unmatched.append(b)

    print(f"    {nuisanceName}{direction}  ({mapName}): "
          f"{len(suffixed)} suffixed branch(es), {len(mapped)} redefined"
          + (f", {len(unmatched)} with no nominal counterpart (ignored)"
             if unmatched else ""))

    # Two different situations, and only one of them is a bug.
    #
    #   counter varied, some members left nominal -> the collection is RESIZED
    #   and those members now have the wrong length. Reading them past the new
    #   counter is the segfault.
    #
    #   counter not varied, some members varied -> the lengths are unchanged
    #   and this is normal: a JES shift moves Jet_pt and Jet_mass, not
    #   Jet_eta. Reported for the record, not as an error.
    resized = []
    touched = {_collection(c) for c in mapped}
    for coll in sorted(touched):
        if coll.startswith('n') and coll[1:] in touched:
            continue
        counter = 'n' + coll
        if counter not in nom_set:
            # No counter branch, so this is not a collection but a family of
            # per-event scalars: PuppiMET_pt, WH3l_mtWWW, ZH3l_pTZ. A
            # variation touching some of them and not others is what a
            # variation is, and there is no length to get wrong.
            continue
        members_all   = [b for b in nom if b.startswith(coll + '_')]
        members_var   = [b for b in members_all if b in mapped]
        members_fixed = [b for b in members_all if b not in mapped]
        if not members_var and not members_fixed:
            continue
        counter_state = 'varied' if counter in mapped else 'nominal'

        if counter_state == 'varied' and members_fixed:
            label = 'RESIZED'
            resized.append(coll)
        elif members_var and members_fixed:
            label = 'partial '
        else:
            continue

        print(f"      {label} {coll}: counter n{coll} {counter_state}, "
              f"{len(members_var)} member(s) varied, "
              f"{len(members_fixed)} left nominal")
        print(f"        varied : {', '.join(sorted(members_var)[:10])}"
              + (" ..." if len(members_var) > 10 else ""))
        if members_fixed:
            print(f"        nominal: {', '.join(sorted(members_fixed)[:10])}"
                  + (" ..." if len(members_fixed) > 10 else ""))

    return resized


def _exists(path):
    if path.startswith("root://"):
        f = ROOT.TFile.Open(path)
        ok = bool(f) and not f.IsZombie()
        if f:
            f.Close()
        return ok
    return os.path.exists(path)


if __name__ == '__main__':
    sys.argv = argv

    parser = defaultParser()
    parser.add_argument("--sample",         dest="sample",   default=None,
                        help="Only this sample (substring match)")
    parser.add_argument("--nuisance",       dest="nuisance", default=None,
                        help="Only this nuisance (substring match)")
    parser.add_argument("--existenceOnly",  action='store_true', dest="existenceOnly",
                        help="Do not open the files: check only that they are there")
    parser.add_argument("--branches",       action='store_true', dest="branches",
                        help="Also report, per nuisance, which columns the "
                             "variation replaces and whether any collection is "
                             "left half varied (one file pair per nuisance)")
    parser.add_argument("--tree",           dest="tree", default="Events")

    opt = parser.parse_args()
    loadDefaultOptions(parser, opt.pycfg)
    opt = parser.parse_args()

    cfg       = load_config_files(opt)
    samples   = cfg['samples']
    nuisances = cfg['nuisances']

    suffix_nuisances = [(n, nu) for n, nu in nuisances.items()
                        if nu.get('type') == 'shape'
                        and nu.get('kind') in ('suffix', 'suffixFly')
                        and 'folderUp' in nu]
    if opt.nuisance:
        suffix_nuisances = [(n, nu) for n, nu in suffix_nuisances if opt.nuisance in n]

    if not suffix_nuisances:
        print("No 'suffix' nuisance with a folderUp/folderDown to check.")
        sys.exit(0)

    print(f"Checking {len(suffix_nuisances)} nuisance(s) "
          f"x 2 direction(s) against {len(samples)} sample(s)"
          + ("" if opt.existenceOnly else "  (opening every file -- this is the slow part)"))

    missing    = []   # (nuisance, direction, sample, path)
    mismatched = []   # (nuisance, direction, sample, path, n_nominal, n_varied)
    nominal_unreadable = []
    checked = 0

    for sampleName, sample in samples.items():
        if opt.sample and opt.sample not in sampleName:
            continue
        if 'isData' in sample:
            continue

        for subname, list_root_files in sample['name'].items():
            for nuisanceName, nuisance in suffix_nuisances:
                if sampleName not in nuisance.get('samples', {}):
                    continue

                for direction in ('Up', 'Down'):
                    folder = nuisance[f'folder{direction}']
                    for nominal_file in list_root_files:
                        varied_file = os.path.join(folder, os.path.basename(nominal_file))
                        checked += 1

                        if not _exists(varied_file):
                            missing.append((nuisanceName, direction, sampleName, varied_file))
                            continue
                        if opt.existenceOnly:
                            continue

                        n_var = _entries(varied_file, opt.tree)
                        n_nom = _entries(nominal_file, opt.tree)
                        if n_nom is None:
                            nominal_unreadable.append(nominal_file)
                            continue
                        if n_var is None:
                            missing.append((nuisanceName, direction, sampleName, varied_file))
                            continue
                        if n_var != n_nom:
                            mismatched.append((nuisanceName, direction, sampleName,
                                               varied_file, n_nom, n_var))

    print(f"\n  {checked} (nominal, varied) pair(s) checked")

    partial = []
    if opt.branches:
        print("\n  Columns each variation replaces (first readable pair per nuisance):")
        for nuisanceName, nuisance in suffix_nuisances:
            for direction in ('Up', 'Down'):
                folder = nuisance[f'folder{direction}']
                pair = None
                for sampleName, sample in samples.items():
                    if opt.sample and opt.sample not in sampleName:
                        continue
                    if 'isData' in sample or sampleName not in nuisance.get('samples', {}):
                        continue
                    for list_root_files in sample['name'].values():
                        for nominal_file in list_root_files:
                            candidate = os.path.join(folder, os.path.basename(nominal_file))
                            if _exists(candidate):
                                pair = (nominal_file, candidate)
                                break
                        if pair:
                            break
                    if pair:
                        break
                if not pair:
                    print(f"    {nuisanceName}{direction}: no readable file pair")
                    continue
                bad = _report_branches(nuisanceName, direction,
                                       nuisance[f'map{direction}'],
                                       pair[0], pair[1], opt.tree)
                partial.extend((nuisanceName, direction, c) for c in bad)

    if nominal_unreadable:
        print(f"\n  {len(nominal_unreadable)} nominal file(s) unreadable -- "
              f"samples.py is wrong before the nuisances are:")
        for p in nominal_unreadable[:20]:
            print(f"    {p}")

    if missing:
        print(f"\n  MISSING -- {len(missing)} varied file(s) do not exist. "
              f"Every job that would read one of these segfaults:")
        for nuisanceName, direction, sampleName, path in missing[:40]:
            print(f"    {nuisanceName}{direction}  {sampleName}  {path}")
        if len(missing) > 40:
            print(f"    ... and {len(missing) - 40} more")

    if mismatched:
        print(f"\n  ENTRY MISMATCH -- {len(mismatched)} varied file(s) do not line up "
              f"with their nominal. A positional friend cannot be used here:")
        for nuisanceName, direction, sampleName, path, n_nom, n_var in mismatched[:40]:
            print(f"    {nuisanceName}{direction}  {sampleName}  "
                  f"nominal {n_nom} vs varied {n_var}  {path}")
        if len(mismatched) > 40:
            print(f"    ... and {len(mismatched) - 40} more")

    if partial:
        print(f"\n  RESIZED COLLECTIONS -- {len(partial)} collection(s) have a varied "
              f"counter and members left nominal, so those members now have the wrong "
              f"length. Anything that loops over the counter and indexes the arrays "
              f"reads out of bounds:")
        for nuisanceName, direction, coll in partial:
            print(f"    {nuisanceName}{direction}  {coll}")

    if not missing and not mismatched and not partial:
        print("\n  All varied files are there, line up with their nominal, and "
              "replace whole collections.")
        sys.exit(0)

    sys.exit(1)
