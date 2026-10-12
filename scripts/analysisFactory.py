#!/usr/bin/env python3

import sys
import os
import logging
import subprocess

import ROOT


# -----------------------------------------------------------------------
#  Shared parser helpers
# -----------------------------------------------------------------------

# directory where analysisFactory.py (and haddfast.py, haddnano.py) live
_SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))


def _base_ns():
    """Base namespace for exec — injects builtins + common modules used in config files."""
    import builtins
    import os as _os
    import sys as _sys
    import re as _re
    import glob as _glob
    return {
        '__builtins__': builtins,
        'os'          : _os,
        'sys'         : _sys,
        're'          : _re,
        'glob'        : _glob,
    }


def load_config_files(opt):
    """
    Load all standard config files (samples, variables, cuts, nuisances, aliases)
    from paths in opt. Returns a dict with all loaded objects.
    Pre-initializes dict variables before exec so files that do e.g.
    cuts['SR'] = '...' don't crash.
    """
    # samples first — needed by aliases.py
    samples = read_exec(opt.samplesFile, 'samples')

    variables = read_exec(opt.variablesFile, 'variables')

    # cuts.py typically does: cuts['SR'] = '...' and supercut = '...'
    cuts_ns = {**_base_ns(), 'cuts': {}, 'supercut': '', 'samples': samples}
    if os.path.exists(opt.cutsFile):
        with open(opt.cutsFile, 'r') as fh:
            exec(fh.read(), cuts_ns)
    supercut = cuts_ns.get('supercut', '')
    cuts     = {k: v for k, v in cuts_ns.get('cuts', {}).items()
                if not (k.startswith('__') and k.endswith('__'))}

    # nuisances.py typically does: nuisances['JES'] = {...}
    nuisances = read_exec(opt.nuisancesFile, 'nuisances',
                          extra_ns={'samples': samples, 'cuts': cuts})

    # aliases.py uses samples (allMC, data, etc.) and sometimes cuts
    aliases = read_exec(opt.aliasesFile, 'aliases',
                        extra_ns={'samples': samples, 'cuts': cuts,
                                  'nuisances': nuisances})

    return {
        'samples'  : samples,
        'variables': variables,
        'cuts'     : cuts,
        'supercut' : supercut,
        'nuisances': nuisances,
        'aliases'  : aliases,
    }


def read_exec(path, varname, extra_ns=None):
    """Execute a python config file and return the named variable as a clean dict.
    extra_ns: optional dict injected into the exec namespace (e.g. {'samples': samples})
    """
    d = {**_base_ns(), varname: {}}
    if extra_ns:
        d.update(extra_ns)
    if os.path.exists(path):
        with open(path, 'r') as handle:
            exec(handle.read(), d)
    result = d.get(varname, {})
    if isinstance(result, dict):
        return {k: v for k, v in result.items()
                if not (k.startswith('__') and k.endswith('__'))}
    return result


def defaultParser():
    import argparse
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--pycfg", dest="pycfg",
                        help="input configuration file", default="configuration.py")
    return parser


def loadDefaultOptions(parser, pycfg=None, quiet=False):
    if pycfg is not None and os.path.exists(pycfg):
        handle = open(pycfg, 'r')
        local_variables = {}
        exec(handle.read(), local_variables)
        handle.close()
        local_variables = {k: v for k, v in local_variables.items()
                           if not (k.startswith('__') and k.endswith('__'))}
        for opt_name, opt_value in local_variables.items():
            parser.add_argument('--' + opt_name, default=opt_value)


def read_pyfile(path, dict_name):
    """Read a python file and return the named dict (cleaned of __dunder__ keys)."""
    d = {}
    if os.path.exists(path):
        handle = open(path, 'r')
        exec(handle.read(), {dict_name: d})   # inject dict so assignments land in it
        handle.close()
        # exec puts __builtins__ etc in the namespace — clean it
        d = {k: v for k, v in d.items() if not (k.startswith('__') and k.endswith('__'))}
    return d


# -----------------------------------------------------------------------
#  AnalysisFactory  —  base class
# -----------------------------------------------------------------------

class AnalysisFactory:
    """
    Base class shared by ShapeFactory (mkShapes) and SnapshotFactory (mkSnapshot).

    Contains everything up to and including current_node + cuts:
      - aliases / functions ordered pipeline
      - weights
      - supercut
      - cuts / categories
      - FilesPerJob chunking
      - submitBatch (proxy check, precompiled vs compile-on-node, single .sub file)
      - checkBatch + --resubmit
      - parallelCompile

    Subclasses must implement:
      - job_prefix()             -> str   e.g. "my_run_analysis" or "my_snapshot"
      - output_file_name(sampleName, subname, i) -> str
      - makeScripts()            -> calls self._make_scripts_loop(...)
      - generate_makefile(...)
      - _build_cpp_output(...)   -> str   the C++ fragment after current_node is ready
      - _build_cpp_header()      -> str   extra #includes before main() (optional)
      - _check_output_exists(sampleName, subname, i) -> bool
      - _build_bash_precompiled(...) -> str
      - _build_bash_compile_on_node(...) -> str
      - _extra_transfer_files_for_submit(precompiled, ...) -> str
    """

    _logger = logging.getLogger('AnalysisFactory')

    # _____________________________________________________________________________
    def __init__(self):
        self._variables  = {}
        self._cuts       = {}
        self._supercut   = ''
        self._samples    = {}
        self._aliases    = {}
        self._nuisances  = {}
        self._lumi       = 1
        self._outputDir  = './test/'
        self._script_batch_location = './scripts_batch/'
        self._scripts_run_folder    = 'scripts_run'
        self._silentMode = False
        self._extra_transfer_files  = []

    def __del__(self):
        pass

    # _____________________________________________________________________________
    def setValues(self, outputDir, variables, cuts, supercut, samples, nuisances,
                  aliases, lumi,
                  scripts_run_folder=None, script_batch_location=None,
                  extra_transfer_files=None):
        self._variables  = variables
        self._samples    = samples
        self._cuts       = cuts
        self._supercut   = supercut
        self._nuisances  = nuisances
        self._aliases    = aliases
        self._outputDir  = os.path.expandvars(outputDir)
        self._lumi       = lumi

        if scripts_run_folder is not None:
            self._scripts_run_folder = os.path.expandvars(scripts_run_folder)
        if script_batch_location is not None:
            self._script_batch_location = os.path.expandvars(script_batch_location)
        if extra_transfer_files is not None:
            self._extra_transfer_files = [os.path.expandvars(f) for f in extra_transfer_files]

    def setConditions(self, silentMode):
        self._silentMode = silentMode

    # _____________________________________________________________________________
    @staticmethod
    def _normalize_sample_name(name_field):
        if isinstance(name_field, dict):
            return name_field
        result = {}
        for item in name_field:
            subname, files = item
            result.setdefault(subname, []).extend(files)
        return result

    def _get_chunks(self, sample, list_root_files, tree="Events"):
        """Return list of chunks for a sample, using FilesPerJob."""
        make_job_every_N = sample.get('FilesPerJob', 1)
        return [{'files': list_root_files[i:i + make_job_every_N]}
                for i in range(0, len(list_root_files), make_job_every_N)]

    # _____________________________________________________________________________
    #  C++ code generation — shared parts
    # _____________________________________________________________________________

    def _build_aliases_ordered(self, sampleName):
        """Build aliases_ordered dict preserving aliases.py order."""
        aliases_to_be_defined  = {}
        functions_to_be_defined = {}
        aliases_ordered = {}
        for aliasName, alias in self._aliases.items():
            if 'samples' not in alias.keys() or sampleName in alias['samples']:
                if 'expr' in alias.keys() and 'linesToAdd' not in alias.keys():
                    aliases_to_be_defined[aliasName] = alias['expr']
                    aliases_ordered[aliasName] = ('expr', alias['expr'])
                if ('function' in alias.keys() or 'linesToAdd' in alias.keys()
                        or 'linesToGlobal' in alias.keys()):
                    functions_to_be_defined[aliasName] = alias
                    aliases_ordered[aliasName] = ('function', alias)
        return aliases_to_be_defined, functions_to_be_defined, aliases_ordered

    def _build_cpp_common(self, tree_name, file_path, sampleName, weight,
                          aliases_ordered):
        """
        Returns a dict of C++ fragments shared between mkShapes and mkSnapshot:
          code_of_function_to_include, lines_to_add, lines_to_global, lines_to_process,
          define_input_files_logic, define_aliases (all in order), define_weights,
          booking_logic, define_variables_and_cuts_logic
        """
        # ---- input files ----
        define_input_files_logic = f"//  Nominal input files\n"
        define_input_files_logic += f'    auto* nominal = new TChain("{tree_name}");\n'
        for file_name in file_path:
            define_input_files_logic += f'    nominal->Add("{file_name}");\n'
        define_input_files_logic += f'    \n'

        # mkShapes needs getBranchNames for suffix nuisances; mkSnapshot doesn't
        define_input_files_logic += self._build_cpp_input_files_extra(file_path, sampleName)

        define_input_files_logic += f'    ROOT::RDataFrame base_df(*nominal);\n'
        define_input_files_logic += f'    auto varied_df = ROOT::RDF::RNode(base_df);\n'

        # ---- weights ----
        define_weights = f"//  weights needed\n"
        define_weights += f'    varied_df = SafeDefine(varied_df, "my_sample_weight", "{weight}");\n'

        # ---- aliases + functions in order ----
        code_of_function_to_include = ""
        _included_external_files    = set()
        lines_to_add     = ""
        lines_to_process = ""
        lines_to_global  = ""
        _included_lines_to_add    = set()
        _included_lines_to_global = set()

        define_aliases = "//  aliases and functions (in aliases.py order)\n"

        for aliasName, (kind, payload) in aliases_ordered.items():
            if kind == 'expr':
                _payload = payload.replace("'", '\\"')
                define_aliases += f'    varied_df = SafeDefine(varied_df, "{aliasName}", "{_payload}");\n'
            else:
                alias = payload

                if 'linesToAdd' in alias:
                    for line in alias['linesToAdd']:
                        if line not in _included_lines_to_add:
                            lines_to_add += line + "\n"
                            _included_lines_to_add.add(line)

                if 'linesToProcess' in alias:
                    for line in alias['linesToProcess']:
                        lines_to_process += "    " + line + "\n"

                if 'linesToGlobal' in alias:
                    for line in alias['linesToGlobal']:
                        if line not in _included_lines_to_global:
                            lines_to_global += line + "\n"
                            _included_lines_to_global.add(line)

                if 'external' in alias:
                    if alias["external"] not in _included_external_files:
                        with open(alias["external"], "r") as fh:
                            code_of_function_to_include += fh.read()
                        code_of_function_to_include += "\n"
                        _included_external_files.add(alias["external"])

                if (('linesToAdd' in alias or 'linesToGlobal' in alias or 'linesToProcess' in alias) and 'expr' in alias):
                    _expr = alias['expr'].replace("'", '\\"')
                    define_aliases += f'    varied_df = SafeDefine(varied_df, "{aliasName}", "{_expr}");\n'
                    continue

                define_aliases += "    varied_df = varied_df"
                if 'xmlfile' not in alias:
                    for i, var in enumerate(alias["variables"]):
                        define_aliases += f'\n                     .Define("_var_{aliasName}_{i}", "{var}")'
                else:
                    for i, var in enumerate(alias["variables"]):
                        define_aliases += f'\n                     .Define("_var_{aliasName}_{i}_f", "(float){var}")'
                define_aliases += ";\n"

                if 'xmlfile' not in alias:
                    new_names = [f"_var_{aliasName}_{i}" for i in range(len(alias["variables"]))]
                    sintax_variables = '{ "' + '", "'.join(new_names) + '" }'
                    define_aliases += f'\n    varied_df = varied_df.Define("{aliasName}", {alias["function"]}, {sintax_variables} );\n\n'
                else:
                    submission_dir = os.getcwd()
                    define_aliases += f'\n    TMVA::Experimental::RReader my_model_{aliasName}("{submission_dir}/{alias["xmlfile"]}");\n'
                    new_names = [f"_var_{aliasName}_{i}_f" for i in range(len(alias["variables"]))]
                    sintax_variables = '{ "' + '", "'.join(new_names) + '" }'
                    define_aliases += f'    varied_df = varied_df.Define("{aliasName}",\n'
                    define_aliases += f'                 TMVA::Experimental::Compute<{len(alias["variables"])}, float>(my_model_{aliasName}),\n'
                    define_aliases += f'                 {sintax_variables});\n\n'

        # ---- nuisance variations (mkShapes only — mkSnapshot returns "" here) ----
        register_variations_logic_for_alternative_trees = self._build_cpp_variations_alt_trees(sampleName, file_path, weight)
        register_variations_logic                       = self._build_cpp_variations_weights(sampleName, weight)
        register_variations_for_special_nuisances       = self._build_cpp_variations_special(sampleName)
        code_for_special_nuisances                      = self._build_cpp_code_special_nuisances(sampleName)
        register_maps_special                           = self._build_cpp_maps_special(sampleName)

        # ---- booking_logic ----
        booking_logic = "//  Initial ...\n"
        booking_logic += f'    auto current_node = ROOT::RDF::RNode(varied_df);\n'

        # ---- define variables and cuts ----
        define_variables_logic = "//  Define variables\n"
        for variableName, variable in self._variables.items():
            if not ('is2d' in variable and variable['is2d'] == 1):
                name = variable['name']
                define_variables_logic += f'    current_node = SafeDefine(current_node, "{variableName}", "{name}");\n'

        if self._supercut != '':
            define_variables_logic += f'    current_node = current_node.Filter("{self._supercut}", "supercut");\n'

        # subclass adds cut-level logic (histos for mkShapes, snapshot for mkSnapshot)
        define_variables_logic += self._build_cpp_cuts_and_output(sampleName)

        return {
            'code_of_function_to_include'    : code_of_function_to_include,
            'code_for_special_nuisances'     : code_for_special_nuisances,
            'lines_to_add'                   : lines_to_add,
            'lines_to_global'                : lines_to_global,
            'lines_to_process'               : lines_to_process,
            'define_input_files_logic'       : define_input_files_logic,
            'define_aliases'                 : define_aliases,
            'define_weights'                 : define_weights,
            'register_variations_alt_trees'  : register_variations_logic_for_alternative_trees,
            'register_variations_weights'    : register_variations_logic,
            'register_variations_special'    : register_variations_for_special_nuisances,
            'register_maps_special'          : register_maps_special,
            'booking_logic'                  : booking_logic,
            'define_variables_logic'         : define_variables_logic,
        }

    # These are overridden by mkShapes (nuisances); mkSnapshot returns ""
    def _build_cpp_input_files_extra(self, file_path, sampleName):
        return ""

    def _build_cpp_variations_alt_trees(self, sampleName, file_path, weight):
        return "\n"

    def _build_cpp_variations_weights(self, sampleName, weight):
        return "//  Register the variations\n"

    def _build_cpp_variations_special(self, sampleName):
        return ""

    def _build_cpp_code_special_nuisances(self, sampleName):
        return ""

    def _build_cpp_maps_special(self, sampleName):
        return ""

    # Must be implemented by subclass
    def _build_cpp_cuts_and_output(self, sampleName):
        raise NotImplementedError

    # _____________________________________________________________________________
    #  Proxy check — shared
    # _____________________________________________________________________________

    def _build_cp_cmd(self, out_file, dest_dir, src_prefix=None):
        """Build copy command: xrdcp -f for EOS paths, cp otherwise.
        src_prefix: optional prefix for the source file (e.g. outputDir for mkShapes)
        """
        import re as _re
        if not os.path.isabs(dest_dir):
            dest_dir = os.path.join(os.getcwd(), dest_dir)
        dest_normalized = _re.sub(r'^/eos/home-([a-z])/',
                                  lambda m: f"/eos/user/{m.group(1)}/", dest_dir)
        src = os.path.join(src_prefix, out_file) if src_prefix else out_file
        if dest_dir.startswith('/eos'):
            return f"xrdcp -f {src} root://eosuser.cern.ch//{dest_normalized}/{out_file}"
        return f"cp {src} {dest_dir}/"

    def _check_proxy(self):
        _proxy_candidates = [
            os.environ.get('X509_USER_PROXY', ''),
            os.path.expanduser('~/.proxy'),
            f"/tmp/x509up_u{os.getuid()}",
        ]
        x509_proxy = next((p for p in _proxy_candidates if p and os.path.exists(p)), None)
        if not x509_proxy:
            print("ERROR: no X509 proxy found.")
            print("       Please run: voms-proxy-init --voms cms --valid 168:00")
            return None, None
        proxy_check = subprocess.run(["voms-proxy-info", "--exists", "--valid", "1:00"],
                                     capture_output=True)
        if proxy_check.returncode != 0:
            print(f"ERROR: proxy found at {x509_proxy} but it is expired or invalid.")
            print("       Please run: voms-proxy-init --voms cms --valid 168:00")
            return None, None
        proxy_basename = os.path.basename(x509_proxy)
        proxy_lines = (f"export X509_USER_PROXY=$(pwd)/{proxy_basename}\n"
                       f"export EOS_MGM_URL=root://eosuser.cern.ch\n"
                       f"voms-proxy-info\n")
        print(f"  -> X509 proxy found and valid: {x509_proxy}")
        return x509_proxy, proxy_lines

    # _____________________________________________________________________________
    #  submitBatch — shared skeleton
    # _____________________________________________________________________________

    def submitBatch(self):
        print("======================")
        print("==== submit Batch ====")
        print("======================")

        os.makedirs(self._script_batch_location, exist_ok=True)
        submission_dir = os.getcwd()

        try:
            root_cflags = subprocess.check_output(["root-config", "--cflags"]).decode().strip()
            root_libs   = subprocess.check_output(["root-config", "--libs"]).decode().strip()
        except Exception:
            root_cflags = root_libs = ""

        self._prepare_submit_batch(root_cflags, root_libs)

        precompiled = self._detect_precompiled()
        print(f"  -> {'Pre-compiled' if precompiled else 'Compile-on-node'} mode.")

        x509_proxy, proxy_lines = self._check_proxy()
        if x509_proxy is None:
            return

        queue_rows = []

        for sampleName, sample in self._samples.items():
            os.makedirs(os.path.join(self._script_batch_location, sampleName), exist_ok=True)

            for subname, list_root_files in self._normalize_sample_name(sample['name']).items():
                os.makedirs(os.path.join(self._script_batch_location, sampleName, subname), exist_ok=True)
                if not self._silentMode:
                    os.makedirs(os.path.join(self._script_batch_location, sampleName, subname, "log"), exist_ok=True)

                chunks = self._get_chunks(sample, list_root_files)

                for i, chunk in enumerate(chunks):
                    job_name   = self.job_prefix() + "_" + sampleName + "_" + subname + "_" + str(i)
                    name_bash  = os.path.join(self._script_batch_location, sampleName, subname, job_name + ".sh")
                    out_file   = self.output_file_name(sampleName, subname, i)
                    name_folder      = os.path.join(self._script_batch_location, sampleName, subname)
                    name_folder_code = os.path.join(self._scripts_run_folder, sampleName, subname)
                    abs_executable   = os.path.join(name_folder_code, job_name)
                    abs_cpp          = os.path.join(name_folder_code, job_name + ".cpp")

                    if precompiled:
                        bash_code = self._build_bash_precompiled(
                            job_name, out_file, sampleName, subname, proxy_lines,
                            submission_dir, root_cflags, root_libs)
                        job_exe = abs_executable
                    else:
                        bash_code = self._build_bash_compile_on_node(
                            job_name, out_file, sampleName, subname, proxy_lines,
                            submission_dir, root_cflags, root_libs)
                        job_exe = abs_cpp

                    with open(name_bash, "w") as f:
                        f.write(bash_code)
                    os.system("chmod +x " + name_bash)

                    if self._silentMode:
                        output_file = error_file = log_file = "/dev/null"
                    else:
                        log_dir     = os.path.join(name_folder, "log")
                        output_file = os.path.join(log_dir, job_name + ".sh.out")
                        error_file  = os.path.join(log_dir, job_name + ".sh.err")
                        log_file    = os.path.join(log_dir, job_name + ".sh.log")

                    queue_rows.append((name_bash, job_exe, output_file, error_file, log_file))

        # single .sub file
        global_submit_path = os.path.join(self._script_batch_location, "submit_all.sub")
        extra_files     = self._extra_transfer_files
        extra_files_str = (", " + ", ".join(extra_files)) if extra_files else ""
        fixed_extra = self._extra_transfer_files_for_submit(precompiled, x509_proxy, extra_files_str)

        queue_lines = "\n".join(
            f"  {exe}, {job_exe}, {out}, {err}, {log}"
            for exe, job_exe, out, err, log in queue_rows
        )

        submit_code = (
            f"# Auto-generated by {self.__class__.__name__}\n"
            f"# Mode: {'pre-compiled' if precompiled else 'compile-on-node'}\n\n"
            "should_transfer_files   = YES\n"
            "when_to_transfer_output = ON_EXIT\n"
            "transfer_output_files   = \n"
            f"transfer_input_files    = $(job_exe), {fixed_extra}\n"
            "getenv                  = True\n"
            "+JobFlavour             = \"longlunch\"\n\n"
            "queue executable, job_exe, output, error, log from (\n"
            + queue_lines + "\n)\n"
        )

        with open(global_submit_path, "w") as f:
            f.write(submit_code)

        print(f"Created: {global_submit_path}  ({len(queue_rows)} jobs)")
        print("Submitting all jobs with a single condor_submit ...")
        subprocess.run(["condor_submit", global_submit_path], cwd=self._script_batch_location)

    # Hooks for submitBatch — override in subclass if needed
    def _prepare_submit_batch(self, root_cflags, root_libs):
        pass

    def _detect_precompiled(self):
        for sampleName, sample in self._samples.items():
            for subname, list_root_files in self._normalize_sample_name(sample['name']).items():
                if list_root_files:
                    exe = os.path.join(self._scripts_run_folder, sampleName, subname,
                                       self.job_prefix() + "_" + sampleName + "_" + subname + "_0")
                    return os.path.isfile(exe)
        return False

    # Abstract — subclass must implement
    def job_prefix(self):
        raise NotImplementedError

    def output_file_name(self, sampleName, subname, i):
        raise NotImplementedError

    def _extra_transfer_files_for_submit(self, precompiled, x509_proxy, extra_files_str):
        raise NotImplementedError

    def _build_bash_precompiled(self, job_name, out_file, sampleName, subname,
                                proxy_lines, submission_dir, root_cflags, root_libs):
        raise NotImplementedError

    def _build_bash_compile_on_node(self, job_name, out_file, sampleName, subname,
                                    proxy_lines, submission_dir, root_cflags, root_libs):
        raise NotImplementedError

    # _____________________________________________________________________________
    #  checkBatch — shared
    # _____________________________________________________________________________

    def checkBatch(self, resubmit=False):
        print("=====================")
        print("==== check Batch ====")
        print("=====================\n")

        stats      = {}
        missing_sh = []

        for sampleName in self._samples:
            sample_batch_dir = os.path.join(self._script_batch_location, sampleName)
            if not os.path.isdir(sample_batch_dir):
                continue
            stats[sampleName] = {'done': 0, 'error': 0, 'missing': 0}

            for subname in os.listdir(sample_batch_dir):
                subname_dir = os.path.join(sample_batch_dir, subname)
                if not os.path.isdir(subname_dir):
                    continue
                log_dir  = os.path.join(subname_dir, "log")
                sh_files = sorted(f for f in os.listdir(subname_dir) if f.endswith('.sh'))

                for sh_file in sh_files:
                    job_name  = sh_file[:-3]
                    abs_sh    = os.path.abspath(os.path.join(subname_dir, sh_file))
                    suffix    = job_name.replace(self.job_prefix() + "_", "", 1)
                    _out      = self._output_file_from_suffix(sampleName, subname, suffix)
                    root_file = _out if os.path.isabs(_out) else os.path.join(self._outputDir, _out)

                    name_err = os.path.join(log_dir, job_name + ".sh.err")
                    name_log = os.path.join(log_dir, job_name + ".sh.log")

                    if os.path.exists(root_file):
                        stats[sampleName]['done'] += 1
                    else:
                        missing_sh.append(abs_sh)
                        if not os.path.exists(name_err):
                            if os.path.exists(name_log):
                                with open(name_log, 'r') as f:
                                    content = f.read()
                                if 'Job was aborted' in content:
                                    stats[sampleName]['error'] += 1
                                else:
                                    stats[sampleName]['missing'] += 1
                            else:
                                stats[sampleName]['missing'] += 1
                        else:
                            with open(name_err, 'r') as f:
                                lines = f.readlines()
                            if len(lines) > 6:
                                stats[sampleName]['error'] += 1
                            else:
                                stats[sampleName]['missing'] += 1

        if not stats:
            print("No batch scripts found. Did you run --submitBatch?")
            return

        col_w = max(len(s) for s in stats.keys()) + 2
        print(f"{'Sample':<{col_w}} {'Done':>8} {'Missing':>8} {'Error':>8}")
        print("-" * (col_w + 28))
        total = {'done': 0, 'error': 0, 'missing': 0}
        for sampleName, s in stats.items():
            print(f"{sampleName:<{col_w}} {s['done']:>8} {s['missing']:>8} {s['error']:>8}")
            for k in total:
                total[k] += s[k]
        print("-" * (col_w + 28))
        print(f"{'TOTAL':<{col_w}} {total['done']:>8} {total['missing']:>8} {total['error']:>8}")

        if not missing_sh:
            print("\nAll jobs done!")
            return

        print(f"\nMissing/failed jobs: {len(missing_sh)}")
        if not resubmit:
            print("Run with --resubmit to generate submit_missing.sub and resubmit.")
            return

        # build submit_missing.sub
        submit_all_path = os.path.join(self._script_batch_location, "submit_all.sub")
        if not os.path.exists(submit_all_path):
            print(f"ERROR: {submit_all_path} not found, cannot resubmit.")
            return

        with open(submit_all_path, 'r') as f:
            all_lines = f.readlines()

        header_lines = []
        for line in all_lines:
            header_lines.append(line)
            if line.strip().startswith("queue executable"):
                break

        exe_to_entry = {}
        in_queue = False
        for line in all_lines:
            if not in_queue:
                if line.strip().startswith("queue executable"):
                    in_queue = True
                continue
            stripped = line.strip()
            if not stripped or stripped == ')':
                continue
            exe_path = os.path.abspath(stripped.split(',')[0].strip())
            exe_to_entry[exe_path] = line

        filtered_entries = []
        unmatched        = []
        missing_sh_set   = set(missing_sh)
        for sh in missing_sh:
            abs_sh = os.path.abspath(sh)
            if abs_sh in exe_to_entry:
                filtered_entries.append(exe_to_entry[abs_sh])
            else:
                unmatched.append(abs_sh)

        for abs_sh in unmatched:
            sh_basename = os.path.basename(abs_sh)
            job_name    = sh_basename[:-3]
            parts       = abs_sh.split(os.sep)
            subname_    = parts[-2]
            sampleName_ = parts[-3]
            name_folder_code = os.path.join(self._scripts_run_folder, sampleName_, subname_)
            abs_exe     = os.path.join(name_folder_code, job_name)
            abs_cpp_    = os.path.join(name_folder_code, job_name + ".cpp")
            job_exe     = abs_exe if os.path.exists(abs_exe) else abs_cpp_
            log_dir_    = os.path.join(os.path.dirname(abs_sh), "log")
            out_f       = os.path.join(log_dir_, job_name + ".sh.out")
            err_f       = os.path.join(log_dir_, job_name + ".sh.err")
            log_f       = os.path.join(log_dir_, job_name + ".sh.log")
            filtered_entries.append(f"  {abs_sh}, {job_exe}, {out_f}, {err_f}, {log_f}\n")

        if not filtered_entries:
            print("WARNING: could not build any entries for submit_missing.sub.")
            return

        submit_missing_path = os.path.join(self._script_batch_location, "submit_missing.sub")
        with open(submit_missing_path, 'w') as f:
            f.writelines(header_lines)
            f.writelines(filtered_entries)
            f.write(")\n")

        print(f"\nCreated: {submit_missing_path}  ({len(filtered_entries)} jobs)")
        print("Submitting missing jobs ...")
        subprocess.run(["condor_submit", submit_missing_path], cwd=self._script_batch_location)

    def _output_file_from_suffix(self, sampleName, subname, suffix):
        """Map job suffix -> expected output root file name. Override in subclass."""
        raise NotImplementedError

    # _____________________________________________________________________________
    #  haddBatch — shared, submit a single condor job to hadd all output files
    # _____________________________________________________________________________

    def _check_all_done(self):
        """
        Silent check: returns (all_done, missing_list).
        Does not print a table — just counts missing files.
        """
        missing = []
        for sampleName in self._samples:
            sample_batch_dir = os.path.join(self._script_batch_location, sampleName)
            if not os.path.isdir(sample_batch_dir):
                continue
            for subname in os.listdir(sample_batch_dir):
                subname_dir = os.path.join(sample_batch_dir, subname)
                if not os.path.isdir(subname_dir):
                    continue
                sh_files = sorted(f for f in os.listdir(subname_dir) if f.endswith('.sh'))
                for sh_file in sh_files:
                    job_name = sh_file[:-3]
                    suffix   = job_name.replace(self.job_prefix() + "_", "", 1)
                    _out     = self._output_file_from_suffix(sampleName, subname, suffix)
                    root_file = _out if os.path.isabs(_out) else os.path.join(self._outputDir, _out)
                    if not os.path.exists(root_file):
                        missing.append(root_file)
        return len(missing) == 0, missing

    def haddBatch(self, tag, output_dir, glob_pattern, output_name,
                  hadd_tool_path, hadd_tool_args="-f -j 8 -E"):
        """
        Submit a single condor job that hadd all files matching glob_pattern in output_dir.
        hadd_tool_path: path to haddfast.py or haddnano.py
        hadd_tool_args: extra args to pass to the tool
        """
        import glob as _glob
        import re as _re

        submission_dir = os.getcwd()

        # silent pre-check: abort if any output files are missing
        all_done, missing = self._check_all_done()
        if not all_done:
            print(f"ERROR: {len(missing)} output files missing — cannot hadd incomplete results.")
            print("Missing files:")
            for f in missing:
                print(f"  {f}")
            print("Run --checkBatch to get the full picture, --resubmit to resubmit missing jobs.")
            return
        print(f"  -> All jobs done, proceeding with hadd ({len(missing)} missing = 0)")

        x509_proxy, proxy_lines = self._check_proxy()
        if x509_proxy is None:
            return

        # resolve absolute output_dir
        if not os.path.isabs(output_dir):
            output_dir = os.path.join(submission_dir, output_dir)
        output_dir_normalized = _re.sub(r'^/eos/home-([a-z])/',
                                         lambda m: f"/eos/user/{m.group(1)}/", output_dir)

        # find all files matching pattern
        all_files = sorted(_glob.glob(os.path.join(output_dir, glob_pattern)))
        if not all_files:
            print(f"ERROR: no files found matching {os.path.join(output_dir, glob_pattern)}")
            return
        print(f"  -> Found {len(all_files)} files to hadd")

        output_eos = f"root://eosuser.cern.ch//{output_dir_normalized}/{output_name}"

        # build input file list with xrootd paths
        input_files_str = " \\\n    ".join(
            f"root://eosuser.cern.ch//{output_dir_normalized}/{os.path.basename(f)}"
            for f in all_files
        )

        os.makedirs(self._script_batch_location, exist_ok=True)

        hadd_sh_path  = os.path.join(self._script_batch_location, f"hadd_{tag}.sh")
        hadd_sub_path = os.path.join(self._script_batch_location, f"hadd_{tag}.sub")

        hadd_tool_basename = os.path.basename(hadd_tool_path)
        # if hadd tool is on EOS, download it in the bash script instead of transfer_input_files
        hadd_tool_abs = os.path.abspath(hadd_tool_path) if not os.path.isabs(hadd_tool_path) else hadd_tool_path
        if hadd_tool_abs.startswith('/eos'):
            import re as _re
            hadd_tool_normalized = _re.sub(r'^/eos/home-([a-z])/',
                                           lambda m: f"/eos/user/{m.group(1)}/", hadd_tool_abs)
            hadd_tool_download = f"xrdcp -f root://eosuser.cern.ch//{hadd_tool_normalized} {hadd_tool_basename}"
            hadd_tool_transfer = None
        else:
            hadd_tool_download = ""
            hadd_tool_transfer = hadd_tool_path

        hadd_sh = f"""#!/bin/bash
set -e
echo "hadd job started at $(date)"
echo "Running on node $(hostname)"
source /cvmfs/sft.cern.ch/lcg/views/LCG_109a/x86_64-el9-gcc13-opt/setup.sh
{proxy_lines}
{hadd_tool_download}
echo "Running {hadd_tool_basename}..."
python {hadd_tool_basename} \\
    {hadd_tool_args} \\
    {output_name} \\
    {input_files_str}

echo "hadd done at $(date)"
xrdcp -f {output_name} {output_eos}
rm -f {output_name}
echo "Done at $(date)"; ls -lh
"""
        with open(hadd_sh_path, "w") as f:
            f.write(hadd_sh)
        os.system(f"chmod +x {hadd_sh_path}")

        transfer_files = x509_proxy
        if hadd_tool_transfer:
            transfer_files = f"{hadd_tool_transfer}, {x509_proxy}"

        hadd_sub = f"""# hadd job for {tag}
executable              = {hadd_sh_path}
transfer_input_files    = {transfer_files}
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
transfer_output_files   =
output                  = {self._script_batch_location}/hadd_{tag}.out
error                   = {self._script_batch_location}/hadd_{tag}.err
log                     = {self._script_batch_location}/hadd_{tag}.log
getenv                  = True
+JobFlavour             = "workday"
queue
"""
        with open(hadd_sub_path, "w") as f:
            f.write(hadd_sub)

        print(f"  -> hadd script : {hadd_sh_path}")
        print(f"  -> submit file : {hadd_sub_path}")
        print(f"  -> output      : {output_eos}")
        print("Submitting hadd job...")
        subprocess.run(["condor_submit", hadd_sub_path],
                       cwd=self._script_batch_location)

    # _____________________________________________________________________________
    #  parallelCompile — shared
    # _____________________________________________________________________________

    def parallelCompile(self):
        print("==========================")
        print("==== compile on batch ====")
        print("==========================")

        where_compile  = "parallel_compile"
        os.makedirs(where_compile, exist_ok=True)
        submission_dir = os.getcwd()

        for sampleName in self._samples:
            name_bash   = os.path.join(where_compile, f"my_script_compile_{sampleName}.sh")
            name_submit = os.path.join(where_compile, f"my_script_compile_{sampleName}.sub")

            bash_code = f"""#!/bin/bash
set -e
echo "Job started at $(date)"
echo "Running on node $(hostname)"
cd {submission_dir}
make {sampleName}
"""
            with open(name_bash, "w") as f:
                f.write(bash_code)
            os.system("chmod +x " + name_bash)

            submit_code = f"""
initialdir            = {submission_dir}
executable            = {name_bash}
error                 = /dev/null
log                   = /dev/null
getenv                = True
+JobFlavour           = "longlunch"
queue
"""
            with open(name_submit, "w") as f:
                f.write(submit_code)

            print("Submitting job to HTCondor: condor_submit " + name_submit)
            subprocess.run(["condor_submit", name_submit])