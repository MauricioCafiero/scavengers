"""Work out how to run Boltz-2 on this machine, rather than assuming where it lives.

The scripts here were written against one particular setup — a Boltz checkout at
`~/python_mac/boltz_local` with an MPS wrapper in it — and hardwiring that makes them useless to
anyone else. This resolves the command to run in order of preference, so a plain `pip install boltz`
works, a separate environment works, and the local wrapper still works if it is there.

Resolution order, first hit wins:

1. `--boltz-cmd` on the command line, or `$PEPTIDEBUILDER_BOLTZ_CMD`: a complete command, used as
   given. The escape hatch for anything unusual.
2. `--boltz-venv` or `$PEPTIDEBUILDER_BOLTZ_VENV`: a virtualenv holding Boltz. If a sibling
   `code/boltz_mps.py` wrapper sits beside it, that is used (it splits the structure and affinity
   passes so they do not share GPU memory); otherwise the venv's own `boltz` entry point.
3. `boltz` importable in the environment running this script — the ordinary `pip install boltz` case.
4. `boltz` on `$PATH`.
5. The original developer's layout, if it happens to exist, so nothing that worked before breaks.

If none match, the error says what to do rather than failing obscurely.
"""
import os
import shutil
import subprocess
import sys

LEGACY_REPO = os.path.expanduser("~/python_mac/boltz_local")


def _wrapper_for(venv):
    """The MPS wrapper beside a venv, if that layout is present."""
    repo = os.path.dirname(os.path.normpath(venv))
    wrapper = os.path.join(repo, "code", "boltz_mps.py")
    return wrapper if os.path.exists(wrapper) else None


def _venv_python(venv):
    python = os.path.join(venv, "bin", "python")
    return python if os.path.exists(python) else None


def resolve_boltz(boltz_cmd=None, boltz_venv=None, verbose=True):
    """The command prefix that runs a Boltz prediction, as a list, plus a label saying how.

    Append `["predict", yaml, ...]` to what comes back.
    """
    explicit = boltz_cmd or os.environ.get("PEPTIDEBUILDER_BOLTZ_CMD")
    if explicit:
        return explicit.split(), "explicit command"

    venv = boltz_venv or os.environ.get("PEPTIDEBUILDER_BOLTZ_VENV")
    if venv:
        venv = os.path.expanduser(venv)
        python = _venv_python(venv)
        if not python:
            sys.exit(f"no python in {venv}: check --boltz-venv / $PEPTIDEBUILDER_BOLTZ_VENV")
        wrapper = _wrapper_for(venv)
        if wrapper:
            return [python, wrapper], f"wrapper {wrapper}"
        return [python, "-m", "boltz"], f"boltz in {venv}"

    try:
        import boltz  # noqa: F401
        return [sys.executable, "-m", "boltz"], "boltz in this environment"
    except ImportError:
        pass

    found = shutil.which("boltz")
    if found:
        return [found], f"boltz on PATH ({found})"

    python = _venv_python(os.path.join(LEGACY_REPO, ".venv"))
    wrapper = os.path.join(LEGACY_REPO, "code", "boltz_mps.py")
    if python and os.path.exists(wrapper):
        return [python, wrapper], f"wrapper {wrapper}"

    sys.exit(
        "cannot find Boltz-2. Any one of these will do:\n"
        "  pip install boltz                          (into this environment)\n"
        "  export PEPTIDEBUILDER_BOLTZ_VENV=/path/to/venv\n"
        "  --boltz-venv /path/to/venv\n"
        "  --boltz-cmd 'python -m boltz'              (or whatever runs it for you)\n"
        "  export PEPTIDEBUILDER_BOLTZ_CMD='...'"
    )


def run_boltz(yaml_path, out_dir, log_path, extra_args=(), boltz_cmd=None, boltz_venv=None,
              verbose=True):
    """Run one Boltz prediction. Returns the process exit code."""
    prefix, how = resolve_boltz(boltz_cmd, boltz_venv, verbose)
    cmd = list(prefix) + ["predict", yaml_path, "--num_workers", "0", "--out_dir", out_dir,
                          *extra_args]
    if verbose:
        print(f"  Boltz via {how}", flush=True)
    with open(log_path, "w") as log:
        return subprocess.call(cmd, stdout=log, stderr=subprocess.STDOUT)


def add_arguments(parser):
    """The two flags every script that folds should accept."""
    parser.add_argument("--boltz-venv",
                        help="virtualenv containing Boltz-2 (or $PEPTIDEBUILDER_BOLTZ_VENV). "
                             "Not needed if boltz is installed here or on PATH")
    parser.add_argument("--boltz-cmd",
                        help="complete command that runs Boltz, e.g. 'python -m boltz' "
                             "(or $PEPTIDEBUILDER_BOLTZ_CMD)")
