#!/bin/zsh
# Guarded launcher for every local `modal run` of code/modal_md.py.
#
# caffeinate holds an idle-sleep assertion for exactly the utility's lifetime,
# so the Mac is held awake for the whole modal client lifetime. The client stays
# attached for the whole leg even with --detach, and Modal cancels its in-flight
# leg the moment this Mac sleeps and the client drops: on 2026-10-04 two
# produce legs were cancelled this way for $1.06 billed with no result.
#
# Usage (from the repo root):
#   code/modal_run.sh --detach code/modal_md.py::produce --structure <structure>
exec caffeinate -i /opt/homebrew/bin/modal run "$@"