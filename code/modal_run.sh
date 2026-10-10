#!/bin/zsh
# Guarded launcher for every local `modal run` of code/modal_md.py.
#
# caffeinate holds sleep assertions for exactly the utility's lifetime,
# so the Mac is held awake for the whole modal client lifetime. The client stays
# attached for the whole leg even with --detach, and Modal cancels its in-flight
# leg the moment this Mac sleeps and the client drops: on 2026-10-04 two
# produce legs were cancelled this way for $1.06 billed with no result.
# -i (prevent idle) alone proved insufficient on 2026-10-10: a 16-minute
# maintenance/dark sleep (03:39:53-03:55:54) dropped the produce client mid-leg
# and its disconnect handler cancelled a 30-minute-in produce. -s additionally
# asserts prevention of system sleep and is valid on AC power (this Mac runs
# on AC for every leg).
#
# Usage (from the repo root):
#   code/modal_run.sh --detach code/modal_md.py::produce --structure <structure>
exec caffeinate -i -s /opt/homebrew/bin/modal run "$@"