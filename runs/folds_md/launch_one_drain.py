import subprocess, os
e = os.environ.copy(); e['SYSNAME'] = 'oxybenzone'
subprocess.Popen(['zsh', 'code/racc_drain.sh', 'runs/folds_md/rep6_oxy.legs'],
  stdout=open('runs/folds_md/rep6_oxy_drain.log','ab'), stderr=subprocess.STDOUT,
  stdin=subprocess.DEVNULL, start_new_session=True,
  cwd='/Users/cafierom/python_mac/peptidebuilder', env=e)
