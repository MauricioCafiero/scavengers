import subprocess, time, os
os.chdir('/Users/cafierom/python_mac/peptidebuilder')
for line in open('runs/folds_md/rep5.chain'):
    parts = line.split()
    if len(parts) != 3:
        continue
    name, mol = parts[0], parts[1]
    open('runs/folds_md/launch_one_drain.py','w').write(
        "import subprocess, os\n"
        f"e = os.environ.copy(); e['SYSNAME'] = '{mol}'\n"
        f"subprocess.Popen(['zsh', 'code/racc_drain.sh', 'runs/folds_md/{name}.legs'],\n"
        f"  stdout=open('runs/folds_md/{name}_drain.log','ab'), stderr=subprocess.STDOUT,\n"
        "  stdin=subprocess.DEVNULL, start_new_session=True,\n"
        "  cwd='/Users/cafierom/python_mac/peptidebuilder', env=e)\n")
    subprocess.run(['/Users/cafierom/miniforge3/bin/python', 'runs/folds_md/launch_one_drain.py'],
                   cwd='/Users/cafierom/python_mac/peptidebuilder')
    time.sleep(3)
print('launched drains for', len(open('runs/folds_md/rep5.chain').read().strip().split(chr(10))), 'batches')