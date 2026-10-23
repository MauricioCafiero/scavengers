import subprocess, os
os.chdir('/Users/cafierom/python_mac/peptidebuilder')
for name, mol in (('rep6_oct_1','octinoxate'), ('rep6_oct_2','octinoxate'),
                  ('rep6_oct_3','octinoxate'), ('rep6_oxy','oxybenzone')):
    open('runs/folds_md/launch_one_drain.py','w').write(
        "import subprocess, os\n"
        f"e = os.environ.copy(); e['SYSNAME'] = '{mol}'\n"
        f"subprocess.Popen(['zsh', 'code/racc_drain.sh', 'runs/folds_md/{name}.legs'],\n"
        f"  stdout=open('runs/folds_md/{name}_drain.log','ab'), stderr=subprocess.STDOUT,\n"
        "  stdin=subprocess.DEVNULL, start_new_session=True,\n"
        "  cwd='/Users/cafierom/python_mac/peptidebuilder', env=e)\n")
    subprocess.run(['/Users/cafierom/miniforge3/bin/python', 'runs/folds_md/launch_one_drain.py'],
                   cwd='/Users/cafierom/python_mac/peptidebuilder')
    subprocess.run(['sleep', '2'])
print('launched 4 rep6 drains')