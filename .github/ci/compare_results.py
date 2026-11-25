"""Compare replay stdout files with results/results.txt (a transcript: '$ python3 -I -B -O <script>.py', stdout, 'exit=0').
Policy (declared before tagging): the three standard-library certificates (verify_certificate, verify_lag_bins_exact,
verify_moderate_speed) must match exactly. The four floating-point scripts must match token by token: non-numeric
tokens exactly; numeric tokens within relative 1e-6 or absolute 1e-9. Interpreter versions are not compared.
usage: compare_results.py RESULTS.txt OUTDIR [SCRIPT ...]  (OUTDIR holds <script>.out files; default: all seven)"""
# SPDX-License-Identifier: MIT
import re, sys
from pathlib import Path
EXACT = ['verify_certificate', 'verify_lag_bins_exact', 'verify_moderate_speed']
FLOAT = ['verify_numerically', 'verify_histogram', 'verify_crossover', 'verify_pumped_current']
REL, ABS = 1e-6, 1e-9
sec, cur = {}, None
for line in Path(sys.argv[1]).read_text().splitlines():
    m = re.match(r'^\$ python3 -I -B -O (\w+)\.py$', line)
    if m: cur = m.group(1); sec[cur] = []; continue
    if cur is None: continue
    if line == 'exit=0': cur = None; continue
    sec[cur].append(line)
if sorted(sec) != sorted(EXACT + FLOAT): sys.exit(f'FAIL sections {sorted(sec)}')
out, bad = Path(sys.argv[2]), []
only = sys.argv[3:] or EXACT + FLOAT
if set(only) - set(EXACT + FLOAT): sys.exit(f'FAIL unknown script {only}')
def num(t):
    try: return float(t)
    except ValueError: return None
for s in [x for x in EXACT + FLOAT if x in only]:
    got = (out / f'{s}.out').read_text().splitlines()
    exp = sec[s]
    if s in EXACT:
        if got != exp: bad.append(f'{s}: exact output differs')
        else: print(f'{s}: identical ({len(exp)} lines)')
        continue
    if len(got) != len(exp): bad.append(f'{s}: line count {len(got)} vs {len(exp)}'); continue
    worst = 0.0
    for a, b in zip(got, exp):
        ta, tb = a.split(), b.split()
        if len(ta) != len(tb): bad.append(f'{s}: token count differs: {a!r} vs {b!r}'); break
        for x, y in zip(ta, tb):
            fx, fy = num(x), num(y)
            if fx is None or fy is None:
                if x != y: bad.append(f'{s}: token {x!r} vs {y!r}')
                continue
            d = abs(fx - fy)
            if d > max(ABS, REL * abs(fy)): bad.append(f'{s}: {x} vs {y}')
            if abs(fy) > ABS: worst = max(worst, d / abs(fy))
    print(f'{s}: within tolerance (largest relative difference {worst:.2e}; {"identical" if got == exp else "not byte-identical"})')
if bad: print('FAIL', bad); sys.exit(1)
print(f'PASS {len(only)} output(s) agree with results/results.txt under the declared policy')
