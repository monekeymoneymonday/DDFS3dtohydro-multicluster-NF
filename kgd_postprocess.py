import json, glob, re, math, os
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------------------
# Collect state files produced by HF_driver (JSON mode)
# -------------------------------------------------------------
files = sorted(glob.glob('state_*.json'))
if not files:
    files = ['state.json']  # single combined file

times, w_max, p_well = [], [], []
for fname in files:
    with open(fname) as f:
        data = json.load(f)
    step = data['step']
    t = step  # assumes dt=1; adjust later if needed
    times.append(t)
    pres = np.array(data['pressure'])
    aper = np.array(data['aperture'])
    w_max.append(aper.max())
    p_well.append(pres[0])

times = np.array(times)
w_max = np.array(w_max)
p_well= np.array(p_well)

# -------------------------------------------------------------
# Analytical KGD (viscous) scaling: w_max = A * t^0.25
# uses user-provided constants
# -------------------------------------------------------------
E   = 2.5e10
nu  = 0.25
E_p = E/(1-nu**2)
mu  = 1.0e-3
Q0  = 1.0e-4  # change as appropriate
A   = 2.0 * (mu*Q0/E_p)**0.25
w_ref = A * times**0.25

plt.figure(figsize=(6,4))
plt.loglog(times, w_max, 'o', label='Simulation')
plt.loglog(times, w_ref, '--', label='KGD viscous')
plt.xlabel('time [s]')
plt.ylabel('maximum aperture [m]')
plt.legend()
plt.title('KGD benchmark')
plt.tight_layout()
plt.savefig('kgd_width.png')

# -------------------------------------------------------------
# Carter leak-off curve (if C_L>0): plot average leak vs sqrt(time)
# assumes leak array saved separately (extend if needed)
# -------------------------------------------------------------
if 'leak' in data:
    leak_rate = np.array([np.mean(json.load(open(f))['leak']) for f in files])
    plt.figure()
    plt.plot(np.sqrt(times), leak_rate, 'o-')
    plt.xlabel('sqrt(time)')
    plt.ylabel('avg leak-off rate [m^2/s]')
    plt.title('Carter leak-off check (should be ~1/sqrt(t))')
    plt.tight_layout()
    plt.savefig('carter_leakoff.png')

print('Post-processing complete: kgd_width.png (and carter_leakoff.png if applicable) saved.') 