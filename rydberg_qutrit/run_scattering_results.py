"""Reproduce intermediate-state scattering sensitivity and export Kraus matrices."""
from pathlib import Path
import json
import numpy as np
from rydberg_scattering import (
    averaged_channel, analyze, ScatteringParameters, export_results
)
from rydberg_full_channel import construct_channel

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'scattering_results'
RESULTS.mkdir(exist_ok=True)

print('Calculating decay + Gaussian detuning baseline ...',flush=True)
baseline=analyze(averaged_channel(include_scattering=False,t2_star_us=3.0,quadrature_order=7))
print('Calculating decay + detuning + intermediate-state scattering ...',flush=True)
model=ScatteringParameters(target='1')
scattered=analyze(averaged_channel(scattering=model,t2_star_us=3.0,quadrature_order=7))

m={'source':'Evered 2023 for laser & Rydberg parameters; 112.8ns Rb 6P3/2 lifetime from earlier lifetime measurements',
   'scattering_model':'leading-order adiabatic elimination, ASSUMED scattering returns to |1>, common-mode Gaussian T2*=3us',
   'intermediate_state_lifetime_us':model.intermediate_lifetime_us,
   'blue_rabi_mhz':model.omega_blue_mhz,'red_rabi_mhz':model.omega_red_mhz,
   'intermediate_detuning_mhz':model.intermediate_detuning_mhz,
   'warning':'This is NOT an experimentally calibrated branching model or a complete experimental fidelity prediction.'}
obj=export_results(scattered,RESULTS,m)
for label,item in [('without scattering',baseline),('with scattering',scattered)]:
    print(label, 'Kraus count',len(item['Kraus']),
          'identity',item['pauli']['II'], 'nonid',item['nonidentity'],
          'rydberg',item['outcomes']['rydberg'],
          'loss',item['outcomes']['any_loss'],flush=True)
    print(' verification',item['normalization'],item['trace_error'],item['reconstruction_error'],flush=True)

simple=lambda x: {'pauli':x['pauli'], 'pauli_nonidentity':x['nonidentity'],
    'rydb_residual':x['outcomes']['rydberg'], 'loss':x['outcomes']['any_loss'],
    'comp':x['outcomes']['computational'], 'kraus_count':len(x['Kraus'])}
summary={'baseline':simple(baseline),'scattering':simple(scattered),
         'model':m}
(RESULTS/'comparison.json').write_text(json.dumps(summary,indent=2))

# Largest norm Kraus: 16 x 4, display computational 4x4 block for compact inspection.
ks=scattered['Kraus']; orders=np.argsort(np.linalg.norm(ks,axis=(1,2)))[::-1]
np.savetxt(RESULTS/'dominant_kraus_computational_block.txt',
           ks[orders[0]][[0,1,4,5],:].view(float).reshape(4,8),fmt='%.12g',
           header='real imag pairs for 4 columns in order |00>, |01>, |10>, |11>')
print('Most significant Kraus index',orders[0], 'norm',np.linalg.norm(ks[orders[0]]))
print('Top comp block',ks[orders[0]][[0,1,4,5],:])
print('Done, exported:',RESULTS,flush=True)
