"""Produce publication-usable comparison figures from saved numerical data."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here=Path(__file__).resolve().parent
baseline=json.loads((here/'noise_diagnostics_results.json').read_text())['decay_plus_T2star']
with_scatter=json.loads((here/'scattering_results/scattering_channel_results.json').read_text())
p0=baseline['q_P'];p1=with_scatter['pauli_probabilities']
fig,ax=plt.subplots(figsize=(11.5,5.3))
labels=[p for p in sorted(p0) if p!='II']
x=np.arange(len(labels));w=0.42
ax.bar(x-w/2, [100*p0[p] for p in labels],w,label='Rydberg decay + T2*')
ax.bar(x+w/2, [100*p1[p] for p in labels],w,label='+ intermediate-state scattering (assumed return to |1>)')
ax.set_yscale('log');ax.set_ylim(5e-7,0.1)
ax.set_xticks(x,labels);ax.set_xlabel('Two-qubit Pauli string')
ax.set_ylabel('Pauli-twirled weight (%)')
ax.set_title('Effect of virtual intermediate-state scattering on the Pauli distribution')
ax.legend(fontsize=8,loc='upper center');ax.grid(axis='y',alpha=.25)
fig.tight_layout();fig.savefig(here/'scattering_results/scattering_pauli_comparison.png',dpi=180)
fig.savefig(here/'scattering_results/scattering_pauli_comparison.svg');plt.close(fig)

fig,ax=plt.subplots(figsize=(9,4.9))
series=[('Decay + T2*',p0['II'], baseline['p_nonidentity_pauli'],baseline['p_rydberg_avg'],baseline['p_loss_avg']),
 ('+ scattering',p1['II'],with_scatter['nonidentity_pauli'],with_scatter['outcomes']['rydberg'],with_scatter['outcomes']['any_loss'])]
x=np.arange(len(series))
for j, lab in enumerate(['Nonidentity Pauli', 'Residual Rydberg', 'Loss state']):
 values=[100*s[j+2] for s in series]
 ax.bar(x+(j-1)*0.24,values,0.23,label=lab)
ax.set_xticks(x,[s[0] for s in series]); ax.set_ylabel('Input-averaged channel weight (%)')
ax.set_title('Error budget excluding the dominant identity contribution')
ax.legend();ax.grid(axis='y',alpha=.25);fig.tight_layout()
fig.savefig(here/'scattering_results/scattering_error_budget.png',dpi=180)
fig.savefig(here/'scattering_results/scattering_error_budget.svg');plt.close(fig)

comparison={
 'baseline_decay_plus_T2star': {
  'identity':p0['II'],'nonidentity':baseline['p_nonidentity_pauli'],
  'rydberg':baseline['p_rydberg_avg'],'loss':baseline['p_loss_avg']},
 'with_intermediate_scattering_return_to_1': {
  'identity':p1['II'],'nonidentity':with_scatter['nonidentity_pauli'],
  'rydberg':with_scatter['outcomes']['rydberg'],'loss':with_scatter['outcomes']['any_loss']},
 'scattering_target_not_calibrated':True,
 'total_change_in_nonidentity_pp':100*(with_scatter['nonidentity_pauli']-baseline['p_nonidentity_pauli'])
}
(here/'scattering_results/comparison.json').write_text(json.dumps(comparison,indent=2))
print(json.dumps(comparison,indent=2))
