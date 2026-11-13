"""Generate full-channel numerical summary and nonidentity Pauli/leakage plot."""
import json
from pathlib import Path
import matplotlib.pyplot as plt
from rydberg_full_channel import channel_report

out=Path(__file__).parent
r=channel_report()
(out/"channel_results.json").write_text(json.dumps(r,indent=2)+"\n")
rows=[(name,100*prob) for name,prob in r["pauli_weights"].items() if name!="II"]
rows.sort(key=lambda x:-x[1])
rows += [("Residual |r>",100*r["p_rydberg_residual"]),
         ("Loss sink",100*r["p_loss"])]
fig,ax=plt.subplots(figsize=(10,5))
ax.bar(range(len(rows)),[p for _,p in rows])
ax.set_xticks(range(len(rows)),[k for k,_ in rows],rotation=60,ha="right")
ax.set_yscale("log")
ax.set_ylim(1e-7,1)
ax.set_ylabel("Probability (%)")
ax.set_title("Reduced two-atom Rydberg model: nonidentity Pauli and leakage weights")
fig.tight_layout()
fig.savefig(out/"rydberg_full_channel_probabilities.png",dpi=180)
plt.close(fig)
print(json.dumps(r,indent=2))
