"""Draw the analytic rate certificate and limiting waiting-time densities."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).resolve().parent
theta = 2 * np.pi * np.arange(2048) / 2048
q = 1 + 16 / 17 * np.cos(theta)
d = (np.sin(2 * theta) - np.cos(2 * theta)) / 5
u = np.linspace(0, 4, 300)
rates = np.geomspace(1500, 120_000, 350)
sigma = 4 * np.pi / 5
error = 2 * np.pi * (18 / 5 + 32 / 17) / rates
bound = rates * np.maximum(0, 17 / 600 - error) ** 2
fig, axs = plt.subplots(1, 2, figsize=(7.12, 2.65), layout="constrained")
fig.patch.set_facecolor("white")
brown, jade = "#923E39", "#177266"
for mass, col, name in [(0.5 + d, jade, "Forward drive"),
                        (0.5 - d, brown, "Reversed drive")]:
    den = (2 * mass * q**2 * np.exp(-u[:, None] * q)).mean(axis=1)
    axs[0].plot(u, den, color=col, lw=2, label=name)
axs[0].set(xlabel=r"Rescaled waiting time $u=R\tau$", ylabel="Probability density",
           xlim=(0, 4), ylim=(0, 1.7), title="A   Waiting times at high activity")
axs[0].legend(frameon=False, fontsize=8)
axs[1].plot(rates, bound, color=jade, lw=2, label="Proved lower bound on estimate")
axs[1].axhline(sigma, color=brown, ls="--", lw=1.6, label="True entropy production")
axs[1].scatter([1e4,1e5], [r * max(0, 17/600 - 2*np.pi*(18/5+32/17)/r)**2
                           for r in [1e4,1e5]], color=jade, s=12)
axs[1].set(xscale="log", yscale="log", xlabel=r"Reporter rate scale $R$ per period",
           ylabel="Entropy rate per period", ylim=(0.05, 110),
           title="B   Guaranteed overstatement")
axs[1].legend(frameon=False, fontsize=7.3, loc="upper left")
for ax in axs:
    ax.spines[["top","right"]].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.title.set_size(9.4)
    ax.set_axisbelow(True)
    ax.grid(axis="y", alpha=.16)
fig.savefig(here/"figure.pdf")
