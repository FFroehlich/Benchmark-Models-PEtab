#!/usr/bin/env python3
"""Reproduce the dynRAS validation figures (Montjourides et al. 2024) from the
PEtab v2 problem.

Every trajectory is obtained by simulating the PEtab v2 problem itself with
AMICI, reusing the collection's simulation code
(``src/python/simulate.py``:``create_v2_simulator``, which encodes the 14-day
harvest resets as experiment periods). The experiment is simulated with dense
output timepoints so that the harvest condition is applied at every ~14-day
boundary. The experimental data shown as points are the problem's own daily
measurements (``measurementData_*.tsv``).

Requirements: petab (v2), amici, matplotlib, numpy (a C++ compiler is needed
the first time, to build the AMICI model).
Run from anywhere: ``python make_figures.py`` (figures are written next to it).
"""
import logging
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import amici.sim.sundials as ass
import petab.v2 as petab

HERE = Path(__file__).resolve().parent
PID = "Montjourides_bioRxiv2024"
DAY = 86400.0

# reuse the collection's AMICI PEtab simulation code
sys.path.insert(0, str(HERE.parents[1] / "src" / "python"))
from simulate import create_v2_simulator  # noqa: E402

problem = petab.Problem.from_yaml(str(HERE / f"{PID}.yaml"))
nominal = problem.get_x_nominal_dict()
meas = problem.measurement_df

sim = create_v2_simulator(problem, verbose=logging.WARNING)
model, em, solver = sim.model, sim.exp_man, sim.solver
solver.set_return_data_reporting_mode(ass.RDataReporting.full)
STATES = list(model.get_state_ids())


def state(y, name):
    return y[:, STATES.index(name)]


def trajectory(n=2801):
    """Dense state trajectory over the 140-day experiment."""
    edata = em.create_edata(PID, problem_parameters=nominal)
    edata.set_timepoints(np.linspace(0, 140 * DAY, n))
    rdata = ass.run_simulation(model, solver, edata)
    assert rdata.status == 0, f"simulation failed with status {rdata.status}"
    return np.asarray(rdata.ts), np.asarray(rdata.x)


t, y = trajectory()
td = t / DAY
CO2 = state(y, "CO2aq_FT") * 44
pH = -np.log10(state(y, "H_FT") * 1e-3)
alk = (state(y, "OH_FT") + state(y, "HCO3_FT") + 2 * state(y, "CO32_FT")
       - state(y, "H_FT")) * 50.04
TAN = (state(y, "NH4_FT") + state(y, "NH3_FT")) * 14.01
NO2 = state(y, "NO2_FT") * 14.01
biomass = state(y, "Fish_Biomass_FT") / 1000.0            # kg
weight = state(y, "Fishweight_FT")


def data(oid):
    m = meas[meas.observableId == oid]
    return m.time.values / DAY, m.measurement.values


def harvest_lines(ax):
    for k in range(0, 141, 14):
        ax.axvline(k, color="grey", ls="--", lw=0.5)


def save(name):
    plt.tight_layout()
    plt.savefig(HERE / name, dpi=120, bbox_inches="tight")
    plt.close("all")
    print("wrote", name)


# --- Figure 1: model vs data (CO2, pH, alkalinity in the fish tank) ---------
fig, axs = plt.subplots(3, 1, figsize=(8, 7), sharex=True)
for ax, (traj, dat, ylab, oid) in zip(axs, [
    (CO2, "observable_CO2", r"CO$_2$ (mg L$^{-1}$)", "observable_CO2"),
    (pH, "observable_pH", "pH", "observable_pH"),
    (alk, "observable_alkalinity", r"Alkalinity (mg L$^{-1}$ CaCO$_3$)",
     "observable_alkalinity"),
]):
    harvest_lines(ax)
    dt, dv = data(oid)
    ax.plot(dt, dv, ".", color="0.5", ms=3, label="data (Jafari et al. 2024)")
    ax.plot(td, traj, "-", color="C0", lw=1.3, label="simulation")
    ax.set_ylabel(ylab)
    ax.spines[["top", "right"]].set_visible(False)
axs[0].legend(loc="upper right", fontsize=8, framealpha=0.9)
axs[0].set_title("Montjourides et al. 2024 (dynRAS) - fish-tank validation")
axs[-1].set_xlabel("Time (days)")
axs[-1].set_xlim(0, 140)
save("fig1_validation.png")

# --- Figure 2: first-20-day zoom showing the diurnal oscillation ------------
fig, axs = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
for ax, (traj, oid, ylab) in zip(axs, [
    (CO2, "observable_CO2", r"CO$_2$ (mg L$^{-1}$)"),
    (pH, "observable_pH", "pH"),
]):
    for start in np.arange(0, 20, 1):                     # feeding periods
        ax.axvspan(start, start + 0.5, color="0.85")
    dt, dv = data(oid)
    ax.plot(dt, dv, ".", color="0.4", ms=5, label="data")
    ax.plot(td, traj, "-", color="C0", lw=1.3, label="simulation")
    ax.set_ylabel(ylab)
    ax.set_xlim(0, 20)
    ax.spines[["top", "right"]].set_visible(False)
axs[0].legend(loc="upper right", fontsize=8)
axs[0].set_title("Diurnal dynamics (first 20 days)")
axs[-1].set_xlabel("Time (days)")
save("fig2_diurnal_zoom.png")

# --- Figure 3: harvest cycles and the nitrogen cascade ----------------------
fig, axs = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
harvest_lines(axs[0])
axs[0].plot(td, biomass, "-", color="C2", lw=1.3, label="standing biomass")
axs[0].axhline(50, color="C3", ls=":", lw=1, label="harvest threshold (50 kg)")
axs[0].set_ylabel("Fish biomass (kg)")
axs[0].legend(loc="upper left", fontsize=8)
axs[0].spines[["top", "right"]].set_visible(False)
axs[0].set_title("14-day harvest cycles (encoded as PEtab experiment periods)")
harvest_lines(axs[1])
axs[1].plot(td, TAN, "-", color="C4", lw=1.3, label="TAN (mg L$^{-1}$ N)")
axs[1].plot(td, NO2, "-", color="C1", lw=1.3, label="NO$_2$ (mg L$^{-1}$ N)")
axs[1].set_ylabel("N species (mg L$^{-1}$)")
axs[1].set_xlabel("Time (days)")
axs[1].set_xlim(0, 140)
axs[1].legend(loc="upper right", fontsize=8)
axs[1].spines[["top", "right"]].set_visible(False)
save("fig3_biomass_nitrogen.png")

print("total log-likelihood at nominal parameters:", sim.simulate().llh)
