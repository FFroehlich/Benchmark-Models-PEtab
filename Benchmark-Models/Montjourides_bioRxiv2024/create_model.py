#!/usr/bin/env python3
"""Build the SBML model for the dynRAS model (Montjourides et al. 2024).

Faithfully encodes the default "Jafari et al. 2024" scenario of the dynRAS
recirculating-aquaculture-system model
(https://github.com/Marizauto/dynRAS) as an SBML rate-rule model.

The 28 ODE states span three compartments -- fish tank (FT, 1000 L),
biofilter (B1, 800 L) and degasser (DGS, 700 L) -- and describe the carbonate
system, ammonia/nitrite, nitrifying bacteria (AOB/NOB) and fish growth.

The reference code drives the diurnal CO2 respiration and TAN excretion with a
``floor``-based time-of-day ``interval``. AMICI supports neither ``floor`` nor
a periodic reset without a per-day event, and 140 daily event
reinitialisations make the very stiff carbonate system numerically
intractable. These diurnal drivers are therefore replaced by smooth continuous
raised-cosine surrogates (``co2_diurnal``, ``mod_prop``) with matching diurnal
phase, peak and daily mean. The 14-day biomass (harvest) resets are encoded as
PEtab v2 experiment periods, not here.
"""
import math
import libsbml

MW = "Montjourides_bioRxiv2024"

# ---------------------------------------------------------------------------
# numeric parameter values (params.py, evaluated)
# ---------------------------------------------------------------------------
F = 60 / 60
k1 = 1.49e-2
k_1 = 1.89e4 * 1e-3
k_H = 5.0e10 * 1e-3
kH = 5.0e10 * 10**-9.3
k3 = 1.40e-3 * 1e3
k_3 = (1.40e-3 / (0.45e-14)) * 1e-3
k_4 = 4.3e10 * 1e-3
k4 = 4.3e10 * 2.73e-10
muAOB = (((0.29 + 0.76) / 2)) / (24 * 60 * 60)
muNOB = ((0.28 + 1.04) / 2) / (24 * 60 * 60)
KNH3 = (1 / 14) * 0.01
KAlk = 0.3
KNO2 = 1 / 14
rhoAOB = ((0.05 + 0.15) / 2) / (24 * 60 * 60)
rhoNOB = ((0.05 + 0.15) / 2) / (24 * 60 * 60)
YAOB = 0.21 * 14
YNOB = 0.04 * 14
TGC = 2.7
T = 14.0
w0 = 70.0

V_FT, V_B1, V_DGS = 1000.0, 800.0, 700.0
Bcap = 0.3 * 80000 * 1e3 * 0.05           # 1.2e6
exchange_rate = (V_B1 + V_DGS + V_FT) * 0.25 / (60 * 60 * 24)
r_deg = 0.0015
K_amm = 2.73e-10
OH_dose = 2000.0
HCO3_dose = 2000.0
fish_max = 50000.0
# HCO3 consumption stoichiometry coefficient for AOB growth in the biofilter
C_HCO3 = ((1 / 14) / (7 / 61.1)) / YAOB

# estimated-parameter nominal values (see parameters table)
co2_resp = 62.5      # fish CO2 respiration coefficient
feed_frac = 0.027    # daily feed fraction of biomass

# ---------------------------------------------------------------------------
# initial conditions (Class_definition_and_solver.py)
# ---------------------------------------------------------------------------
init = {
    "CO2aq_FT": 10 / 44.1, "HCO3_FT": 200 / 61.01, "CO32_FT": 1e-6 * 1e3,
    "H_FT": 10**-7.6 * 1e3, "OH_FT": 10**-6.4 * 1e3,
    "NH4_FT": 0.0, "NH3_FT": 0.0, "NO2_FT": 0.0,
    "Fishweight_FT": w0, "Fish_Biomass_FT": 525 * w0,
    "CO2aq_B1": 10 / 44.1, "HCO3_B1": 200 / 61.01, "CO32_B1": 1e-6 * 1e3,
    "H_B1": 10**-7.6 * 1e3, "OH_B1": 10**-6.5 * 1e3,
    "NH4_B1": 0.0, "NH3_B1": 0.0, "NO2_B1": 0.0,
    "AOB_B1": 960000.0, "NOB_B1": 480000.0,
    "CO2aq_DGS": 10 / 44.1, "HCO3_DGS": 200 / 61.01, "CO32_DGS": 1e-6 * 1e3,
    "H_DGS": 10**-7.6 * 1e3, "OH_DGS": 10**-6.4 * 1e3,
    "NH4_DGS": 0.0, "NH3_DGS": 0.0, "NO2_DGS": 0.0,
}
compartment_of = {}
for s in init:
    if s.endswith("_FT"):
        compartment_of[s] = "FT"
    elif s.endswith("_B1"):
        compartment_of[s] = "B1"
    else:
        compartment_of[s] = "DGS"

# Sump (fixed reservoir) concentrations
sump = {
    "CO2aq": 0.0, "HCO3": 70 / 61.01, "CO32": 1e-6 * 1e3,
    "H": 10**-7.5 * 1e3, "OH": 10**-6.5 * 1e3,
    "NH4": 0.0, "NH3": 0.0, "NO2": 0.0,
}

# scalar parameters (name -> value)
scalar_params = dict(
    F=F, k1=k1, k_1=k_1, kH=kH, k_H=k_H, k3=k3, k_3=k_3, k4=k4, k_4=k_4,
    muAOB=muAOB, muNOB=muNOB, KNH3=KNH3, KAlk=KAlk, KNO2=KNO2,
    rhoAOB=rhoAOB, rhoNOB=rhoNOB, YAOB=YAOB, YNOB=YNOB, C_HCO3=C_HCO3,
    TGC=TGC, Temp=T, w0=w0, V_FT=V_FT, V_B1=V_B1, V_DGS=V_DGS, Bcap=Bcap,
    exchange_rate=exchange_rate, r_deg=r_deg, K_amm=K_amm,
    OH_dose=OH_dose, HCO3_dose=HCO3_dose, fish_max=fish_max,
    co2_resp=co2_resp, feed_frac=feed_frac,
    CO2aq_Sump=sump["CO2aq"], HCO3_Sump=sump["HCO3"], CO32_Sump=sump["CO32"],
    H_Sump=sump["H"], OH_Sump=sump["OH"], NH4_Sump=sump["NH4"],
    NH3_Sump=sump["NH3"], NO2_Sump=sump["NO2"],
)

# assignment-rule expressions (helper quantities)
maxdose_HCO3 = "0.0001 * HCO3_dose / V_B1"
maxdose_OH = "0.0001 * OH_dose / V_B1"
T200 = "(200/50.04)"
assignments = {
    "day": "time / 86400",
    # diurnal phase (rad); noon = pi, midnight = 0/2pi.
    "ph": "2 * 3.141592653589793 * time / 86400",
    # Smooth continuous surrogate for the reference's floor-based diurnal drivers
    # (AMICI supports neither floor nor a periodic reset without a per-day event,
    # and 140 daily event reinitialisations make the stiff carbonate system
    # numerically intractable). The raised cosine matches the reference tent's
    # midnight value (1) and noon peak (1.037^12 = 1.548).
    "co2_diurnal": "1 + 0.274 * (1 - cos(ph))",
    "co2_base": ("co2_resp * (Fishweight_FT*1e-3)^(-0.3) * 1.06^14 * (44/32) / 3600"
                 " * Fish_Biomass_FT/1000 / 44 / V_FT"),
    "CO2P": "co2_base * co2_diurnal",
    # TAN diurnal factor: matches the reference beta(3,3) driver's daily mean (1),
    # midnight floor (0.1) and noon peak; total daily N load is preserved.
    "mod_prop": "1 - 0.9 * cos(ph)",
    "TAN_total": "feed_frac * Fish_Biomass_FT * 0.46 * 0.092 / 14 * 1000",
    "TAN_excr": "mod_prop * TAN_total",
    "NH3_frac": "1 / (1 + (H_FT*1e-3)/K_amm)",
    "NH3_excr": "NH3_frac * TAN_excr",
    "NH4_excr": "(1 - NH3_frac) * TAN_excr",
    "dweight": ("3 * TGC * Temp / (86400*1000)"
                " * (w0^(1.0/3.0) + TGC*Temp*time/(86400*1000))^2"),
    # AOB kinetics in the biofilter (mmol/(L.s) style, before /V)
    "aob_kin": "muAOB * (NH3_B1/(KNH3+NH3_B1)) * (HCO3_B1/(KAlk+HCO3_B1)) * AOB_B1",
    "nob_kin": "muNOB * (NO2_B1/(NO2_B1+KNO2)) * NOB_B1",
    # biofilter alkalinity dosing (piecewise in days)
    "add_HCO3": (f"piecewise({maxdose_HCO3}/(1 + (HCO3_B1/{T200})^10),"
                 " (day <= 14) || (day > 56 && day <= 70)"
                 " || (day > 98 && day <= 112) || (day > 140 && day <= 154), 0)"),
    "T_HCO3_OH": ("piecewise((100/50.04), day > 14 && day <= 42,"
                  " (70/50.04), day > 42 && day <= 56,"
                  " (100/50.04), day > 70 && day <= 84,"
                  " (70/50.04), day > 84 && day <= 98,"
                  " (70/50.04), day > 112 && day <= 140, 1e6)"),
    "T_H_OH": "0.79e-6 * 1e3 * CO2aq_B1 / T_HCO3_OH",
    "add_OH": (f"piecewise({maxdose_OH}/(1 + (T_H_OH/H_B1)^10),"
               " (day > 14 && day <= 56) || (day > 70 && day <= 98)"
               " || (day > 112 && day <= 140), 0)"),
}

# ---------------------------------------------------------------------------
# rate rules (dX/dt) -- transcribed from ChemODE_*.py
# ---------------------------------------------------------------------------
def flow(src, x, V):
    """Recirculation inflow term (src - x)*F/V."""
    return f"({src} - {x})*F/{V}"


rates = {}
# --- fish tank (inflow from degasser) ---
rates["CO2aq_FT"] = f"k_1*H_FT*HCO3_FT - k1*CO2aq_FT + CO2P + {flow('CO2aq_DGS','CO2aq_FT','V_FT')}"
rates["HCO3_FT"] = ("k1*CO2aq_FT - k_1*H_FT*HCO3_FT + k_H*H_FT*CO32_FT - kH*HCO3_FT"
                    f" + {flow('HCO3_DGS','HCO3_FT','V_FT')}")
rates["CO32_FT"] = f"kH*HCO3_FT - k_H*H_FT*CO32_FT + {flow('CO32_DGS','CO32_FT','V_FT')}"
rates["H_FT"] = ("k4*NH4_FT - k_4*NH3_FT*H_FT + k1*CO2aq_FT - k_1*H_FT*HCO3_FT"
                 " + kH*HCO3_FT - k_H*H_FT*CO32_FT + k3 - k_3*H_FT*OH_FT"
                 f" + {flow('H_DGS','H_FT','V_FT')}")
rates["OH_FT"] = f"k3 - k_3*H_FT*OH_FT + {flow('OH_DGS','OH_FT','V_FT')}"
rates["NH4_FT"] = ("-k4*NH4_FT + k_4*NH3_FT*H_FT + NH4_excr/86400/V_FT"
                   f" + {flow('NH4_DGS','NH4_FT','V_FT')}")
rates["NH3_FT"] = ("k4*NH4_FT - k_4*NH3_FT*H_FT + NH3_excr/86400/V_FT"
                   f" + {flow('NH3_DGS','NH3_FT','V_FT')}")
rates["NO2_FT"] = flow("NO2_DGS", "NO2_FT", "V_FT")
rates["Fishweight_FT"] = "dweight"
rates["Fish_Biomass_FT"] = "(Fish_Biomass_FT/Fishweight_FT) * dweight"

# --- biofilter (inflow from fish tank, nitrification + dosing) ---
rates["CO2aq_B1"] = f"k_1*H_B1*HCO3_B1 - k1*CO2aq_B1 + {flow('CO2aq_FT','CO2aq_B1','V_B1')}"
rates["HCO3_B1"] = ("k1*CO2aq_B1 - k_1*H_B1*HCO3_B1 + k_H*H_B1*CO32_B1 - kH*HCO3_B1"
                    " - aob_kin*C_HCO3/V_B1 + add_HCO3"
                    f" + {flow('HCO3_FT','HCO3_B1','V_B1')}")
rates["CO32_B1"] = f"kH*HCO3_B1 - k_H*H_B1*CO32_B1 + {flow('CO32_FT','CO32_B1','V_B1')}"
rates["H_B1"] = ("k4*NH4_B1 - k_4*NH3_B1*H_B1 + k1*CO2aq_B1 - k_1*H_B1*HCO3_B1"
                 " + kH*HCO3_B1 - k_H*H_B1*CO32_B1 + k3 - k_3*H_B1*OH_B1"
                 f" + {flow('H_FT','H_B1','V_B1')}")
rates["OH_B1"] = f"k3 - k_3*H_B1*OH_B1 + add_OH + {flow('OH_FT','OH_B1','V_B1')}"
rates["NH4_B1"] = f"-k4*NH4_B1 + k_4*NH3_B1*H_B1 + {flow('NH4_FT','NH4_B1','V_B1')}"
rates["NH3_B1"] = ("k4*NH4_B1 - k_4*NH3_B1*H_B1 - aob_kin/V_B1/YAOB"
                   f" + {flow('NH3_FT','NH3_B1','V_B1')}")
rates["NO2_B1"] = ("aob_kin/V_B1/YAOB - nob_kin/V_B1/YNOB"
                   f" + {flow('NO2_FT','NO2_B1','V_B1')}")
rates["AOB_B1"] = "aob_kin*(Bcap - AOB_B1 - NOB_B1) - rhoAOB*AOB_B1"
rates["NOB_B1"] = "nob_kin*(Bcap - AOB_B1 - NOB_B1) - rhoNOB*NOB_B1"

# --- degasser (inflow from biofilter + sump exchange + CO2 stripping) ---
def dgs_flow(x):
    src_b1 = x.replace("_DGS", "_FT").replace("DGS", "B1") if False else None
    return None


rates["CO2aq_DGS"] = ("k_1*H_DGS*HCO3_DGS - k1*CO2aq_DGS - CO2aq_DGS*r_deg"
                      " + (CO2aq_Sump - CO2aq_DGS)*exchange_rate/V_DGS"
                      f" + {flow('CO2aq_B1','CO2aq_DGS','V_DGS')}")
rates["HCO3_DGS"] = ("k1*CO2aq_DGS - k_1*H_DGS*HCO3_DGS + k_H*H_DGS*CO32_DGS - kH*HCO3_DGS"
                     " + (HCO3_Sump - HCO3_DGS)*exchange_rate/V_DGS"
                     f" + {flow('HCO3_B1','HCO3_DGS','V_DGS')}")
rates["CO32_DGS"] = ("kH*HCO3_DGS - k_H*H_DGS*CO32_DGS"
                     " + (CO32_Sump - CO32_DGS)*exchange_rate/V_DGS"
                     f" + {flow('CO32_B1','CO32_DGS','V_DGS')}")
rates["H_DGS"] = ("k4*NH4_DGS - k_4*NH3_DGS*H_DGS + k1*CO2aq_DGS - k_1*H_DGS*HCO3_DGS"
                  " + kH*HCO3_DGS - k_H*H_DGS*CO32_DGS + k3 - k_3*H_DGS*OH_DGS"
                  " + (H_Sump - H_DGS)*exchange_rate/V_DGS"
                  f" + {flow('H_B1','H_DGS','V_DGS')}")
rates["OH_DGS"] = ("k3 - k_3*H_DGS*OH_DGS + (OH_Sump - OH_DGS)*exchange_rate/V_DGS"
                   f" + {flow('OH_B1','OH_DGS','V_DGS')}")
rates["NH4_DGS"] = ("-k4*NH4_DGS + k_4*NH3_DGS*H_DGS"
                    " + (NH4_Sump - NH4_DGS)*exchange_rate/V_DGS"
                    f" + {flow('NH4_B1','NH4_DGS','V_DGS')}")
rates["NH3_DGS"] = ("k4*NH4_DGS - k_4*NH3_DGS*H_DGS"
                    " + (NH3_Sump - NH3_DGS)*exchange_rate/V_DGS"
                    f" + {flow('NH3_B1','NH3_DGS','V_DGS')}")
rates["NO2_DGS"] = ("(NO2_Sump - NO2_DGS)*exchange_rate/V_DGS"
                    f" + {flow('NO2_B1','NO2_DGS','V_DGS')}")

# ---------------------------------------------------------------------------
# assemble SBML
# ---------------------------------------------------------------------------
def L3(expr):
    ast = libsbml.parseL3Formula(expr)
    if ast is None:
        raise ValueError(f"parse error in: {expr}\n{libsbml.getLastParseL3Error()}")
    return ast


doc = libsbml.SBMLDocument(2, 4)
model = doc.createModel()
model.setId(MW)
model.setName(MW)
model.setMetaId(MW)

for cid in ("FT", "B1", "DGS"):
    c = model.createCompartment()
    c.setId(cid)
    c.setName(cid)
    c.setSize(1.0)
    c.setSpatialDimensions(3)
    c.setConstant(True)

for sid, val in init.items():
    s = model.createSpecies()
    s.setId(sid)
    s.setName(sid)
    s.setCompartment(compartment_of[sid])
    s.setInitialConcentration(float(val))
    s.setBoundaryCondition(False)
    s.setConstant(False)
    s.setHasOnlySubstanceUnits(False)

for pid, val in scalar_params.items():
    p = model.createParameter()
    p.setId(pid)
    p.setName(pid)
    p.setValue(float(val))
    p.setConstant(True)

for pid, expr in assignments.items():
    p = model.createParameter()
    p.setId(pid)
    p.setName(pid)
    p.setConstant(False)
    r = model.createAssignmentRule()
    r.setVariable(pid)
    r.setMath(L3(expr))

for sid, expr in rates.items():
    r = model.createRateRule()
    r.setVariable(sid)
    r.setMath(L3(expr))

# publication reference annotation (is-described-by the preprint DOI)
cv = libsbml.CVTerm()
cv.setQualifierType(libsbml.BIOLOGICAL_QUALIFIER)
cv.setBiologicalQualifierType(libsbml.BQB_IS_DESCRIBED_BY)
cv.addResource("https://doi.org/10.1101/2024.06.28.600787")
model.addCVTerm(cv)

if doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR):
    doc.printErrors()
    raise SystemExit("SBML build errors")

import sys
out = sys.argv[1] if len(sys.argv) > 1 else str(__import__("pathlib").Path(__file__).with_name(f"model_{MW}.xml"))
libsbml.writeSBMLToFile(doc, out)
print("wrote", out, "states:", len(rates), "assignments:", len(assignments))
