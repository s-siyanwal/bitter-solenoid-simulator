"""Physical constants and material properties (PDF values where given)."""
import math

MU_0 = 4.0 * math.pi * 1e-7          # vacuum permeability [T m/A] (PDF)
C_LIGHT = 299792458.0                 # speed of light [m/s]
GAMMA_PROTON = 2.6752218744e8         # proton gyromagnetic ratio [rad/s/T]

# Copper (PDF: rho = 1.68e-8 Ohm m at 20 C)
RHO_CU_20 = 1.68e-8                   # [Ohm m] (PDF)
ALPHA_CU = 3.93e-3                    # temperature coefficient [1/K] (ASSUMPTION, standard annealed Cu)
K_CU = 390.0                          # thermal conductivity [W/m/K] (ASSUMPTION)
DENS_CU = 8960.0                      # density [kg/m^3] (ASSUMPTION)
CP_CU = 385.0                         # specific heat [J/kg/K] (ASSUMPTION)
YIELD_CU_HARD = 250e6                 # 0.2% yield of hard-drawn Cu [Pa] (ASSUMPTION)

# Water at 20 C (PDF values; used as reference / for tests)
K_WATER_20 = 0.598
PR_WATER_20 = 7.0
NU_KIN_WATER_20 = 1.004e-6
RHO_WATER_20 = 998.2
CP_WATER = 4184.0

# --- Literature-audit P0 additions (AUDIT_REPORT §5) ---
# Contact / joint resistance per plate interface. Default 0 keeps the published
# continuum path bit-stable; HAO Table 2.4 measured ~65 µΩ on a small stack.
R_C_OHM_DEFAULT = 0.0
R_C_HAO_EXAMPLE_OHM = 65e-6

# Stacked axial channels: MON p.49, 104 report friction 10–20× smooth-tube
# values (f ≈ 0.1 "normal"). Apply the multiplier to Darcy f for Δp only;
# keep the conventional (smooth) h / Nu (MON design rule).
FRICTION_MULTIPLIER_SMOOTH = 1.0
FRICTION_MULTIPLIER_STACK_MON = 15.0

# Overlap sector shared with helical.py (ASSUMPTION unless overridden).
OVERLAP_DEG_DEFAULT = 30.0

# Supply voltage ceiling used only to report R_c headroom (PDF: 8 V).
V_SUPPLY_MAX_DEFAULT = 8.0

# Soft Re floor for correlation validity (BIR p.4–5 uses ~5500; Dittus–Boelter
# claimed from 2300). Hard Re ≥ 1e4 floor removed from the optimiser.
RE_MIN_CORRELATION = 5500.0
