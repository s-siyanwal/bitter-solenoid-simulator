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
