"""Engineering-practice material catalog.

Callers may only select ids from this allow-list. Free-text materials are
rejected. Factors that are not a specific heat certificate are labeled
catalog approximations in the source note.
"""

CONDUCTORS = {
    "ofhc_cu": {
        "id": "ofhc_cu",
        "name": "OFHC copper",
        "t_min_C": -50.0,
        "t_max_C": 200.0,
        "rho20": 1.68e-8,
        "alpha": 3.93e-3,
        "k": 390.0,
        "density": 8960.0,
        "cp": 385.0,
        "yield_annealed_Pa": 70e6,
        "yield_hard_Pa": 250e6,
        "n_m3": 8.47e28,
        "v_fermi": 1.57e6,
        "note": "Default conductor. rho20 = 1.68e-8 ohm m, alpha = 3.93e-3 /K, "
                "k = 390 W/m/K, density 8960, cp 385. n is about 8.47e28 1/m^3 "
                "(one free electron per atom). Yield stores both tempers: "
                "annealed about 70 MPa, hard about 250 MPa (mid of 200-300 MPa).",
    },
    "cu_ag": {
        "id": "cu_ag",
        "name": "Cu-Ag magnet alloy (catalog approximation)",
        "t_min_C": -50.0,
        "t_max_C": 200.0,
        "rho20": 1.68e-8 * 1.08,
        "alpha": 3.93e-3,
        "k": 370.0,
        "density": 8960.0,
        "cp": 385.0,
        "yield_annealed_Pa": 150e6,
        "yield_hard_Pa": 340e6,
        "n_m3": 8.45e28,
        "v_fermi": 1.57e6,
        "note": "Catalog approximation for a typical 0.03-0.1 percent Ag magnet "
                "copper: resistivity about 8 percent above OFHC, thermal "
                "conductivity slightly lower, higher yield. Not a heat certificate.",
    },
    "cu_zr": {
        "id": "cu_zr",
        "name": "Cu-Zr strength alloy (catalog approximation)",
        "t_min_C": -50.0,
        "t_max_C": 200.0,
        "rho20": 1.68e-8 * 1.15,
        "alpha": 3.90e-3,
        "k": 350.0,
        "density": 8940.0,
        "cp": 385.0,
        "yield_annealed_Pa": 200e6,
        "yield_hard_Pa": 420e6,
        "n_m3": 8.40e28,
        "v_fermi": 1.55e6,
        "note": "Catalog approximation for a strength-oriented Cu-Zr magnet "
                "alloy: higher yield, resistivity about 15 percent above OFHC. "
                "Not a heat certificate.",
    },
    "al_1350": {
        "id": "al_1350",
        "name": "Aluminium 1350",
        "t_min_C": -50.0,
        "t_max_C": 150.0,
        "rho20": 2.82e-8,
        "alpha": 4.03e-3,
        "k": 230.0,
        "density": 2705.0,
        "cp": 900.0,
        "yield_annealed_Pa": 28e6,
        "yield_hard_Pa": 80e6,
        "n_m3": 6.02e28,
        "v_fermi": 2.02e6,
        "note": "Allowed, with a warning. Higher resistivity and lower strength "
                "than copper. Water compatibility differs (oxide film, not the "
                "copper/deionized-water system). Not recommended for this Bitter stack. "
                "n uses one conduction electron per atom, a catalog convention.",
    },
}

COOLANTS = {
    "di_water": {
        "id": "di_water",
        "name": "Deionized water",
        "t_min_C": 1.0,
        "t_max_C": 95.0,
        "dielectric": False,
        "freeze_C": 0.0,
        "note": "Existing water correlations, valid roughly 1-95 C. Real Bitter "
                "water is deionized (order 1 Mohm cm) and preferably deaerated. "
                "Aerated water can mark copper at a few m/s. Deaerated high-purity "
                "water is commonly run near 8-15 m/s in high-field labs. This "
                "0.5 T design is far below that.",
    },
    "water_glycol_30": {
        "id": "water_glycol_30",
        "name": "Water-glycol 30 percent (catalog approximation)",
        "t_min_C": -10.0,
        "t_max_C": 95.0,
        "dielectric": False,
        "freeze_C": -14.0,
        "note": "Catalog approximation: freeze point depressed, viscosity about "
                "2.4 times water, conductivity about 0.80 times water, cp about "
                "0.90 times water, density about 1.04 times water. Electrically "
                "conductive. Not a dielectric. Do not treat it as insulation "
                "inside a Bitter hole.",
    },
    "galden_ht135": {
        "id": "galden_ht135",
        "name": "Galden HT135-class fluorocarbon (catalog approximation)",
        "t_min_C": -20.0,
        "t_max_C": 120.0,
        "dielectric": True,
        "freeze_C": -70.0,
        "note": "Dielectric coolant, h will be lower. Catalog approximation for "
                "a Galden HT135-class perfluoropolyether: higher cost, worse "
                "heat transfer than water. Not a product certificate.",
    },
}

INSULATORS = {
    "polyimide_kapton": {
        "id": "polyimide_kapton",
        "name": "Polyimide (Kapton stand-in)",
        "t_min_C": -50.0,
        "t_max_C": 220.0,
        "t_min_m": 2.5e-5,
        "note": "Stand-in for the 0.25 mm insulator assumption. Continuous rating "
                "about 220 C. Dielectric strength is high. Water uptake is a known "
                "Bitter-plate failure mode.",
    },
    "mica": {
        "id": "mica",
        "name": "Mica",
        "t_min_C": -50.0,
        "t_max_C": 500.0,
        "t_min_m": 1.0e-4,
        "note": "High temperature, historically used between Bitter plates. "
                "Water compatibility is better than organic films if the sheets "
                "stay intact. Minimum thickness in this catalog is 0.1 mm.",
    },
    "ptfe": {
        "id": "ptfe",
        "name": "PTFE",
        "t_min_C": -50.0,
        "t_max_C": 260.0,
        "t_min_m": 1.0e-4,
        "note": "Chemically compatible with water. Cold-flow under the axial "
                "compressive load is the practical limit. Minimum thickness 0.1 mm.",
    },
    "g10": {
        "id": "g10",
        "name": "G-10",
        "t_min_C": -40.0,
        "t_max_C": 130.0,
        "t_min_m": 2.0e-4,
        "note": "Poor wet insulation choice. Glass-epoxy wicks water and the "
                "continuous rating is only about 130 C. Warn before using it "
                "between cooled Bitter plates.",
    },
}

HOUSINGS = {
    "ss304": {
        "id": "ss304",
        "name": "304 stainless bore wall",
        "k": 16.0,
        "density": 8000.0,
        "cp": 500.0,
        "conductive": True,
        "note": "Non-current-carrying housing. Used only for a secondary "
                "external-convection term and a mass estimate. Conductivity does "
                "not create a shorted turn unless conductive bore is enabled.",
    },
    "g10_bore": {
        "id": "g10_bore",
        "name": "G-10 bore wall",
        "k": 0.29,
        "density": 1800.0,
        "cp": 1000.0,
        "conductive": False,
        "note": "Insulating bore wall. External convection only. Wet G-10 is a "
                "poor long-term choice if coolant can reach it.",
    },
    "aluminum_6061": {
        "id": "aluminum_6061",
        "name": "6061 aluminium bore wall",
        "k": 167.0,
        "density": 2700.0,
        "cp": 896.0,
        "conductive": True,
        "note": "Conductive bore tube. Eddy and shorted-turn warning applies "
                "only if conductive bore is enabled. This is not a finite-element "
                "eddy-current solution.",
    },
}

ENVIRONMENT = {
    "T_amb_C": {"default": 20.0, "min": -20.0, "max": 45.0, "note": "Ambient air temperature."},
    "p_site_Pa": {"default": 101325.0, "min": 7.0e4, "max": 2.0e5,
                  "note": "Site pressure. Saturation temperature uses Antoine over about 0.07-0.2 MPa."},
    "aeration": {"choices": ["aerated", "deaerated"], "default": "deaerated",
                 "note": "Dissolved-oxygen flag. Aerated water: erosion caution above about 8 m/s."},
    "temper": {"choices": ["annealed", "hard"], "default": "hard",
               "note": "Conductor temper used for the yield check."},
}


def _require(table, kind, ident):
    if ident not in table:
        raise ValueError("unknown %s %r; catalog ids: %s" % (kind, ident, ", ".join(sorted(table))))
    return table[ident]


def get_conductor(ident):
    return _require(CONDUCTORS, "conductor", ident)


def get_coolant(ident):
    return _require(COOLANTS, "coolant", ident)


def get_insulator(ident):
    return _require(INSULATORS, "insulator", ident)


def get_housing(ident):
    return _require(HOUSINGS, "housing", ident)


def yield_pa(conductor, temper):
    if temper not in ("annealed", "hard"):
        raise ValueError("unknown temper %r; choices: annealed, hard" % temper)
    return conductor["yield_annealed_Pa"] if temper == "annealed" else conductor["yield_hard_Pa"]


def coolant_props(ident, T):
    """Return rho, mu, k, cp, Pr, nu at temperature T [C]."""
    spec = get_coolant(ident)
    from .materials import water_props
    T = float(T)
    if ident == "di_water":
        return water_props(T)
    if ident == "water_glycol_30":
        w = water_props(T)
        rho = w["rho"] * 1.04
        mu = w["mu"] * 2.4
        k = w["k"] * 0.80
        cp = w["cp"] * 0.90
        return {"rho": rho, "mu": mu, "k": k, "cp": cp, "Pr": cp * mu / k, "nu": mu / rho,
                "catalog_approximation": True, "freeze_C": spec["freeze_C"]}
    # galden_ht135 catalog point at 25 C, mild viscosity temperature factor
    rho = 1720.0
    k = 0.065
    cp = 1000.0
    mu = 1.72e-3 * (10.0 ** (80.0 * (1.0 / (T + 273.15) - 1.0 / 298.15)))
    return {"rho": rho, "mu": mu, "k": k, "cp": cp, "Pr": cp * mu / k, "nu": mu / rho,
            "catalog_approximation": True, "freeze_C": spec["freeze_C"]}


def describe():
    """Plain-text allow-list for the CLI."""
    lines = ["Material catalog (allow-list). Unknown ids raise ValueError.", ""]
    lines.append("Conductors:")
    for key in sorted(CONDUCTORS):
        c = CONDUCTORS[key]
        lines.append("  %s: %s" % (key, c["name"]))
        lines.append("    T window %.0f to %.0f C; rho20 %.3g; yield annealed/hard %.0f/%.0f MPa" % (
            c["t_min_C"], c["t_max_C"], c["rho20"],
            c["yield_annealed_Pa"] / 1e6, c["yield_hard_Pa"] / 1e6))
        lines.append("    %s" % c["note"])
    lines.append("Coolants (liquid only in v1):")
    for key in sorted(COOLANTS):
        c = COOLANTS[key]
        lines.append("  %s: %s" % (key, c["name"]))
        lines.append("    T window %.0f to %.0f C. %s" % (c["t_min_C"], c["t_max_C"], c["note"]))
    lines.append("Insulators (plate-to-plate dielectric, not a pressure vessel):")
    for key in sorted(INSULATORS):
        c = INSULATORS[key]
        lines.append("  %s: %s" % (key, c["name"]))
        lines.append("    max T %.0f C, min thickness %.3g m. %s" % (c["t_max_C"], c["t_min_m"], c["note"]))
    lines.append("Housing / bore wall (non-current-carrying unless conductive bore is enabled):")
    for key in sorted(HOUSINGS):
        c = HOUSINGS[key]
        lines.append("  %s: %s. %s" % (key, c["name"], c["note"]))
    lines.append("Environment:")
    for key in sorted(ENVIRONMENT):
        lines.append("  %s: %s" % (key, ENVIRONMENT[key]["note"]))
    lines.append("Warnings that always apply: round graded holes are not a Florida-Bitter plate; "
                 "contact resistance is sampled only in the emulation; Swiss-roll SNR is heuristic; "
                 "this catalog is not a certified magnet design.")
    return "\n".join(lines)
