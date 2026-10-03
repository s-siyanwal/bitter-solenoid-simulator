// bittersim.js - browser port of the analytic core of the Python package
// (same equations E4, E5, E11-E18, E25-E27; exact loop fields via AGM elliptic integrals).
(function (root) {
  "use strict";
  var MU0 = 4e-7 * Math.PI, RHO20 = 1.68e-8, ALPHA = 3.93e-3, KCU = 390, DCU = 8960, CPW = 4184;
  var GAMMA_P = 2.6752218744e8, C0 = 299792458;

  function rhoCu(T) { return RHO20 * (1 + ALPHA * (T - 20)); }
  function water(T) {
    var rho = 1000 * (1 - (T + 288.9414) / (508929.2 * (T + 68.12963)) * Math.pow(T - 3.9863, 2));
    var mu = 2.414e-5 * Math.pow(10, 247.8 / (T + 273.15 - 140));
    var k = 0.5706 + 1.756e-3 * T - 6.46e-6 * T * T;
    return { rho: rho, mu: mu, k: k, cp: CPW, Pr: CPW * mu / k };
  }
  function ellipke(m) {
    var a = 1, b = Math.sqrt(1 - m), s = 0.5 * m, w = 0.5;
    for (var i = 0; i < 60; i++) {
      var c = 0.5 * (a - b), an = 0.5 * (a + b); b = Math.sqrt(a * b); a = an; w *= 2; s += w * c * c;
      if (Math.abs(c) < 1e-17) break;
    }
    var K = Math.PI / (2 * a); return [K, K * (1 - s)];
  }
  function loopField(a, I, r, z) {
    if (r < 1e-12 * a) return [0, MU0 * I * a * a / (2 * Math.pow(a * a + z * z, 1.5))];
    var ap2 = (a + r) * (a + r) + z * z, am2 = (a - r) * (a - r) + z * z, m = 4 * a * r / ap2;
    if (m >= 1) return [0, 0];
    var KE = ellipke(m), pre = MU0 * I / (2 * Math.PI * Math.sqrt(ap2));
    return [pre * z / r * (-KE[0] + (a * a + r * r + z * z) / am2 * KE[1]),
            pre * (KE[0] + (a * a - r * r - z * z) / am2 * KE[1])];
  }
  // Gauss-Legendre nodes by Newton iteration
  function leggauss(n) {
    var x = [], w = [];
    for (var i = 1; i <= n; i++) {
      var t = Math.cos(Math.PI * (i - 0.25) / (n + 0.5)), pp = 0;
      for (var it = 0; it < 100; it++) {
        var p0 = 1, p1 = t;
        for (var k = 2; k <= n; k++) { var p2 = ((2 * k - 1) * t * p1 - (k - 1) * p0) / k; p0 = p1; p1 = p2; }
        pp = n * (t * p1 - p0) / (t * t - 1);
        var dt = p1 / pp; t -= dt; if (Math.abs(dt) < 1e-16) break;
      }
      x.push(t); w.push(2 / ((1 - t * t) * pp * pp));
    }
    return [x, w];
  }
  function gauss(n, lo, hi, panels) {
    var g = leggauss(n), xs = [], ws = [];
    for (var p = 0; p < panels; p++) {
      var a = lo + (hi - lo) * p / panels, b = lo + (hi - lo) * (p + 1) / panels;
      for (var i = 0; i < n; i++) { xs.push(0.5 * (b - a) * g[0][i] + 0.5 * (a + b)); ws.push(0.5 * (b - a) * g[1][i]); }
    }
    return [xs, ws];
  }
  function bitterLoops(R1, R2, L, C) {
    var gr = gauss(12, R1, R2, 1), gz = gauss(16, -L / 2, L / 2, 4), loops = [];
    for (var i = 0; i < gr[0].length; i++)
      for (var j = 0; j < gz[0].length; j++) loops.push([gr[0][i], gz[0][j], C / gr[0][i] * gr[1][i] * gz[1][j]]);
    return loops;
  }
  function fieldLoops(loops, r, z) {
    var br = 0, bz = 0;
    for (var i = 0; i < loops.length; i++) { var f = loopField(loops[i][0], loops[i][2], r, z - loops[i][1]); br += f[0]; bz += f[1]; }
    return [br, bz];
  }
  function homogeneity(loops, rad, n) {
    n = n || 400;
    var vals = [], B0, i, f;
    for (i = 0; i < n; i++) {
      var phi = Math.acos(1 - 2 * (i + 0.5) / n), th = Math.PI * (1 + Math.sqrt(5)) * (i + 0.5);
      var x = rad * Math.cos(th) * Math.sin(phi), y = rad * Math.sin(th) * Math.sin(phi), zz = rad * Math.cos(phi);
      f = fieldLoops(loops, Math.hypot(x, y), zz); vals.push(Math.hypot(f[0], f[1]));
    }
    [[0, rad], [0, -rad], [rad, 0]].forEach(function (p) { var g = fieldLoops(loops, p[0], p[1]); vals.push(Math.hypot(g[0], g[1])); });
    f = fieldLoops(loops, 0, 0); B0 = Math.hypot(f[0], f[1]); vals.push(B0);
    return [(Math.max.apply(null, vals) - Math.min.apply(null, vals)) / B0 * 1e6, B0];
  }
  function bitterAxis(R1, R2, L, C, z) {
    function t(zeta) { return Math.asinh(zeta / R1) - Math.asinh(zeta / R2); }
    return 0.5 * MU0 * C * (t(L / 2 - z) + t(L / 2 + z));
  }
  function friction(Re) { return Re < 2300 ? 64 / Re : Math.pow(0.790 * Math.log(Re) - 1.64, -2); }
  function channel(v, D, L, T, Kminor) {
    var w = water(T), Re = w.rho * v * D / w.mu, Pr = w.Pr, Nu, Nug;
    if (Re < 2300) { Nu = Nug = 4.36; } else {
      var f8 = friction(Re) / 8;
      Nug = f8 * (Re - 1000) * Pr / (1 + 12.7 * Math.sqrt(f8) * (Math.pow(Pr, 2 / 3) - 1));
      Nu = 0.023 * Math.pow(Re, 0.8) * Math.pow(Pr, 0.4);
    }
    var f = friction(Re), A = Math.PI * D * D / 4;
    return { Re: Re, Pr: Pr, Nu: Nu, h: Nu * w.k / D, hg: Nug * w.k / D, f: f,
             dp: (f * L / D + Kminor) * 0.5 * w.rho * v * v, mdot: w.rho * v * A, Q: v * A };
  }
  function annulusDT(q, a, b, k) { return b <= a ? 0 : q / (4 * k) * (2 * b * b * Math.log(b / a) - (b * b - a * a)); }

  function evaluate(p) {
    var d = Object.assign({ R1: 0.05, R2: 0.15, L: 0.8, d_plate: 2e-3, d_ins: 0.25e-3, D_hole: 5e-3, v_flow: 2.5,
      pitch_factor: 2.0, T_in: 20, B0: 0.5, dsv: 0.03, K_minor: 1.5, eta_pump: 0.7 }, p || {});
    var R1 = d.R1, R2 = d.R2, L = d.L, lnr = Math.log(R2 / R1), pitch = d.pitch_factor * d.D_hole;
    var nRows = Math.max(1, Math.floor((R2 - R1) / pitch + 1e-9)), drRow = (R2 - R1) / nRows, nHoles = 0, nPer = [];
    for (var k = 0; k < nRows; k++) { var rk = R1 + (k + 0.5) * drRow; var nk = Math.max(1, Math.floor(2 * Math.PI * rk / pitch + 1e-9)); nPer.push(nk); nHoles += nk; }
    var fh = nHoles * Math.PI * d.D_hole * d.D_hole / 4 / (Math.PI * (R2 * R2 - R1 * R1));
    var lamAx = d.d_plate / (d.d_plate + d.d_ins), lam = lamAx * (1 - fh), nT = L / (d.d_plate + d.d_ins);
    var C = d.B0 / (MU0 * (Math.asinh(L / (2 * R1)) - Math.asinh(L / (2 * R2)))), NI = C * L * lnr, I = NI / nT;
    var Awet = Math.PI * d.D_hole * lamAx, Tcu = d.T_in + 10, Tb = d.T_in + 5, P, fl, mtot, dTmix;
    for (var it = 0; it < 30; it++) {
      P = 2 * Math.PI * rhoCu(Tcu) * C * C * L * lnr / lam;
      fl = channel(d.v_flow, d.D_hole, L, Tb, d.K_minor); mtot = fl.mdot * nHoles; dTmix = P / (mtot * CPW);
      var Tb2 = d.T_in + 0.5 * dTmix, Tc2 = Tb2 + P / (fl.h * Awet * nHoles * L);
      var done = Math.abs(Tc2 - Tcu) < 1e-10 && Math.abs(Tb2 - Tb) < 1e-10; Tcu = Tc2; Tb = Tb2; if (done) break;
    }
    P = 2 * Math.PI * rhoCu(Tcu) * C * C * L * lnr / lam;
    var Rturn = 2 * Math.PI * rhoCu(Tcu) / (d.d_plate * (1 - fh) * lnr), V = I * Rturn * nT;
    var n0 = nPer[0], cell = drRow * 2 * Math.PI * (R1 + 0.5 * drRow) / n0, b = Math.sqrt(cell / Math.PI), a = d.D_hole / 2;
    var Th = Tcu, dTw, dTf, dTc;
    for (it = 0; it < 30; it++) {
      var rh = rhoCu(Th), ql = (2 * Math.PI / n0) * rh * C * C / lam * Math.log((R1 + drRow) / R1);
      dTw = ql * L / (fl.mdot * CPW); dTf = ql / Awet / fl.h; dTc = annulusDT(rh * Math.pow(C / (lam * R1), 2), a, b, KCU);
      var Tn = d.T_in + dTw + dTf + dTc, dn = Math.abs(Tn - Th) < 1e-10; Th = Tn; if (dn) break;
    }
    var Qtot = fl.Q * nHoles, loops, hom;
    if (d._fast) { hom = [NaN, d.B0]; loops = []; }
    else { loops = bitterLoops(R1, R2, L, C); hom = homogeneity(loops, d.dsv / 2); }
    return { d: d, C: C, NI: NI, I: I, nTurns: nT, V: V, P: P, Ppump: fl.dp * Qtot / d.eta_pump, Re: fl.Re, h: fl.h, hg: fl.hg,
      dp: fl.dp, flowLmin: Qtot * 60000, nHoles: nHoles, Tout: d.T_in + dTmix, Thot: Th, Tcu: Tcu, dTw: dTw, dTf: dTf, dTc: dTc,
      ppm: hom[0], B0num: hom[1], Jin: C / (lam * R1) / 1e6, hoop: C * d.B0 / lam / 1e6, lam: lam, loops: loops,
      mass: Math.PI * (R2 * R2 - R1 * R1) * L * lam * DCU, fL: GAMMA_P * d.B0 / (2 * Math.PI) / 1e6,
      Cth: Math.PI * (R2 * R2 - R1 * R1) * L * lam * DCU * 385, Rth: 1 / (fl.h * Awet * nHoles * L),
      P20: 2 * Math.PI * RHO20 * C * C * L * lnr / lam };
  }
  // Swiss roll (E25-E27), defaults identical to the Python SwissRoll class
  function swissRoll(fTune, o) {
    o = Object.assign({ r: 5e-3, a: 12e-3, N: 28, eps_r: 3, t: 10e-6, lossMult: 50 }, o || {});
    var F = Math.PI * o.r * o.r / (o.a * o.a), w0 = 2 * Math.PI * fTune;
    var delta = Math.sqrt(2 * RHO20 / (w0 * MU0)), sig = RHO20 / (delta * (1 - Math.exp(-o.t / delta)));
    var G = o.lossMult * 2 * sig / (MU0 * o.r * (o.N - 1));
    return { F: F, w0: w0, Gamma: G, Q: w0 / G, mu: function (f) {
      var w = 2 * Math.PI * f, re = w * w - w0 * w0, im = w * G, den = re * re + im * im;
      // 1 - F w^2 / (re + j im)
      return [1 - F * w * w * re / den, F * w * w * im / den];
    } };
  }
  function snrGain(mu, depth, rc, kappa) {
    depth = depth || 0.05; rc = rc || 0.03; kappa = kappa === undefined ? 0.05 : kappa;
    var eta0 = Math.pow(rc * rc / (rc * rc + depth * depth), 1.5);
    var eta = mu[0] > 0 ? 1 / (1 + (1 / eta0 - 1) / mu[0]) : 0;
    return (eta / eta0) / Math.sqrt(1 + kappa * Math.abs(mu[1]));
  }
// Appended emulation API. The continuum evaluate() above is unchanged.
  var E_CHARGE = 1.60217662e-19, M_E = 9.1093837e-31, K_B = 1.38064852e-23;
  var LABEL = "mesoscopic emulation, not molecular dynamics";

  function coolantProps(id, T) {
    var ids = ["di_water", "galden_ht135", "water_glycol_30"];
    if (ids.indexOf(id) < 0) throw new Error("unknown coolant '" + id + "'; catalog ids: " + ids.join(", "));
    var w = water(T);
    if (id === "di_water") return { rho: w.rho, mu: w.mu, k: w.k, cp: w.cp, Pr: w.Pr };
    if (id === "water_glycol_30") {
      var rho = w.rho * 1.04, mu = w.mu * 2.4, k = w.k * 0.80, cp = w.cp * 0.90;
      return { rho: rho, mu: mu, k: k, cp: cp, Pr: cp * mu / k };
    }
    var muG = 1.72e-3 * Math.pow(10, 80 * (1 / (T + 273.15) - 1 / 298.15));
    return { rho: 1720, mu: muG, k: 0.065, cp: 1000, Pr: 1000 * muG / 0.065 };
  }

  function catalogIds() {
    return {
      conductor: ["ofhc_cu", "cu_ag", "cu_zr", "al_1350"],
      coolant: ["di_water", "water_glycol_30", "galden_ht135"],
      insulator: ["polyimide_kapton", "mica", "ptfe", "g10"],
      housing: ["ss304", "g10_bore", "aluminum_6061"]
    };
  }

  function requireMaterial(kind, id) {
    var ids = catalogIds()[kind];
    if (!ids || ids.indexOf(id) < 0) {
      throw new Error("unknown " + kind + " '" + id + "'; catalog ids: " + (ids ? ids.join(", ") : "none"));
    }
    return id;
  }

  function conductorSpec(id) {
    requireMaterial("conductor", id);
    var table = {
      ofhc_cu: { rho20: 1.68e-8, alpha: 3.93e-3, k: 390, n: 8.47e28, vF: 1.57e6, yA: 70e6, yH: 250e6 },
      cu_ag: { rho20: 1.68e-8 * 1.08, alpha: 3.93e-3, k: 370, n: 8.45e28, vF: 1.57e6, yA: 150e6, yH: 340e6 },
      cu_zr: { rho20: 1.68e-8 * 1.15, alpha: 3.90e-3, k: 350, n: 8.40e28, vF: 1.55e6, yA: 200e6, yH: 420e6 },
      al_1350: { rho20: 2.82e-8, alpha: 4.03e-3, k: 230, n: 6.02e28, vF: 2.02e6, yA: 28e6, yH: 80e6 }
    };
    return table[id];
  }

  function drude(o) {
    var J = o.J, rho = o.rho, n = o.n, B = o.B, vF = o.v_fermi || 1.57e6;
    var kohler = !!o.kohler, aK = o.kohler_a == null ? 1 : o.kohler_a;
    var tau = M_E / (n * E_CHARGE * E_CHARGE * rho);
    var v_d = J / (n * E_CHARGE);
    var hall = (E_CHARGE * B / M_E) * tau;
    var frac = kohler ? aK * hall * hall : 0;
    var JE = J * J * rho * (1 + frac);
    var q = rho * J * J;
    return {
      tau: tau, v_d: v_d, mean_free_path: vF * tau, omega_c_tau: hall,
      kohler_drho_over_rho: frac, J_dot_E: JE, q_continuum: q,
      label: LABEL, model_grade: kohler ? "approximate" : "established"
    };
  }

  function johnson(o) {
    var TK = o.T_C + 273.15, df = o.f_hi - o.f_lo;
    var Sv = 4 * K_B * TK * o.R;
    return {
      S_v: Sv, V_rms: Math.sqrt(Sv * df),
      note: "Supply-sense Johnson voltage. Coil Johnson noise is not the MRI noise floor.",
      model_grade: "established"
    };
  }

  function tsatC(pPa) {
    var pBar = pPa / 1e5;
    return 1730.63 / (5.1962 - Math.log10(pBar)) - 233.426;
  }

  function mulberry32(seed) {
    var a = seed >>> 0;
    return function () {
      a |= 0; a = a + 0x6D2B79F5 | 0;
      var t = Math.imul(a ^ a >>> 15, 1 | a);
      t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
      return ((t ^ t >>> 14) >>> 0) / 4294967296;
    };
  }

  function gaussRand(rng) {
    var u = Math.max(rng(), 1e-12), v = rng();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  function percentile(arr, p) {
    var a = arr.slice().sort(function (x, y) { return x - y; });
    var idx = (a.length - 1) * p;
    var lo = Math.floor(idx), hi = Math.ceil(idx);
    if (lo === hi) return a[lo];
    return a[lo] * (hi - idx) + a[hi] * (idx - lo);
  }

  function emulateFast(p, opt) {
    opt = opt || {};
    var n = Math.max(1, Math.min(300, opt.realizations || 40));
    var seed = opt.seed;
    if (seed == null) throw new Error("seed is required");
    var rng = mulberry32(seed);
    var cond = conductorSpec(p.conductor || "ofhc_cu");
    requireMaterial("coolant", p.coolant || "di_water");
    requireMaterial("insulator", p.insulator || "polyimide_kapton");
    requireMaterial("housing", p.housing || "ss304");
    var baseIn = Object.assign({}, p, { _fast: true });
    var base = evaluate(baseIn);
    var J = base.Jin * 1e6;
    var rho = cond.rho20 * (1 + cond.alpha * (base.Tcu - 20));
    var dru = drude({ J: J, rho: rho, n: cond.n, B: p.B0 || 0.5, v_fermi: cond.vF, kohler: false });
    var R = base.V / base.I;
    var john = johnson({ T_C: base.Tcu, R: R, f_lo: opt.f_lo == null ? 1 : opt.f_lo, f_hi: opt.f_hi == null ? 1e4 : opt.f_hi });
    var store = { B0: [], V: [], P: [], Thot: [], vd: [], mfp: [], Vj: [] };
    var med = opt.contact_median == null ? 1e-6 : opt.contact_median;
    for (var i = 0; i < n; i++) {
      var q = Object.assign({}, p, { _fast: true });
      q.R1 = Math.max(0.03, p.R1 + gaussRand(rng) * 2e-4);
      q.R2 = Math.max(q.R1 + 0.02, p.R2 + gaussRand(rng) * 2e-4);
      q.L = Math.max(0.2, p.L + gaussRand(rng) * 5e-4);
      q.d_plate = p.d_plate * (1 + (rng() * 2 - 1) * 0.01);
      q.D_hole = p.D_hole * (1 + (rng() * 2 - 1) * 0.01);
      q.T_in = p.T_in + gaussRand(rng) * 0.2;
      q.v_flow = Math.max(0.05, p.v_flow * (1 + (rng() * 2 - 1) * 0.10));
      var lot = 1 + (rng() * 2 - 1) * 0.02;
      var nTurnsGuess = q.L / (q.d_plate + (p.d_ins || 2.5e-4));
      var nIf = Math.max(1, Math.round(nTurnsGuess) - 1);
      var Rc = 0;
      for (var k = 0; k < nIf; k++) {
        var z = gaussRand(rng);
        Rc += Math.exp(Math.log(med) + z);
      }
      var nT = q.L / (q.d_plate + (q.d_ins || 2.5e-4));
      var NI = base.I * nT;
      var lnr = Math.log(q.R2 / q.R1);
      var Cc = NI / (q.L * lnr);
      q.B0 = MU0 * Cc * (Math.asinh(q.L / (2 * q.R1)) - Math.asinh(q.L / (2 * q.R2)));
      var ev = evaluate(q);
      var Rbulk = ev.V / ev.I;
      var contact = (Rbulk + Rc) / Rbulk;
      var Ji = ev.Jin * 1e6;
      var rhoi = cond.rho20 * (1 + cond.alpha * (ev.Tcu - 20)) * lot;
      var di = drude({ J: Ji, rho: rhoi, n: cond.n, B: q.B0, v_fermi: cond.vF, kohler: !!opt.kohler });
      var ji = johnson({ T_C: ev.Tcu, R: (Rbulk + Rc) * lot, f_lo: opt.f_lo == null ? 1 : opt.f_lo, f_hi: opt.f_hi == null ? 1e4 : opt.f_hi });
      store.B0.push(q.B0);
      store.V.push(ev.V * contact * lot);
      store.P.push(ev.P * contact * lot);
      store.Thot.push(q.T_in + (ev.Thot - q.T_in) * contact * lot);
      store.vd.push(di.v_d);
      store.mfp.push(di.mean_free_path);
      store.Vj.push(ji.V_rms);
    }
    function pack(arr, nom) {
      return { p05: percentile(arr, 0.05), p50: percentile(arr, 0.5), p95: percentile(arr, 0.95), nominal: nom };
    }
    return {
      label: LABEL,
      model_grade: { johnson: "established", kohler: "approximate", contact_resistance: "approximate", hooge_1f: "speculative" },
      note: "Browser subset, cap 300. Same formulas as Python. The Python Monte Carlo is the reference. Loop homogeneity is the continuum value; this subset does not resample ppm.",
      n: n, seed: seed,
      drude: dru, johnson: john,
      ppm_nominal: base.ppm,
      percentiles: {
        B0: pack(store.B0, base.B0num),
        V: pack(store.V, base.V),
        P: pack(store.P, base.P),
        T_hot: pack(store.Thot, base.Thot),
        v_d: pack(store.vd, dru.v_d),
        mfp: pack(store.mfp, dru.mean_free_path),
        V_johnson: pack(store.Vj, john.V_rms)
      }
    };
  }

  function presets() {
    return {
      pdf_initial_guess: {
        R1: 0.05, R2: 0.15, L: 0.8, d_plate: 0.002, d_ins: 0.00025, D_hole: 0.005,
        v_flow: 2.5, pitch_factor: 2, T_in: 20, B0: 0.5, dsv: 0.03
      },
      seed1_optimum: {
        R1: 0.05000755934750425, R2: 0.29999726655942893, L: 1.1019511068456953,
        d_plate: 0.00598979191409142, d_ins: 0.00025, D_hole: 0.005915140047654358,
        v_flow: 1.6973400741803273, pitch_factor: 7.961978105433783,
        T_in: 20, B0: 0.5, dsv: 0.03
      }
    };
  }

  var PUBLISHED = {
    optimum_to_85C_s: 4922,
    pdf_guess_to_tsat_s: 479,
    optimum_caption: "Published lumped time to 85 C at the seed-1 optimum. Not a boiling model.",
    guess_caption: "PDF-guess adiabatic path to saturation, about 479 s. Invalid after boiling. Separate from the 4922 s figure."
  };

  function coarseLoops(R1, R2, L, C) {
    var gr = gauss(4, R1, R2, 1), gz = gauss(6, -L / 2, L / 2, 2), loops = [], i, j;
    for (i = 0; i < gr[0].length; i++)
      for (j = 0; j < gz[0].length; j++) loops.push([gr[0][i], gz[0][j], C / gr[0][i] * gr[1][i] * gz[1][j]]);
    return loops;
  }

  function fieldGrid(R1, R2, L, C, nRho, nZ) {
    nRho = nRho || 22;
    nZ = nZ || 16;
    var loops = coarseLoops(R1, R2, L, C);
    var rho = [], z = [], Bz = [], iz, ir, zz, rr, row;
    var rMax = R2 * 1.05;
    for (iz = 0; iz < nZ; iz++) {
      zz = -0.55 * L + (1.1 * L) * iz / (nZ - 1);
      z.push(zz);
      row = [];
      for (ir = 0; ir < nRho; ir++) {
        rr = rMax * ir / (nRho - 1);
        row.push(fieldLoops(loops, Math.max(rr, 0), zz)[1]);
      }
      Bz.push(row);
    }
    for (ir = 0; ir < nRho; ir++) rho.push(rMax * ir / (nRho - 1));
    return { rho: rho, z: z, Bz: Bz, note: "Coarse Gauss-loop map for the figure. Homogeneity ppm still comes from evaluate()." };
  }

  function langevinCloud(o) {
    var tau = o.tau, vd = o.v_d, T = o.T_C, seed = o.seed;
    if (seed == null) throw new Error("seed is required");
    var n = Math.max(8, Math.min(400, o.n || 80));
    var steps = o.steps || 30;
    var rng = mulberry32(seed >>> 0);
    var dt = tau / 20;
    var sig = Math.sqrt(2 * K_B * (T + 273.15) / (M_E * tau) * dt);
    var v = [], i, s, mean = 0, m2 = 0;
    for (i = 0; i < n; i++) v.push(gaussRand(rng) * Math.sqrt(K_B * (T + 273.15) / M_E));
    for (s = 0; s < steps; s++) {
      for (i = 0; i < n; i++) v[i] += -(v[i] / tau) * dt + (vd / tau) * dt + sig * gaussRand(rng);
    }
    for (i = 0; i < n; i++) { mean += v[i]; }
    mean /= n;
    for (i = 0; i < n; i++) m2 += (v[i] - mean) * (v[i] - mean);
    return {
      v: v, mean: mean, v_d: vd, stderr: Math.sqrt(m2 / n) / Math.sqrt(n),
      label: LABEL,
      note: "Seeded 1-D Langevin cloud inside one representative volume. Not a sample of the magnet."
    };
  }

  function adiabaticTrace(r, pSite, nSteps) {
    nSteps = nSteps || 240;
    var alpha = 3.93e-3, Cth = r.mass * 385, T = r.Tcu;
    var scale = r.P / (1 + alpha * (r.Tcu - 20));
    var Tsat = tsatC(pSite || 101325);
    var tEnd = 8000, dt = tEnd / 4000, t = 0, i, P;
    var series = [], t85 = null, tTsat = null;
    for (i = 0; i <= 4000; i++) {
      if (i % Math.max(1, Math.floor(4000 / nSteps)) === 0) series.push([t, T]);
      if (t85 === null && T >= 85) t85 = t;
      if (T >= Tsat) { tTsat = t; series.push([t, T]); break; }
      P = scale * (1 + alpha * (T - 20));
      T += dt * P / Cth;
      t += dt;
    }
    return {
      series: series, t85_s: t85, tsat_s: tTsat, Tsat_C: Tsat,
      published_optimum_to_85C_s: PUBLISHED.optimum_to_85C_s,
      published_pdf_guess_to_tsat_s: PUBLISHED.pdf_guess_to_tsat_s,
      note: "Adiabatic lumped trace for the sparkline. Interpretation stops at Tsat. The published 4922 s optimum figure is not this curve."
    };
  }

  // ------------------------------------------------------------------ fMRI layer (Python: harmonics.py, stability.py, rfsnr.py)
  function bitterLoopsN(R1, R2, L, C, nr, nz, panels) {
    var gr = gauss(nr, R1, R2, 1), gz = gauss(nz, -L / 2, L / 2, panels), loops = [];
    for (var i = 0; i < gr[0].length; i++)
      for (var j = 0; j < gz[0].length; j++) loops.push([gr[0][i], gz[0][j], C / gr[0][i] * gr[1][i] * gz[1][j]]);
    return loops;
  }
  function lstsq(V, f) {            // Householder QR least squares (V: m x n array of rows)
    var m = V.length, n = V[0].length, A = V.map(function (r) { return r.slice(); }), b = f.slice(), j, i, k;
    for (j = 0; j < n; j++) {
      var nrm = 0; for (i = j; i < m; i++) nrm += A[i][j] * A[i][j]; nrm = Math.sqrt(nrm);
      var alpha = A[j][j] > 0 ? -nrm : nrm, v = [];
      for (i = 0; i < m; i++) v.push(i < j ? 0 : A[i][j]); v[j] -= alpha;
      var vv = 0; for (i = j; i < m; i++) vv += v[i] * v[i]; if (vv === 0) continue;
      for (k = j; k < n; k++) { var s = 0; for (i = j; i < m; i++) s += v[i] * A[i][k]; s = 2 * s / vv; for (i = j; i < m; i++) A[i][k] -= s * v[i]; }
      var sb = 0; for (i = j; i < m; i++) sb += v[i] * b[i]; sb = 2 * sb / vv; for (i = j; i < m; i++) b[i] -= sb * v[i];
    }
    var x = new Array(n);
    for (j = n - 1; j >= 0; j--) { var t = b[j]; for (k = j + 1; k < n; k++) t -= A[j][k] * x[k]; x[j] = t / A[j][j]; }
    return x;
  }
  var N_AXIS = 12, J_SHIM = 2.0e6;
  function axisNodes(r0) { var z = []; for (var k = 0; k < 33; k++) z.push(r0 * Math.cos(Math.PI * (k + 0.5) / 33)); return z; }
  function polyFitAxis(z, f, r0) {
    var V = z.map(function (zz) { var row = [], x = zz / r0, p = 1; for (var n = 0; n <= N_AXIS; n++) { row.push(p); p *= x; } return row; });
    return lstsq(V, f);
  }
  function bitterAxisZ(R1, R2, L, C, z) {
    function t(zeta) { return Math.asinh(zeta / R1) - Math.asinh(zeta / R2); }
    return 0.5 * MU0 * C * (t(L / 2 - z) + t(L / 2 + z));
  }
  function loopAxis(a, I, z) { return MU0 * I * a * a / (2 * Math.pow(a * a + z * z, 1.5)); }
  function shims(R1, R2, L, C, dsv, o) {
    o = Object.assign({ gap: 0.01, zf: [0.5, 1.5], nr: 12, nz: 32, panels: 8 }, o || {});
    var r0 = dsv / 2, a = R1 - o.gap, zs = o.zf.map(function (f) { return f * a; }), z = axisNodes(r0);
    var main = polyFitAxis(z, z.map(function (zz) { return bitterAxisZ(R1, R2, L, C, zz); }), r0);
    var cols = zs.map(function (zk) { return polyFitAxis(z, z.map(function (zz) { return loopAxis(a, 1, zz - zk) + loopAxis(a, 1, zz + zk); }), r0); });
    var A = [[cols[0][2], cols[1][2]], [cols[0][4], cols[1][4]]], det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    var I0 = (-main[2] * A[1][1] + main[4] * A[0][1]) / det, I1 = (-main[4] * A[0][0] + main[2] * A[1][0]) / det;
    var coil = bitterLoopsN(R1, R2, L, C, o.nr, o.nz, o.panels);
    var tot = coil.concat([[a, zs[0], I0], [a, -zs[0], I0], [a, zs[1], I1], [a, -zs[1], I1]]);
    var h0 = homogeneity(coil, r0), h1 = homogeneity(tot, r0);
    return { ppm_unshimmed: h0[0], ppm_shimmed: h1[0], shim_radius_m: a, shim_z_m: zs, shim_NI_A: [I0, I1],
             shim_power_W: 2 * RHO20 * J_SHIM * (Math.abs(I0) + Math.abs(I1)) * 2 * Math.PI * a,
             zonal_ppm_unshimmed: main.map(function (c) { return c / main[0] * 1e6; }) };
  }
  var ALPHA_L = 16.5e-6;
  var STAB_SPEC = { mode: "current", ripple_ppm: 1.0, ripple_hz: 300.0, drift_ppm_per_h: 2.0, water_amp_K: 0.1,
                    water_period_s: 300.0, water_drift_K_per_h: 0.5, duration_s: 600.0, target_ppm: 1.0 };
  function stability(res, tau, LR, spec) {
    var sp = Object.assign({}, STAB_SPEC, spec || {}), dt = 0.1, n = Math.floor(sp.duration_s / dt + 0.5) + 1;
    var w = 2 * Math.PI * sp.ripple_hz, fL = GAMMA_P * res.B0 / (2 * Math.PI), rise = res.T_cu_mean_C - res.T_in;
    var a = isFinite(tau) ? Math.exp(-dt / tau) : 1, Tc = res.T_cu_mean_C, T0 = res.T_cu_mean_C;
    var mn = { d: 1e300, e: 1e300, r: 1e300, s: 1e300 }, mx = { d: -1e300, e: -1e300, r: -1e300, s: -1e300 }, series = [];
    for (var i = 0; i < n; i++) {
      var t = i * dt;
      if (i > 0) { var tp = (i - 1) * dt; var Tin = res.T_in + sp.water_amp_K * Math.sin(2 * Math.PI * tp / sp.water_period_s) + sp.water_drift_K_per_h * tp / 3600; Tc = a * Tc + (1 - a) * (Tin + rise); }
      var d = sp.drift_ppm_per_h * t / 3600, e = -ALPHA_L * (Tc - T0) * 1e6;
      var r = sp.mode === "current" ? 0 : ((1 + ALPHA * (T0 - 20)) / (1 + ALPHA * (Tc - 20)) - 1) * 1e6, s = d + e + r;
      [["d", d], ["e", e], ["r", r], ["s", s]].forEach(function (q) { mn[q[0]] = Math.min(mn[q[0]], q[1]); mx[q[0]] = Math.max(mx[q[0]], q[1]); });
      if (i % 10 === 0) series.push([t, s]);
    }
    var ra = sp.mode === "current" ? sp.ripple_ppm : sp.ripple_ppm / Math.sqrt(1 + Math.pow(w * LR, 2));
    var b = { psu_ripple: 2 * ra, psu_drift: mx.d - mn.d, thermal_expansion: mx.e - mn.e, resistance_drift: mx.r - mn.r,
              total: mx.s - mn.s + 2 * ra };
    return { budget_ppm: b, total_Hz: b.total * fL / 1e6, f_larmor_Hz: fL, meets_target: b.total <= sp.target_ppm, series: series, spec: sp };
  }
  // Bessel J1 (rational approximations, |err| < 1e-8), complex helpers
  function besselJ1(x) {
    var ax = Math.abs(x), y, ans1, ans2;
    if (ax < 8) {
      y = x * x;
      ans1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1 + y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
      ans2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y))));
      return ans1 / ans2;
    }
    var z = 8 / ax, xx = ax - 2.356194491; y = z * z;
    ans1 = 1 + y * (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 + y * (-0.240337019e-6))));
    ans2 = 0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
    var r = Math.sqrt(0.636619772 / ax) * (Math.cos(xx) * ans1 - z * Math.sin(xx) * ans2);
    return x < 0 ? -r : r;
  }
  function cx(re, im) { return [re, im || 0]; }
  function cadd(a, b) { return [a[0] + b[0], a[1] + b[1]]; }
  function csub(a, b) { return [a[0] - b[0], a[1] - b[1]]; }
  function cmul(a, b) { return [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]; }
  function cdiv(a, b) { var d = b[0] * b[0] + b[1] * b[1]; return [(a[0] * b[0] + a[1] * b[1]) / d, (a[1] * b[0] - a[0] * b[1]) / d]; }
  function cexp(a) { var e = Math.exp(a[0]); return [e * Math.cos(a[1]), e * Math.sin(a[1])]; }
  function csqrt(a) { var r = Math.hypot(a[0], a[1]), re = Math.sqrt((r + a[0]) / 2), im = Math.sqrt(Math.max(0, (r - a[0]) / 2)); return [re, a[1] < 0 ? -im : im]; }
  function cabs2(a) { return a[0] * a[0] + a[1] * a[1]; }
  var RF_DEFAULTS = { a: 0.03, d_wire: 2e-3, gap: 2e-3, t_slab: 0.05, depth: 0.02, sigma: 0.5, T_coil: 293, T_tissue: 310, T_slab: 293 };
  function rfEvaluate(mu, p, f) {
    var q = Object.assign({}, RF_DEFAULTS, p || {}), a = q.a, g = q.gap, t = q.t_slab, d = q.depth, w = 2 * Math.PI * f;
    var N = 6000, ks = [], B1i = [], Rti = [], Rsi = [], i, j, NS = 41;
    var smu = mu ? csqrt(cx(mu[0], mu[1])) : null;
    for (i = 0; i < N; i++) {
      var u = i / (N - 1), k = (60 / a) * u * u + 1e-9, J = besselJ1(k * a), tau, hz2 = 0;
      if (!mu) { tau = cx(Math.exp(-k * (g + t))); }
      else {
        var qk = cdiv(cx(k), smu), E = cexp(cmul(qk, cx(-t))), inc = Math.exp(-k * g), eta = smu;
        var one = cx(1), BpA = cdiv(cmul(E, csub(eta, one)), cadd(eta, one));
        var A = cdiv(cx(2 * inc), cadd(cadd(one, eta), cmul(cmul(BpA, E), csub(one, eta)))), Bp = cmul(BpA, A);
        tau = cadd(cmul(A, E), Bp);
        var prev = null;
        for (j = 0; j < NS; j++) {
          var s = t * j / (NS - 1);
          var h = cmul(qk, csub(cmul(A, cexp(cmul(qk, cx(-s)))), cmul(Bp, cexp(cmul(qk, cx(s - t))))));
          var v = cabs2(h); if (prev !== null) hz2 += 0.5 * (v + prev) * t / (NS - 1); prev = v;
        }
      }
      ks.push(k);
      B1i.push(cmul(cx(k * J * Math.exp(-k * d)), tau));
      Rti.push(cabs2(tau) * J * J / (k * k));
      Rsi.push(J * J * hz2 / k);
    }
    var B1 = cx(0), Rt = 0, Rs = 0;
    for (i = 1; i < N; i++) {
      var dk = ks[i] - ks[i - 1];
      B1 = cadd(B1, cmul(cx(0.5 * dk), cadd(B1i[i], B1i[i - 1]))); Rt += 0.5 * dk * (Rti[i] + Rti[i - 1]); Rs += 0.5 * dk * (Rsi[i] + Rsi[i - 1]);
    }
    var b1 = MU0 * (a / 2) * Math.hypot(B1[0], B1[1]);
    var Rtis = q.sigma * w * w * MU0 * MU0 * Math.PI * a * a / 4 * Rt;
    var Rsl = mu ? w * MU0 * Math.abs(mu[1]) * 2 * Math.PI * (a / 2) * (a / 2) * Rs : 0;
    var delta = Math.sqrt(2 * RHO20 / (w * MU0)), Rc = 2 * a * RHO20 / (q.d_wire * delta);
    var noise = Math.sqrt(q.T_coil * Rc + q.T_tissue * Rtis + q.T_slab * Rsl);
    return { B1_T_per_A: b1, R_coil: Rc, R_tissue: Rtis, R_slab: Rsl, snr_metric: b1 / noise };
  }
  function rfCompare(B0, detune, lossMult, p) {
    var fL = GAMMA_P * B0 / (2 * Math.PI), mu = swissRoll(fL * detune, { lossMult: lossMult === undefined ? 50 : lossMult }).mu(fL);
    var q = Object.assign({}, RF_DEFAULTS, p || {});
    var s = rfEvaluate(mu, q, fL), air = rfEvaluate(null, q, fL), con = rfEvaluate(null, Object.assign({}, q, { gap: 0, t_slab: 0 }), fL);
    return { mu: mu, slab: s, air: air, contact: con, gain_vs_air: s.snr_metric / air.snr_metric, gain_vs_contact: s.snr_metric / con.snr_metric };
  }
  function lumpedTau(P20, Cth, Rth) { var k = 1 / Rth - ALPHA * P20; return k <= 0 ? Infinity : Cth / k; }

  var api = { shims: shims, zonalFit: polyFitAxis, stability: stability, rfEvaluate: rfEvaluate, rfCompare: rfCompare,
              besselJ1: besselJ1, lumpedTau: lumpedTau, STAB_SPEC: STAB_SPEC, RF_DEFAULTS: RF_DEFAULTS,
              evaluate: evaluate, bitterAxis: bitterAxis, fieldLoops: fieldLoops, swissRoll: swissRoll, snrGain: snrGain,
              ellipke: ellipke, loopField: loopField, MU0: MU0, drude: drude, johnson: johnson,
              coolantProps: coolantProps, catalogIds: catalogIds, requireMaterial: requireMaterial,
              emulateFast: emulateFast, tsatC: tsatC, conductorSpec: conductorSpec, LABEL: LABEL,
              presets: presets, PUBLISHED: PUBLISHED, fieldGrid: fieldGrid, langevinCloud: langevinCloud,
              adiabaticTrace: adiabaticTrace };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.BitterSim = api;
})(this);
