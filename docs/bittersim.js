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
    var Qtot = fl.Q * nHoles, loops = bitterLoops(R1, R2, L, C), hom = homogeneity(loops, d.dsv / 2);
    return { d: d, C: C, NI: NI, I: I, nTurns: nT, V: V, P: P, Ppump: fl.dp * Qtot / d.eta_pump, Re: fl.Re, h: fl.h, hg: fl.hg,
      dp: fl.dp, flowLmin: Qtot * 60000, nHoles: nHoles, Tout: d.T_in + dTmix, Thot: Th, Tcu: Tcu, dTw: dTw, dTf: dTf, dTc: dTc,
      ppm: hom[0], B0num: hom[1], Jin: C / (lam * R1) / 1e6, hoop: C * d.B0 / lam / 1e6, lam: lam, loops: loops,
      mass: Math.PI * (R2 * R2 - R1 * R1) * L * lam * DCU, fL: GAMMA_P * d.B0 / (2 * Math.PI) / 1e6 };
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
  var api = { evaluate: evaluate, bitterAxis: bitterAxis, fieldLoops: fieldLoops, swissRoll: swissRoll, snrGain: snrGain,
              ellipke: ellipke, loopField: loopField, MU0: MU0 };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.BitterSim = api;
})(this);
