"""The static web demo (docs/bittersim.js) must reproduce the Python model."""
import json, os, shutil, subprocess
import pytest
from bittersim import BitterDesign, evaluate_design
from bittersim.emulation import E_CHARGE, M_ELECTRON, K_BOLTZMANN, drude_rve, johnson_voltage
from bittersim.materials import water_props
from bittersim.catalog import coolant_props

DOCS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs")
JS = os.path.join(DOCS, "bittersim.js")


def _node(body):
    return json.loads(subprocess.check_output(["node", "-e", body]).decode())


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
@pytest.mark.parametrize("p", [{}, {"R2": 0.3, "L": 1.1, "d_plate": 6e-3, "D_hole": 5.9e-3, "v_flow": 1.7, "pitch_factor": 7.96}])
def test_js_matches_python(p):
    js = "const b=require(%r);const r=b.evaluate(%s);console.log(JSON.stringify({I:r.I,V:r.V,P:r.P,Thot:r.Thot,ppm:r.ppm,dp:r.dp}))" % (
        JS, json.dumps(p))
    out = _node(js)
    r = evaluate_design(BitterDesign(**p))
    for k_js, k_py in (("I", "I_A"), ("V", "V_total_V"), ("P", "P_elec_W"), ("Thot", "T_hot_C"),
                       ("ppm", "homogeneity_ppm"), ("dp", "dp_Pa")):
        assert out[k_js] == pytest.approx(r[k_py], rel=1e-9), k_js


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_js_drude_johnson_and_water():
    rho, n, J = 1.68e-8, 8.47e28, 1.0e6
    py = drude_rve(J, rho, 20.0, 0.5, n, 0.002, 1.57e6, kohler=False)
    w = water_props(20.0)
    cat = coolant_props("di_water", 20.0)
    assert cat["rho"] == pytest.approx(w["rho"], rel=1e-12)
    Sv = 4.0 * K_BOLTZMANN * (20.0 + 273.15) * 0.02
    js = (
        "const b=require(%r);"
        "const d=b.drude({J:%r,rho:%r,n:%r,B:0.5,d_plate:0.002,v_fermi:1.57e6,kohler:false});"
        "const j=b.johnson({T_C:20,R:0.02,f_lo:1,f_hi:10000});"
        "const w=b.coolantProps('di_water',20);"
        "console.log(JSON.stringify({tau:d.tau,vd:d.v_d,Sv:j.S_v,rho:w.rho,mu:w.mu,k:w.k}));"
    ) % (JS, J, rho, n)
    out = _node(js)
    assert out["tau"] == pytest.approx(py["tau_s"], rel=1e-9)
    assert out["vd"] == pytest.approx(J / (n * E_CHARGE), rel=1e-9)
    assert out["Sv"] == pytest.approx(Sv, rel=1e-9)
    assert out["rho"] == pytest.approx(w["rho"], rel=1e-9)
    assert out["mu"] == pytest.approx(w["mu"], rel=1e-9)
    assert out["k"] == pytest.approx(w["k"], rel=1e-9)
    assert johnson_voltage(20.0, 0.02, 1.0, 10000.0)["S_v_V2_per_Hz"] == pytest.approx(Sv, rel=1e-12)


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_js_rejects_unknown_material():
    js = (
        "const b=require(%r);let msg='';"
        "try{b.coolantProps('unobtainium',20);}catch(e){msg=e.message;}"
        "console.log(JSON.stringify({msg:msg}));"
    ) % JS
    out = _node(js)
    assert "di_water" in out["msg"]
    assert "unobtainium" in out["msg"]


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_js_presets_match_printed_precision():
    js = (
        "const b=require(%r);"
        "const a=b.evaluate(b.presets().pdf_initial_guess);"
        "const c=b.evaluate(b.presets().seed1_optimum);"
        "console.log(JSON.stringify({pdfV:a.V,pdfP:a.P/1e3,pdfT:a.Thot,pdfPpm:a.ppm,"
        "oV:c.V,oP:c.P/1e3,oT:c.Thot,oPpm:c.ppm,oFlow:c.flowLmin}));"
    ) % JS
    out = _node(js)
    assert round(out["pdfV"] + 1e-12, 3) == 19.537
    assert round(out["pdfP"] + 1e-12, 3) == 17.974
    assert round(out["pdfT"] + 1e-12, 3) == 21.153
    assert round(out["pdfPpm"] + 1e-12, 2) == 154.15
    assert round(out["oV"] + 1e-12, 3) == 4.543
    assert round(out["oP"] + 1e-12, 3) == 11.704
    assert round(out["oT"] + 1e-12, 3) == 25.946
    assert round(out["oPpm"] + 1e-12, 2) == 99.95
    assert round(out["oFlow"] + 1e-12, 1) == 321.8


def test_landing_page_links_demo_and_explains_basics():
    html = open(os.path.join(DOCS, "index.html"), encoding="utf-8").read()
    for needle in ('href="demo.html"', "1/r", "slit", "cooling holes", "Lorentz force"):
        assert needle in html, needle
    assert "bittersim.js" not in html          # the landing page stays static and light


def test_demo_page_labels():
    html = open(os.path.join(DOCS, "demo.html"), encoding="utf-8").read()
    assert "Advanced (work in progress)" in html and '<details class="adv" id="advanced">' in html
    for needle in (
        "Mesoscopic emulation, not molecular dynamics",
        "PDF initial guess",
        "seed-1 optimum",
        "PROJECT_REPORT.md",
        "csch",
        "cfield",
        "chist",
        "clang",
        "cspark",
        "Not Florida-Bitter",
        "Not an MRI noise floor",
    ):
        assert needle in html, needle


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
def test_js_fmri_layer_matches_python():
    from bittersim import harmonics as H, stability as S, rfsnr as Q
    geo = (0.05, 0.3, 1.1, 230440.0)
    res = {"B0": 0.5, "T_in": 20.0, "T_cu_mean_C": 20.9}
    js = ("const b=require(%r);const s=b.shims(0.05,0.3,1.1,230440.0,0.03);"
          "const st=b.stability({B0:0.5,T_in:20,T_cu_mean_C:20.9},54.3,0.657,{mode:'voltage'});"
          "const rf=b.rfCompare(0.5,1.03,50);"
          "console.log(JSON.stringify({p0:s.ppm_unshimmed,p1:s.ppm_shimmed,I:s.shim_NI_A,W:s.shim_power_W,"
          "tot:st.budget_ppm.total,res:st.budget_ppm.resistance_drift,g:rf.gain_vs_air,gc:rf.gain_vs_contact,"
          "rs:rf.slab.R_slab,j1:b.besselJ1(3.7)}))") % JS
    out = _node(js)
    s = H.shimmed_homogeneity(*geo, dsv=0.03)
    assert out["p0"] == pytest.approx(s["ppm_unshimmed"], rel=1e-8)
    assert out["p1"] == pytest.approx(s["ppm_shimmed"], rel=1e-6)
    assert out["I"] == pytest.approx(s["shim_NI_A"], rel=1e-6)
    assert out["W"] == pytest.approx(s["shim_power_W"], rel=1e-6)
    st = S.simulate(res, 54.3, 0.657, {"mode": "voltage"})
    assert out["tot"] == pytest.approx(st["budget_ppm"]["total"], rel=1e-9)
    assert out["res"] == pytest.approx(st["budget_ppm"]["resistance_drift"], rel=1e-9)
    c = Q.compare(B0=0.5, detune=1.03, loss_multiplier=50.0)
    from scipy.special import j1
    assert out["j1"] == pytest.approx(float(j1(3.7)), abs=1e-7)
    assert out["g"] == pytest.approx(c["gain_vs_air"], rel=1e-6)      # J1 approximation ~1e-8
    assert out["gc"] == pytest.approx(c["gain_vs_contact"], rel=1e-6)
    assert out["rs"] == pytest.approx(c["slab"]["R_slab"], rel=1e-6)
