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
