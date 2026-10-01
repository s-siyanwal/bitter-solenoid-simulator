"""The static web demo (docs/bittersim.js) must reproduce the Python model."""
import json, os, shutil, subprocess
import pytest
from bittersim import BitterDesign, evaluate_design

DOCS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs")


@pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
@pytest.mark.parametrize("p", [{}, {"R2": 0.3, "L": 1.1, "d_plate": 6e-3, "D_hole": 5.9e-3, "v_flow": 1.7, "pitch_factor": 7.96}])
def test_js_matches_python(p):
    js = "const b=require(%r);const r=b.evaluate(%s);console.log(JSON.stringify({I:r.I,V:r.V,P:r.P,Thot:r.Thot,ppm:r.ppm,dp:r.dp}))" % (
        os.path.join(DOCS, "bittersim.js"), json.dumps(p))
    out = json.loads(subprocess.check_output(["node", "-e", js]).decode())
    r = evaluate_design(BitterDesign(**p))
    for k_js, k_py in (("I", "I_A"), ("V", "V_total_V"), ("P", "P_elec_W"), ("Thot", "T_hot_C"),
                       ("ppm", "homogeneity_ppm"), ("dp", "dp_Pa")):
        assert out[k_js] == pytest.approx(r[k_py], rel=1e-9), k_js
