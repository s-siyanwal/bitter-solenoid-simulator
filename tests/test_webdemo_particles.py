"""Browser particle panel (docs/bittersim.js) vs bittersim.particles."""
import json, os, shutil, subprocess
import numpy as np
import pytest
from bittersim import particles as Pt

JS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "docs", "bittersim.js")
pytestmark = pytest.mark.skipif(shutil.which("node") is None, reason="node not installed")
G = (0.05000755934750425, 0.29999726655942893, 1.1019511068456953, 230440.62137427254)


def _node(body):
    return json.loads(subprocess.check_output(["node", "-e", body]).decode())


def test_js_P1_matches_python_at_same_N():
    pts = [[0, 0, 0], [0, 0, 0.2], [0.01, 0.005, 0.01], [0, 0, -0.5], [0.012, -0.003, -0.008]]
    for N in (256, 2048):
        js = _node("const b=require(%r);console.log(JSON.stringify(b.particleField(%s,%r,%r,%r,%r,%d)))" % ((JS, json.dumps(pts)) + G + (N,)))
        py = Pt.element_field(np.array(pts, float), *G, N=N, parallel=True)
        assert np.max(np.abs(np.array(js) - py)) < 1e-12 * np.max(np.abs(py))


def test_js_axis_summary_and_bloch_gruneisen():
    out = _node("const b=require(%r);const r=b.particleAxis(%r,%r,%r,%r,1024);"
                "console.log(JSON.stringify({e:r.max_err_inner,c:r.Bz_continuum[10],bg:[b.blochGruneisen(20),b.blochGruneisen(60)]}))" % ((JS,) + G))
    assert out["c"] == pytest.approx(0.5, rel=1e-12)
    from bittersim.fields import B_bitter_axis
    z = np.linspace(-0.5 * G[2], 0.5 * G[2], 21)
    bz = Pt.element_field(np.column_stack([0 * z, 0 * z, z]), *G, N=1024, parallel=True)[:, 2]
    e_py = np.max(np.abs(bz / B_bitter_axis(*(G + (z,))) - 1)[np.abs(z) <= 0.4 * G[2] + 1e-12])
    assert out["e"] == pytest.approx(e_py, rel=1e-8)
    assert out["bg"][0] == pytest.approx(Pt.bloch_gruneisen_rho(20.0), rel=1e-12)
    assert out["bg"][1] == pytest.approx(Pt.bloch_gruneisen_rho(60.0), rel=1e-12)


@pytest.mark.parametrize("seed", [1, 2])
def test_js_P2_statistically_unbiased(seed):
    g = _node("const b=require(%r);console.log(JSON.stringify(b.greenKubo(40,3000,%d)))" % (JS, seed))
    assert abs(g["rho"] / g["rho_target"] - 1) < 4 * g["stderr_rel"]
    assert g["rho_target"] == pytest.approx(Pt.bloch_gruneisen_rho(40.0), rel=1e-12)
