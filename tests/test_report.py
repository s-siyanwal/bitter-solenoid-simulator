"""The living report must carry the published figures and the changelog rule."""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def test_update_report_check():
    script = os.path.join(ROOT, "examples", "update_report.py")
    subprocess.check_call([sys.executable, script, "--check"], cwd=ROOT)


def test_tex_has_theory_and_pseudocode():
    tex = open(os.path.join(ROOT, "report", "bitter_solenoid_report.tex"), encoding="utf-8").read()
    for needle in ("procedure EMULATE", "procedure OPTIMISE", "asinh", "mesoscopic emulation"):
        assert needle in tex, needle


def test_readme_and_summary_are_rendered_from_results():
    """Every README / DESIGN_SUMMARY number comes from results/*.json via make_docs.py."""
    script = os.path.join(ROOT, "examples", "make_docs.py")
    subprocess.check_call([sys.executable, script, "--check"], cwd=ROOT)
