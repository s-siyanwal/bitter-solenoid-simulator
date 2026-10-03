# Hosting the static demo

The site is `docs/index.html` (overview) plus the interactive demo, the pair `docs/demo.html` and `docs/bittersim.js`. It has no build step and no server. Open `demo.html` in a browser. Advanced panels are in a collapsed section on the demo page. The Python Monte Carlo (`python -m bittersim emulate`) is the reference. The page runs a capped, seeded subset (at most 300 draws) of the same formulas and says so on the emulation panel.

The owner authorized a public repository so GitHub Pages can serve `docs/` from `main`. On GitHub Free and on Pro or Team the published site is public. Private Pages would need Enterprise Cloud. The interactive page is the demo. Do not treat `results/` figures as a substitute.

## GitHub Free

GitHub Pages builds only from a public repository on a Free account. Publishing this demo from the private repo will fail, or a site that was public will be unpublished, on a Free account.

## GitHub Pro or Team

Pages can build from a private repository, but the site is still public at `https://s-siyanwal.github.io/bitter-solenoid-simulator/` unless the account is Enterprise Cloud with private Pages visibility.

## Recommended path when the source stays private

Use a separate public repository `bitter-solenoid-demo` that contains only `docs/`. A workflow on a tag can copy that directory across. Do not copy `results/` or `notebooks/` if those are private.

## Manual workflow

Pages is configured with **Source: GitHub Actions**. `.github/workflows/pages.yml` uploads `docs/` and deploys it on every push to `main` that touches `docs/` (or on manual dispatch). The old setting served branch `pages-demo` at `/`, so Jekyll rendered the README and `bittersim.js` was missing (404); that was the cause of the broken demo, fixed 2026-10-03.
