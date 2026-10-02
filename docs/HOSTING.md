# Hosting the static demo

The interactive demo is the pair `docs/index.html` and `docs/bittersim.js`. It has no build step and no server. Open `index.html` in a browser. The Python Monte Carlo (`python -m bittersim emulate`) is the reference. The page runs a capped, seeded subset (at most 300 draws) of the same formulas and says so on the emulation panel.

The owner authorized a public repository so GitHub Pages can serve `docs/` from `main`. On GitHub Free and on Pro or Team the published site is public. Private Pages would need Enterprise Cloud. The interactive page is the demo. Do not treat `results/` figures as a substitute.

## GitHub Free

GitHub Pages builds only from a public repository on a Free account. Publishing this demo from the private repo will fail, or a site that was public will be unpublished, on a Free account.

## GitHub Pro or Team

Pages can build from a private repository, but the site is still public at `https://s-siyanwal.github.io/bitter-solenoid-simulator/` unless the account is Enterprise Cloud with private Pages visibility.

## Recommended path when the source stays private

Use a separate public repository `bitter-solenoid-demo` that contains only `docs/`. A workflow on a tag can copy that directory across. Do not copy `results/` or `notebooks/` if those are private.

## Manual workflow

`.github/workflows/pages-demo.yml` runs only on `workflow_dispatch`. If the repository is not public it prints the Free-plan limitation and does not upload the site. If the repository is public it uploads `docs/` as an artifact named `bitter-solenoid-docs`. Turning on Pages in the repository settings is a separate manual step and is not done here.
