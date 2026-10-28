# Testing

CI builds [libneo](https://github.com/itpplasma/libneo) from a pinned public
commit, installs GECKO with the `dev` and `step` extras and runs

```bash
python -m pytest python/tests -v -ra
```

on Python 3.12, including the STEP tests, which require CadQuery/OCP. The
system and Python build dependencies are listed in the workflow
`.github/workflows/tests.yml`. A local setup equivalent to CI, after installing
those system dependencies:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install numpy scipy netCDF4 scikit-build-core f90wrap==0.3.0
git clone https://github.com/itpplasma/libneo.git .ci-deps/libneo
git -C .ci-deps/libneo checkout f6090e8150c2f5b1c0afde09baadf79b50b1ffe7
CMAKE_BUILD_PARALLEL_LEVEL=8 FC=gfortran python -m pip install \
  --no-build-isolation .ci-deps/libneo -Ccmake.define.LIBNEO_BUILD_TESTING=OFF
python -m pip install -e '.[dev,step]'
python -m pytest python/tests -v -ra
```

## What is tested

The tests check that offset walls are closed and outward-oriented, the RZ
clearance of walls from a densely sampled last closed flux surface, coil tube
profiles, and detection of separated and overlapping coil pairs. They use the
open ALPES geometry in `python/tests/data/alpes/`, which is attributed to The
ALPES Team and licensed CC BY 4.0, separately from the MIT-licensed code.

## Private data

Tests on private geometry are skipped by default. To run them, set
`GECKO_PRIVATE_WOUT` and `GECKO_PRIVATE_COILS` to a VMEC `wout` file and a
coils file and pass `--private-data` to pytest. The private viewer test also
needs a working off-screen rendering backend.
