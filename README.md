
## To-do's
- [X] cut-off in kariba + here (Ruby?) 
- [X] single zone particle cutoff shape test/proof 
- [X] shock: temperature jump (Ruby?)
- [X] f_pl scaling (Ruby?)
- [ ] photon time scale example plot (e.g. particle cooling plot) with sa break
- [X] move DEFAULTs all to one file (default_values.hpp)
- [ ] disk + AGN target fields (Ruby?)
- [X] adding agnjet jet dynamics (Ruby?)
- [ ] double check Compton switch (Marc/Ruby?)
- [ ] hadronic interactions (Ruby/Dimitris)
- [ ] integrate MLjet (Dimitris?)
- [X] switch to nanobind + litgen (= remove the pybind file)

## Installation

BHJet builds a C++ extension and requires the Kariba radiation library. Run the
following commands in a terminal. If you already have kariba, start from
its root directory and skip the first two commands.

### 1. Clone BHJet and Kariba

The default build expects another (!!) directory named `kariba`:

```shell
git clone --branch main https://github.com/rubyduncan/modBHJet.git
cd modBHJet
git clone https://github.com/rubyduncan/kariba.git ../kariba
git -C ../kariba checkout 6ca6968adf8568e99b2010212fd1419a79d94d35
```

```shell
CMAKE_ARGS="-Dkariba_SOURCE_DIR=/absolute/path/to/kariba" python -m pip install --no-build-isolation -v .
```

### 2. Create and activate the environment

The root [`environment.yml`](environment.yml) includes the compiler toolchain,
CMake, GSL, OpenMP, Nanobind, and the Python dependencies:

```shell
micromamba env create -f environment.yml
micromamba activate modbhjet_env
```

If the environment already exists and the YAML has changed, update it with
`micromamba env update -n modbhjet_env -f environment.yml` (or use
`conda env update` with the same arguments).

### 3. Compile and install

With the environment active, run this from the BHJet repository root:

```shell
python -m pip install --no-build-isolation -v .
python -c "import bhjet; from bhjet import BHJet, BLJet; print(bhjet.__file__)"
```

The installation compiles the C++ extension and installs the Python package into
the active environment. `--no-build-isolation` uses the build dependencies
already installed from the YAML.

### 4. Recompile after changes

After editing BHJet or Kariba, rerun the installation from the BHJet repository
root in the same active environment:

```shell
python -m pip install --no-build-isolation --force-reinstall --no-deps -v .
```

This rebuilds and reinstalls the package without reinstalling its dependencies.
Restart any Python process or Jupyter kernel that has already imported BHJet so
it loads the rebuilt extension. If you change the public C++ interface, also
regenerate the Python bindings using `litgen_bindings/generate_bindings.py`
before rebuilding; this requires the optional `litgen` package.

### OpenMP threads

OpenMP support is required at build time and is included through `llvm-openmp`
in the YAML environment. Zone calculations run in parallel automatically;
existing Python model scripts need no changes. To select a thread count, set
`OMP_NUM_THREADS` before starting Python or Jupyter:

```shell
OMP_NUM_THREADS=4 python my_model.py
OMP_NUM_THREADS=4 jupyter lab
```

Use `OMP_NUM_THREADS=1` for a single thread. If it is unset, the OpenMP runtime
chooses the thread count. When using multiprocessing, choose the threads per
process to fit the available CPU cores.
