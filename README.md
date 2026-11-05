
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
following commands in a terminal. If you already have kariba, skip that part.

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

File for the environment setup is provided: [`environment.yml`](environment.yml)

```shell
micromamba env create -f environment.yml
micromamba activate modbhjet_env
```

If the environment already exists and the YAML has changed, update it with
`micromamba env update -n modbhjet_env -f environment.yml` 

### 3. Compile and install

With the environment active, run this from the BHJet repository:

```shell
python -m pip install --no-build-isolation -v .
```

The installation command will compile the C++ extension and install the Python package. `--no-build-isolation` uses the build dependencies already installed from the YAML.

### 4. Recompile after changes

If you need to make edits for BHJet or Kariba (or are here because you forgot how to run your own code), rerun the installation from the BHJet repository root in the same active environment:

```shell
python -m pip install --no-build-isolation --force-reinstall --no-deps -v .
```
If you change any of the cpp you also need to regenerate the bindings with: `litgen_bindings/generate_bindings.py`



### OpenMP threads (!!)

OpenMP support is required at build time and is included through `llvm-openmp`
in the yaml. Zone calculations run in parallel; existing Python model scripts need no changes. To select a thread count, set
`OMP_NUM_THREADS` before starting Python or Jupyter:

```shell
OMP_NUM_THREADS=4 jupyter lab
```

