# Step-by-Step Guide to Learning pyCycle

This guide walks you through learning pyCycle by **running scripts in a logical order**: thermodynamics first, then flow/atmosphere, then individual elements, and finally full engine cycles. Run commands from the **repo root** (`pyCycle/`).

**Prerequisites:** Python 3.x, OpenMDAO (see README for version), and `pip install -e .` (or `pip install om-pycycle`) so that `import pycycle` works.

---

## Phase 1: Thermodynamics (foundation)

pyCycle uses either **CEA** (chemical equilibrium) or **TABULAR** (precomputed tables) for gas properties. Start with CEA building blocks.

### 1.1 Species data (NASA polynomials)

See how species properties (H0, S0, Cp0) are computed from temperature using the small CO2/CO/O2 dataset:

```bash
python -m pycycle.thermo.cea.species_data
```

**You’ll see:** Products list and sample H0, S0, Cp0 at 500 K and 1200 K.

### 1.2 Thermo data module (co2_co_o2)

Confirm that the thermo_data module loads and builds a `Properties` object:

```bash
python -m pycycle.thermo.cea.thermo_data.co2_co_o2
```

### 1.3 Chemical equilibrium + properties (SetTotalTP)

Run the CEA chemical equilibrium and property calculation at fixed T, P (single point):

```bash
python -m pycycle.thermo.cea.chem_eq
```

**You’ll see:** Equilibrium `n_moles`, mixture `h`, and `gamma` for the CO2/CO/O2 mixture.

### 1.4 Linear-solve RHS (PropsRHS)

See how the right-hand sides for the T/P linear systems are built (used inside ThermoCalcs):

```bash
python -m pycycle.thermo.cea.props_rhs
```

### 1.5 Property calculations (PropsCalcs)

Run the component that computes h, S, gamma, Cp, Cv, rho, R from equilibrium mole numbers and linear-solve results:

```bash
python -m pycycle.thermo.cea.props_calcs
```

### 1.6 Full Thermo group (total T, P)

Run the top-level Thermo group in `total_TP` mode with CEA:

```bash
python -m pycycle.thermo.thermo
```

**You’ll see:** Flow total enthalpy and gamma from the Thermo group.

### 1.7 Static pressure calculator (PsCalc)

When static pressure is known, see how Mach number, velocity, and area are computed:

```bash
python -m pycycle.thermo.static_ps_calc
```

### 1.8 CEA ThermoAdd (mixing composition)

Run the CEA mixer that adds a reactant (e.g. fuel) to an inflow and outputs mixed composition and mass-averaged enthalpy:

```bash
python -m pycycle.thermo.cea.thermo_add
```

### 1.9 Unit conversion (flow station outputs)

Run the component that passes thermo outputs to flow-station names in English units:

```bash
python -m pycycle.thermo.unit_comps
```

### 1.10 Tabular thermo (optional)

If you want to compare with tabular interpolation (uses precomputed `air_jetA.pkl`):

```bash
python -m pycycle.thermo.tabular.tabular_thermo
python -m pycycle.thermo.tabular.thermo_add
python -m pycycle.thermo.tabular.tab_cea_comparison
```

---

## Phase 2: Atmosphere and flow start

### 2.1 Ambient (standard atmosphere at altitude)

Get Ps, Ts, rhos at an altitude (e.g. 30,000 ft):

```bash
python -m pycycle.elements.ambient
```

### 2.2 US1976 atmosphere

Run the US 1976 standard atmosphere component by itself:

```bash
python -m pycycle.elements.US1976
```

### 2.3 Flight conditions (ambient + flow at MN)

Combines altitude and Mach to set freestream flow (used as inlet inlet conditions):

```bash
python -m pycycle.elements.flight_conditions
```

### 2.4 Flow start (total + static at a given MN)

A flow station defined by P, T, MN, W (uses Thermo in total_TP and static_MN modes):

```bash
python -m pycycle.elements.flow_start
```

---

## Phase 3: Engine elements (one at a time)

Each element script builds a tiny OpenMDAO model: independent variables → element → outputs. Run them to see inputs/outputs and get a feel for the API.

### 3.1 Inlet

```bash
python -m pycycle.elements.inlet
```

### 3.2 Compressor

```bash
python -m pycycle.elements.compressor
```

### 3.3 Combustor

```bash
python -m pycycle.elements.combustor
```

### 3.4 Turbine

```bash
python -m pycycle.elements.turbine
```

### 3.5 Nozzle

```bash
python -m pycycle.elements.nozzle
```

### 3.6 Duct

```bash
python -m pycycle.elements.duct
```

### 3.7 Mixer (two flows → one)

```bash
python -m pycycle.elements.mixer
```

### 3.8 Splitter (one flow → two)

```bash
python -m pycycle.elements.splitter
```

### 3.9 Shaft (power balance)

```bash
python -m pycycle.elements.shaft
```

### 3.10 Bleed out

```bash
python -m pycycle.elements.bleed_out
```

### 3.11 CFD start (Ps, V, area, W → flow station)

```bash
python -m pycycle.elements.cfd_start
```

### 3.12 Turbine cooling (optional, more advanced)

```bash
python -m pycycle.elements.cooling
```

### 3.13 Compressor / turbine maps (optional)

These need map data (e.g. NCP01, LPT2269):

```bash
python -m pycycle.elements.compressor_map
python -m pycycle.elements.turbine_map
```

---

## Phase 4: Simple turbojet cycle

Run a full design-point + off-design turbojet (FlightConditions → Inlet → Compressor → Combustor → Turbine → Nozzle, with shaft and performance):

```bash
python example_cycles/simple_turbojet.py
```

**What to do:** Open `example_cycles/simple_turbojet.py` and read how the cycle is built (`Turbojet` class), how flow is connected (`pyc_connect_flow`), and how design/off-design balances are set up. Then run it and inspect design/off-design outputs (e.g. thrust, SFC).

---

## Phase 5: Other example cycles

Run these in order of increasing complexity.

### 5.1 Afterburning turbojet

```bash
python example_cycles/afterburning_turbojet.py
```

### 5.2 Wet simple turbojet

```bash
python example_cycles/wet_simple_turbojet.py
```

### 5.3 High-bypass turbofan

```bash
python example_cycles/high_bypass_turbofan.py
```
*(Can be slow; see example_cycles/tests for benchmark timeouts.)*

### 5.4 Mixed-flow turbofan

```bash
python example_cycles/mixedflow_turbofan.py
```

### 5.5 Single-spool turboshaft

```bash
python example_cycles/single_spool_turboshaft.py
```

### 5.6 Multi-spool turboshaft

```bash
python example_cycles/multi_spool_turboshaft.py
```

### 5.7 Electric propulsor

```bash
python example_cycles/electric_propulsor.py
```

### 5.8 Wet propulsor

```bash
python example_cycles/wet_propulsor.py
```

---

## Phase 6: Unit tests (optional)

Run tests to confirm everything matches expected values. From repo root:

```bash
# Thermo (CEA)
python -m pytest pycycle/thermo/cea/test/ -v

# Thermo (tabular)
python -m pytest pycycle/thermo/tabular/test/ -v

# Thermo (total/static)
python -m pytest pycycle/thermo/test/ -v

# Elements
python -m pytest pycycle/elements/test/ -v

# Run all example cycles (may skip slow ones)
python -m pytest example_cycles/tests/test_all_examples.py -v
```

Individual element tests, e.g.:

```bash
python -m pycycle.elements.test.test_ambient
python -m pycycle.elements.test.test_inlet
python -m pycycle.elements.test.test_compressor
python -m pycycle.elements.test.test_combustor
python -m pycycle.elements.test.test_turbine
python -m pycycle.elements.test.test_nozzle
```

---

## Quick reference: run order summary

| Step | Script | Purpose |
|------|--------|--------|
| 1 | `python -m pycycle.thermo.cea.species_data` | Species H0, S0, Cp0 |
| 2 | `python -m pycycle.thermo.cea.thermo_data.co2_co_o2` | Load thermo data |
| 3 | `python -m pycycle.thermo.cea.chem_eq` | Equilibrium + props |
| 4 | `python -m pycycle.thermo.cea.props_rhs` | RHS for T/P solves |
| 5 | `python -m pycycle.thermo.cea.props_calcs` | h, S, gamma, Cp, Cv, rho, R |
| 6 | `python -m pycycle.thermo.thermo` | Full Thermo group |
| 7 | `python -m pycycle.thermo.static_ps_calc` | Ps → MN, V, area |
| 8 | `python -m pycycle.elements.ambient` | Ps, Ts at altitude |
| 9 | `python -m pycycle.elements.flow_start` | P, T, MN, W → flow station |
| 10 | `python -m pycycle.elements.inlet` | Inlet element |
| 11 | `python -m pycycle.elements.compressor` | Compressor |
| 12 | `python -m pycycle.elements.combustor` | Combustor |
| 13 | `python -m pycycle.elements.turbine` | Turbine |
| 14 | `python -m pycycle.elements.nozzle` | Nozzle |
| 15 | `python example_cycles/simple_turbojet.py` | Full turbojet |

---

## Tips

- **Repo root:** Always run from the directory that contains `pycycle/` and `example_cycles/` (e.g. `c:\...\pyCycle`).
- **OpenMDAO:** If you see `ModuleNotFoundError: No module named 'openmdao'`, install OpenMDAO and pyCycle (e.g. `pip install openmdao om-pycycle` or `pip install -e .` in the repo).
- **Reading code:** After running a script, open the corresponding `.py` file and read the `if __name__ == "__main__":` block to see how the Problem is built.
- **Paper:** The [pyCycle paper](https://www.mdpi.com/2226-4310/6/8/87/pdf) explains cycle modeling and design/off-design setup in detail.
