# Coupled Mission Scaffold

This package is the boundary between:

- pyCycle engine models in `propulsion/power_arc/duality.py`
- the sibling HyTank checkout at `../HyTank`
- AeroSandbox aircraft sizing and mission optimization

The first-pass workflow is deck-based:

1. Use pyCycle to generate an engine deck in SI units.
2. Generate or provide a tank deck if the old deck-based tank workflow is still needed.
3. Use `EngineDeck` and tank deck CSVs inside an AeroSandbox mission model.

This avoids calling OpenMDAO inside every AeroSandbox optimizer evaluation.

## AeroSandbox Mission

`coupled_mission/aerosandbox_mission.py` builds an actual `asb.Opti`
problem using the same point-mass dynamics pattern as the AeroSandbox
`04 - 2D Aircraft Dynamics for Mission Performance Analysis.ipynb` tutorial:

- `asb.DynamicsPointMass2DSpeedGamma`
- `dyn.add_gravity_force()`
- `dyn.add_force(..., axes="wind")`
- `dyn.constrain_derivatives(opti, time)`

It includes wing, tank, fuel, trajectory, flight-path angle, angle of attack,
mass-state, and throttle decision variables. The default setup constrains the
five altitude/speed waypoints from
`C:\Users\AlexanderAmos\Downloads\mission_profile_5pt.txt`:

| Point | Altitude | Speed |
| --- | ---: | ---: |
| Takeoff | 0 kft | 180 kt |
| Transonic accel | 40 kft | 688 kt |
| Begin cruise | 95 kft | 2860 kt |
| End cruise | 95 kft | 2860 kt |
| Landing | 0 kft | 145 kt |

Thrust, fuel flow, electric power, Mach, dynamic pressure, and aero forces are
outputs of those states and controls, not prescribed profiles. It couples:

- engine thrust and fuel flow from a pyCycle-derived engine deck
- tank capacity and dry mass from a HyTank/LNGTank-derived tank deck
- mission-node lift/drag and mass-burn constraints

Fuel burn follows the AeroSandbox rocket tutorial pattern:

```python
opti.constrain_derivative(
    derivative=-fuel_flow,
    variable=dyn.mass_props.mass,
    with_respect_to=time,
    method="trapezoidal",
)
```

Run it after installing AeroSandbox:

```powershell
pip install aerosandbox
python coupled_mission\aerosandbox_mission.py
```

With an AeroSandbox aircraft file that defines `airplane`, following the
tutorial notebook's `from cessna152 import airplane` pattern:

```powershell
python coupled_mission\aerosandbox_mission.py
```

The current engine model is a differentiable mode-wise scaling of the CSV deck.
Once a pyCycle grid is generated, replace this with a smoother fitted surrogate.

## Tank Deck

From the pyCycle repo root:

```powershell
python coupled_mission\tank_deck.py
```

`TankCase` currently supports constant fuel extraction, heater power, and
ambient temperature across one transient segment. Extend `run_tank_case()` when
the mission model needs vector-valued profiles.

## Engine Deck

`propulsion/data/example_engine_deck.csv` shows the required columns:

```text
mode,mach,altitude_m,throttle,thrust_N,fuel_flow_kg_s,electric_power_W,inlet_area_m2,nozzle_throat_area_m2
```

The current `EngineDeck` uses nearest-neighbor lookup. Replace this with a
smooth AeroSandbox-compatible surrogate after the pyCycle grid is defined.
