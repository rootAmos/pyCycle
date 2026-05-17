"""Weight summation, method selection, ATA rows, and aircraft CG states."""

from dataclasses import dataclass
from pathlib import Path
import csv

import aerosandbox as asb
import aerosandbox.tools.units as u

try:
    from . import comp_weights_flops as flops
    from . import comp_weights_gasp as gasp
    from . import comp_weights_raymer as raymer
except ImportError:
    import comp_weights_flops as flops
    import comp_weights_gasp as gasp
    import comp_weights_raymer as raymer


WeightInputs = raymer.WeightInputs


@dataclass(frozen=True)
class WeightItem:
    ata: object
    component: str
    qty: object
    unit_mass_kg: object
    total_mass_kg: object
    x_m: object
    y_m: object
    z_m: object
    notes: str = ""


def ata_by_component():
    ata = raymer.ata_by_component()
    ata.update(
        {
            "anti_ice": 30,
            "recorders_displays": 31,
            "lighting": 33,
            "navigation": 34,
            "oxygen": 35,
            "pneumatics": 36,
            "water_waste": 38,
            "apu": 49,
            "doors": 52,
            "pylons_nacelles": 54,
            "windows": 56,
            "duality": 71,
            "propulsive_motors": 71,
            "generators": 24,
            "turbines": 71,
            "burners": 71,
            "hv_cables": 73,
            "lv_system": 24,
        }
    )
    return ata


def selected_component_methods(inputs, component_methods=None, ata_methods=None):
    ata = ata_by_component()
    component_methods = component_methods or inputs.component_methods or {}
    ata_methods = ata_methods or inputs.ata_methods or {}
    selected = {}
    for component, component_ata in ata.items():
        method = ata_methods.get(str(component_ata), ata_methods.get(component_ata))
        if method:
            selected[component] = method
    selected.update(component_methods)
    return selected


def selected_component_weights_lb(inputs, component_methods=None, ata_methods=None):
    component_methods = selected_component_methods(inputs, component_methods, ata_methods)
    method_weights = {
        "raymer": raymer.RaymerWeights.component_weights_lb(inputs),
        "flops": flops.FlopsWeights.component_weights_lb(inputs),
        "gasp": gasp.GaspWeights.component_weights_lb(inputs),
    }
    components = raymer.RaymerWeights.component_weights_lb(inputs)
    for component, method in component_methods.items():
        key = f"{component}_lb"
        if key in method_weights.get(method, {}):
            components[key] = method_weights[method][key]
    return components


def structure_components():
    return raymer.structure_components()


def system_components():
    return raymer.system_components()


def payload_components():
    return raymer.payload_components()


def oew_components():
    return structure_components() + system_components() + tuple(propulsion_item_names()) + ("tank_dry",)


def mtow_components():
    return oew_components() + payload_components() + ("fuel",)


def propulsion_item_names():
    return ("propulsive_motors", "generators", "turbines", "burners", "hv_cables", "lv_system")


def _component_masses_kg(inputs):
    return {
        key.replace("_lb", "_kg"): value * u.lbm
        for key, value in selected_component_weights_lb(inputs).items()
    }


def _merged_component_locations_m(inputs, component_locations_m=None):
    locations = {key: tuple(value) for key, value in inputs.component_locations_m.items()}
    locations.update(component_locations_m or {})
    return locations


def _item(component, mass_kg, locations, ata=None, qty=1, notes=""):
    x_m, y_m, z_m = locations.get(component, (0.0, 0.0, 0.0))
    return WeightItem(
        ata=ata if ata is not None else ata_by_component().get(component, ""),
        component=component,
        qty=qty,
        unit_mass_kg=mass_kg / qty if qty else mass_kg,
        total_mass_kg=mass_kg,
        x_m=x_m,
        y_m=y_m,
        z_m=z_m,
        notes=notes,
    )


def propulsion_items(inputs, locations):
    items_kg = dict(inputs.propulsion_items_kg)
    if not items_kg and inputs.duality_weight_lb:
        items_kg = {"duality": inputs.duality_weight_lb * u.lbm}
    return [_item(name, mass_kg, locations) for name, mass_kg in items_kg.items()]


def weight_items(inputs, component_locations_m=None, component_methods=None, ata_methods=None, include_ata_placeholders=True):
    locations = _merged_component_locations_m(inputs, component_locations_m)
    components = selected_component_weights_lb(inputs, component_methods, ata_methods)
    names = list(structure_components() + system_components())
    items = [_item(name, components[f"{name}_lb"] * u.lbm, locations) for name in names]
    items += [
        _item("crew", inputs.number_crew * inputs.crew_weight_lb * u.lbm, locations),
        _item("passengers", inputs.number_passengers * inputs.passenger_weight_lb * u.lbm, locations),
        _item("cargo", inputs.cargo_weight_lb * u.lbm, locations),
        _item("fuel", inputs.fuel_weight_lb * u.lbm, locations),
        _item("tank_dry", inputs.tank_dry_weight_lb * u.lbm, locations),
    ]
    items += propulsion_items(inputs, locations)
    if include_ata_placeholders:
        items += ata_placeholder_items(items, locations)
    return items


def target_atas_from_mass_props(path=Path("outputs/weights/mass_props_v0.1.8_Draft_2026-05-17.csv")):
    if not path.exists():
        return []
    with path.open(newline="") as stream:
        return sorted({int(row["ata"]) for row in csv.DictReader(stream) if row.get("ata", "").isdigit()})


def ata_placeholder_items(items, locations):
    covered = {int(item.ata) for item in items if str(item.ata).isdigit()}
    return [
        _item(f"ata_{ata}_allowance", 0.0, locations, ata=ata, notes="No component equation selected yet.")
        for ata in target_atas_from_mass_props()
        if ata not in covered
    ]


def _cg_from_items(items, component_names=None, weight_overrides_kg=None):
    component_names = set(component_names or [item.component for item in items])
    weight_overrides_kg = weight_overrides_kg or {}
    selected = [
        (weight_overrides_kg.get(item.component, item.total_mass_kg), item)
        for item in items
        if item.component in component_names
    ]
    mass_kg = sum(mass_kg for mass_kg, _ in selected)
    if float(mass_kg) == 0.0:
        return {"mass_kg": 0.0, "x_m": 0.0, "y_m": 0.0, "z_m": 0.0}
    return {
        "mass_kg": mass_kg,
        "x_m": sum(mass_kg * item.x_m for mass_kg, item in selected) / mass_kg,
        "y_m": sum(mass_kg * item.y_m for mass_kg, item in selected) / mass_kg,
        "z_m": sum(mass_kg * item.z_m for mass_kg, item in selected) / mass_kg,
    }


def aircraft_cg_states(inputs, component_locations_m=None, component_methods=None, ata_methods=None):
    items = weight_items(inputs, component_locations_m, component_methods, ata_methods, include_ata_placeholders=False)
    oew = _cg_from_items(items, oew_components())
    zero_fuel = _cg_from_items(items, oew_components() + payload_components())
    mtow = _cg_from_items(items, mtow_components())
    mlw_fuel_kg = max(0.0, inputs.landing_design_gross_weight_lb * u.lbm - zero_fuel["mass_kg"])
    mlw = _cg_from_items(items, mtow_components(), {"fuel": mlw_fuel_kg})
    return {"OEW": oew, "ZFW": zero_fuel, "MTOW": mtow, "MLW": mlw}


def _weight_breakdown(inputs, component_methods=None, ata_methods=None):
    components = selected_component_weights_lb(inputs, component_methods, ata_methods)
    structure_lb = sum(components[f"{key}_lb"] for key in structure_components())
    systems_lb = sum(components[f"{key}_lb"] for key in system_components())
    crew_payload_lb = inputs.number_crew * inputs.crew_weight_lb
    passenger_payload_lb = inputs.number_passengers * inputs.passenger_weight_lb
    payload_lb = crew_payload_lb + passenger_payload_lb + inputs.cargo_weight_lb
    duality_lb = sum(inputs.propulsion_items_kg.values()) / u.lbm or inputs.duality_weight_lb
    operating_empty_without_propulsion_lb = structure_lb + systems_lb + inputs.tank_dry_weight_lb
    operating_empty_lb = operating_empty_without_propulsion_lb + duality_lb
    total_lb = operating_empty_lb + payload_lb + inputs.fuel_weight_lb
    return {
        "components_lb": components,
        "components_kg": _component_masses_kg(inputs),
        "structure_lb": structure_lb,
        "systems_lb": systems_lb,
        "crew_payload_lb": crew_payload_lb,
        "passenger_payload_lb": passenger_payload_lb,
        "payload_lb": payload_lb,
        "operating_empty_without_engine_lb": operating_empty_without_propulsion_lb,
        "operating_empty_lb": operating_empty_lb,
        "duality_weight_lb": duality_lb,
        "tank_dry_weight_lb": inputs.tank_dry_weight_lb,
        "fuel_weight_lb": inputs.fuel_weight_lb,
        "total_aircraft_weight_lb": total_lb,
        "total_aircraft_mass_kg": total_lb * u.lbm,
    }


def _mass_properties(inputs, component_locations_m=None, include_raymer_propulsion_accessories=False):
    mass_props = {}
    for item in weight_items(inputs, component_locations_m, include_ata_placeholders=False):
        mass_props[item.component] = asb.mass_properties_from_radius_of_gyration(
            mass=item.total_mass_kg,
            x_cg=item.x_m,
            y_cg=item.y_m,
            z_cg=item.z_m,
        )
    total_mass_props = asb.MassProperties(mass=0)
    for mass_prop in mass_props.values():
        total_mass_props = total_mass_props + mass_prop
    return mass_props, total_mass_props


def _weight_item_csv_row(item, group="item"):
    return {
        "group": group,
        "ata": item.ata,
        "component": item.component,
        "qty": item.qty,
        "unit_mass_kg": float(item.unit_mass_kg),
        "total_mass_kg": float(item.total_mass_kg),
        "x_m": float(item.x_m),
        "y_m": float(item.y_m),
        "z_m": float(item.z_m),
        "mx_kgm": float(item.total_mass_kg * item.x_m),
        "my_kgm": float(item.total_mass_kg * item.y_m),
        "mz_kgm": float(item.total_mass_kg * item.z_m),
        "notes": item.notes,
    }


def write_weight_breakdown_csv(inputs, output_csv, component_locations_m=None, component_methods=None, ata_methods=None):
    breakdown = _weight_breakdown(inputs, component_methods, ata_methods)
    items = weight_items(inputs, component_locations_m, component_methods, ata_methods)
    states = aircraft_cg_states(inputs, component_locations_m, component_methods, ata_methods)
    rows = [_weight_item_csv_row(item) for item in sorted(items, key=lambda item: (int(item.ata or 0), item.component))]
    for component, weight_lb in (
        ("crew_payload", breakdown["crew_payload_lb"]),
        ("passenger_payload", breakdown["passenger_payload_lb"]),
        ("payload_total", breakdown["payload_lb"]),
        ("fuel", breakdown["fuel_weight_lb"]),
        ("tank_dry", breakdown["tank_dry_weight_lb"]),
        ("duality", breakdown["duality_weight_lb"]),
        ("structure_total", breakdown["structure_lb"]),
        ("systems_total", breakdown["systems_lb"]),
        ("operating_empty_without_engine", breakdown["operating_empty_without_engine_lb"]),
        ("operating_empty", breakdown["operating_empty_lb"]),
        ("total_aircraft", breakdown["total_aircraft_weight_lb"]),
    ):
        rows.append({"group": "total", "ata": "", "component": component, "qty": "", "unit_mass_kg": "", "total_mass_kg": float(weight_lb * u.lbm), "x_m": "", "y_m": "", "z_m": "", "mx_kgm": "", "my_kgm": "", "mz_kgm": "", "notes": ""})
    for state, cg in states.items():
        rows.append({"group": "cg_state", "ata": "", "component": state, "qty": "", "unit_mass_kg": "", "total_mass_kg": float(cg["mass_kg"]), "x_m": float(cg["x_m"]), "y_m": float(cg["y_m"]), "z_m": float(cg["z_m"]), "mx_kgm": "", "my_kgm": "", "mz_kgm": "", "notes": ""})
    output_csv = Path(output_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=("group", "ata", "component", "qty", "unit_mass_kg", "total_mass_kg", "x_m", "y_m", "z_m", "mx_kgm", "my_kgm", "mz_kgm", "notes"))
        writer.writeheader()
        writer.writerows(rows)
    return output_csv


def calculate__weight(inputs):
    return _weight_breakdown(inputs)["total_aircraft_weight_lb"]


def weight_breakdown_from_aircraft(aircraft):
    return _weight_breakdown(aircraft.to_weight_inputs())


def write_aircraft_weight_breakdown_csv(aircraft, output_csv):
    return write_weight_breakdown_csv(aircraft.to_weight_inputs(), output_csv)


def main():
    output_csv = Path("outputs/weights/weight_breakdown.csv")
    try:
        from .aircraft import Aircraft
    except ImportError:
        from aircraft import Aircraft

    aircraft = Aircraft.from_json(
        mass_kg=4000.0,
        fuel_mass_kg=1200.0,
        tank_dry_mass_kg=350.0,
        propulsion_mass_kg=450.0,
    )
    inputs = aircraft.to_weight_inputs()
    breakdown = _weight_breakdown(inputs)
    write_weight_breakdown_csv(inputs, output_csv)
    print(" weight breakdown")
    print(f"Total aircraft weight: {breakdown['total_aircraft_weight_lb']:.3f} lb")
    print(f"Total aircraft mass: {breakdown['total_aircraft_mass_kg']:.3f} kg")
    print(f"OEW without engine: {breakdown['operating_empty_without_engine_lb']:.3f} lb")
    print(f"Payload: {breakdown['payload_lb']:.3f} lb")
    print(f"Fuel: {breakdown['fuel_weight_lb']:.3f} lb")
    print(f"Tank dry: {breakdown['tank_dry_weight_lb']:.3f} lb")
    print(f"Duality: {breakdown['duality_weight_lb']:.3f} lb")
    print("Major groups:")
    print(f"  Structure: {breakdown['structure_lb']:.3f} lb")
    print(f"  Systems: {breakdown['systems_lb']:.3f} lb")
    print(f"Weight breakdown CSV: {output_csv}")


if __name__ == "__main__":
    main()
