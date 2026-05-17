"""Hybrid aircraft object with AeroSandbox geometry and sizing adapters."""

from dataclasses import dataclass, field
import json
from pathlib import Path
import sys

import aerosandbox as asb
import aerosandbox.numpy as np
import aerosandbox.tools.units as u

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

DEFAULT_AIRCRAFT_JSON = Path(__file__).with_name("data") / "aircraft.json"

try:
    from .volume import AircraftVolumeInputs
except ImportError:
    from volume import AircraftVolumeInputs

from weight.comp_weights_raymer import WeightInputs


@dataclass(frozen=True)
class PropulsionSystem:
    """Propulsion assumptions used by weight and volume adapters."""

    mass_kg: object = 0.0
    volume_m3: object = 3.0
    number_engines: object = 2.0
    number_propulsive_motors: object = 4.0
    number_generators: object = 4.0
    number_turbines: object = 2.0
    generator_turbine_engine_face_mach: object = 0.35
    generator_turbine_subsonic_pressure_recovery: object = 0.98
    generator_turbine_min_pressure_recovery: object = 0.35
    generator_turbine_supersonic_recovery_coefficient: object = 0.075
    generator_turbine_supersonic_recovery_exponent: object = 1.35
    total_engine_thrust_lb: object = 12000.0
    engine_diameter_ft: object = 2.0
    engine_front_to_cockpit_length_m: object = 10.5
    electric_propulsor_mach_limit: object = 1.2
    propulsive_efficiency: object = 0.75
    motor_controller_power_density_W_kg: object = 8000.0
    generator_power_density_W_kg: object = 8000.0
    turbine_power_density_W_kg: object = 6000.0
    generator_efficiency: object = 0.96
    turbine_mechanical_efficiency: object = 0.98
    mass_fractions: object = field(default_factory=dict)


@dataclass(frozen=True)
class FuelSystem:
    """Fuel and tank assumptions in SI units."""

    mass_kg: object = 0.0
    volume_m3: object = 0.0
    fuel_density_kg_m3: object = 422.0
    tank_dry_mass_kg: object = 0.0


@dataclass(frozen=True)
class Payload:
    """Payload assumptions used by the volume adapter."""

    mass_kg: object = 0.0
    volume_m3: object = 5.0


@dataclass(frozen=True)
class Systems:
    uninstalled_avionics_weight_lb: object = 800.0
    electrical_rating_kva: object = 40.0
    electrical_routing_distance_ft: object = 1.0
    number_mechanical_functions: object = 1.0
    number_hydraulic_utility_functions: object = 5.0
    number_control_functions: object = 4.0


@dataclass(frozen=True)
class Interiors:
    number_crew: object = 1.0
    number_passengers: object = 5.0
    crew_weight_lb: object = 200.0
    passenger_weight_lb: object = 200.0
    cargo_weight_lb: object = 0.0
    furnishings_weight_lb: object = 217.6


@dataclass(frozen=True)
class Aircraft:
    """Top-level sizing object with native AeroSandbox geometry."""

    airplane: object
    mass_kg: object
    propulsion: PropulsionSystem
    fuel: FuelSystem
    payload: Payload
    landing_mass_kg: object
    systems: Systems = field(default_factory=Systems)
    interiors: Interiors = field(default_factory=Interiors)
    ata_methods: object = field(default_factory=dict)
    component_methods: object = field(default_factory=dict)
    component_locations_m: object = field(default_factory=dict)

    @property
    def geometry(self):
        return airplane_geometry(self.airplane)

    def to_weight_inputs(self):
        """Adapt aircraft-level assumptions to Raymer correlation inputs."""
        g = self.geometry
        p = self.propulsion
        f = self.fuel
        s = self.systems
        i = self.interiors

        return WeightInputs(
            design_gross_weight_lb=self.mass_kg / u.lbm,
            landing_design_gross_weight_lb=self.landing_mass_kg / u.lbm,
            ultimate_load_factor=7.5,
            landing_ultimate_load_factor=4.5,
            mach=5.0,
            dynamic_pressure_lb_ft2=500.0,
            wing_area_ft2=g["planform_area_m2"] / u.foot**2,
            aspect_ratio=g["aspect_ratio"],
            taper_ratio=g["taper_ratio"],
            sweep_25_rad=g["sweep_25_rad"],
            root_thickness_to_chord=g["root_thickness_to_chord"],
            wing_mounted_control_area_ft2=g["wing_mounted_control_area_m2"]
            / u.foot**2,
            horizontal_tail_area_ft2=g["horizontal_tail_area_m2"] / u.foot**2,
            horizontal_tail_span_ft=g["horizontal_tail_span_m"] / u.foot,
            fuselage_width_at_htail_ft=g["fuselage_width_m"] / u.foot,
            vertical_tail_area_ft2=g["vertical_tail_area_m2"] / u.foot**2,
            vertical_tail_aspect_ratio=g["vertical_tail_aspect_ratio"],
            vertical_tail_height_ft=g["vertical_tail_height_m"] / u.foot,
            horizontal_tail_height_ft=0.0,
            tail_length_ft=g["tail_length_m"] / u.foot,
            rudder_area_ft2=g["rudder_area_m2"] / u.foot**2,
            fuselage_structural_length_ft=g["fuselage_length_m"] / u.foot,
            fuselage_structural_depth_ft=g["fuselage_height_m"] / u.foot,
            fuselage_structural_width_ft=g["fuselage_width_m"] / u.foot,
            main_gear_length_in=self.airplane.main_gear_length_in,
            nose_gear_length_in=self.airplane.nose_gear_length_in,
            number_engines=p.number_engines,
            total_engine_thrust_lb=p.total_engine_thrust_lb,
            thrust_per_engine_lb=p.total_engine_thrust_lb / p.number_engines,
            engine_diameter_ft=p.engine_diameter_ft,
            engine_front_to_cockpit_length_ft=p.engine_front_to_cockpit_length_m
            / u.foot,
            total_fuel_volume_gal=f.volume_m3 / u.gallon,
            number_mechanical_functions=s.number_mechanical_functions,
            number_hydraulic_utility_functions=s.number_hydraulic_utility_functions,
            number_control_functions=s.number_control_functions,
            electrical_rating_kva=s.electrical_rating_kva,
            electrical_routing_distance_ft=s.electrical_routing_distance_ft,
            number_generators=p.number_generators,
            uninstalled_avionics_weight_lb=s.uninstalled_avionics_weight_lb,
            number_crew=i.number_crew,
            number_passengers=i.number_passengers,
            crew_weight_lb=i.crew_weight_lb,
            passenger_weight_lb=i.passenger_weight_lb,
            cargo_weight_lb=i.cargo_weight_lb,
            furnishings_weight_lb=i.furnishings_weight_lb,
            fuel_weight_lb=f.mass_kg / u.lbm,
            tank_dry_weight_lb=f.tank_dry_mass_kg / u.lbm,
            duality_weight_lb=p.mass_kg / u.lbm,
            propulsion_items_kg={
                name: p.mass_kg * fraction
                for name, fraction in p.mass_fractions.items()
            },
            ata_methods=self.ata_methods,
            component_methods=self.component_methods,
            component_locations_m=self.component_locations_m,
        )

    def to_volume_inputs(self):
        """Adapt aircraft-level assumptions to aircraft volume inputs."""
        return AircraftVolumeInputs(
            planform_area_m2=self.geometry["planform_area_m2"],
            fuel_mass_kg=self.fuel.mass_kg,
            propulsion_volume_m3=self.propulsion.volume_m3,
            payload_volume_m3=self.payload.volume_m3,
            fuel_1_density_kg_m3=self.fuel.fuel_density_kg_m3,
        )

    @classmethod
    def from_json(
        cls,
        path=DEFAULT_AIRCRAFT_JSON,
        *,
        mass_kg=0.0,
        fuel_mass_kg=0.0,
        tank_dry_mass_kg=0.0,
        propulsion_mass_kg=0.0,
    ):
        data = load_aircraft_definition(path)
        geometry = data["geometry"]
        fuel = dict(data["fuel"])
        payload = data["payload"]
        propulsion = dict(data["propulsion"])
        systems = data.get("systems", {})
        interiors = data.get("interiors", {})
        aircraft = data.get("aircraft", {})
        ata_methods = data.get("weight_methods", {}).get("ata_methods", {})
        component_methods = data.get("weight_methods", {}).get("component_methods", {})
        component_locations_m = data.get("mass_properties", {}).get("component_locations_m", {})
        landing_mass_fraction = aircraft.get("landing_mass_fraction", 0.85)
        fuel_density_kg_m3 = fuel["fuel_density_kg_m3"]
        propulsion["mass_kg"] = propulsion_mass_kg
        airplane = build_geometric_asb_airplane(**geometry)
        return cls(
            airplane=airplane,
            mass_kg=mass_kg,
            landing_mass_kg=mass_kg * landing_mass_fraction,
            propulsion=PropulsionSystem(**propulsion),
            fuel=FuelSystem(
                mass_kg=fuel_mass_kg,
                volume_m3=fuel_mass_kg / fuel_density_kg_m3,
                fuel_density_kg_m3=fuel_density_kg_m3,
                tank_dry_mass_kg=tank_dry_mass_kg,
            ),
            payload=Payload(**payload),
            systems=Systems(**systems),
            interiors=Interiors(**interiors),
            ata_methods=ata_methods,
            component_methods=component_methods,
            component_locations_m=component_locations_m,
        )


def _flatten_nested_geometry(geometry):
    """Convert aircraft.json geometry subsections into build_airplane kwargs."""
    if any(key in geometry for key in ("wing", "fuselage", "vtail", "landing_gear", "controls")):
        wing = geometry.get("wing", {})
        fuselage = geometry.get("fuselage", {})
        vtail = geometry.get("vtail", {})
        landing_gear = geometry.get("landing_gear", {})
        controls = geometry.get("controls", {})
        geometry = {
            "name": geometry.get("name", "Sizing Aircraft"),
            "planform_area_m2": wing.get("planform_area_m2", 80.0),
            "aspect_ratio": wing.get("aspect_ratio", 3.0),
            "taper_ratio": wing.get("taper_ratio", 0.25),
            "sweep_25_deg": wing.get("sweep_25_deg", 60.0),
            "root_thickness_to_chord": wing.get("root_thickness_to_chord", 0.06),
            "main_wing_tip_le_x_m": wing.get("tip_le_x_m", 2.066),
            "fuselage_length_m": fuselage.get("length_m", 30.0),
            "fuselage_height_m": fuselage.get("height_m", 3.6),
            "fuselage_width_m": fuselage.get("width_m", 3.0),
            "vtail_area_m2": vtail.get("area_m2", 20.8),
            "vtail_span_m": vtail.get("span_m", 5.12),
            "vtail_dihedral_angle_deg": vtail.get("dihedral_angle_deg", 37.0),
            "vtail_le_x_m": vtail.get("le_x_m", 16.5),
            "tail_length_m": vtail.get("tail_length_m", 16.5),
            "main_gear_length_in": landing_gear.get("main_gear_length_in", 42.0),
            "nose_gear_length_in": landing_gear.get("nose_gear_length_in", 30.0),
            "rudder_area_m2": controls.get("rudder_area_m2", 2.0),
            "wing_mounted_control_area_m2": controls.get("wing_mounted_control_area_m2", 6.4),
        }
    return dict(geometry)


def load_aircraft_definition(path=DEFAULT_AIRCRAFT_JSON):
    with Path(path).open("r", encoding="utf-8") as stream:
        data = json.load(stream)
    geometry = _flatten_nested_geometry(dict(data.get("geometry", {})))
    if "sweep_25_deg" in geometry and "sweep_25_rad" not in geometry:
        geometry["sweep_25_rad"] = np.radians(geometry.pop("sweep_25_deg"))
    data["geometry"] = geometry
    return data


def propulsion_from_definition(path=DEFAULT_AIRCRAFT_JSON, **overrides):
    data = dict(load_aircraft_definition(path)["propulsion"])
    data.update({key: value for key, value in overrides.items() if value is not None})
    return PropulsionSystem(**data)


def fuel_from_definition(path=DEFAULT_AIRCRAFT_JSON, **overrides):
    data = dict(load_aircraft_definition(path)["fuel"])
    data.update({key: value for key, value in overrides.items() if value is not None})
    mass_kg = data.get("mass_kg", 0.0)
    tank_dry_mass_kg = data.get("tank_dry_mass_kg", 0.0)
    return FuelSystem(
        mass_kg=mass_kg,
        volume_m3=mass_kg / data["fuel_density_kg_m3"],
        fuel_density_kg_m3=data["fuel_density_kg_m3"],
        tank_dry_mass_kg=tank_dry_mass_kg,
    )


def payload_from_definition(path=DEFAULT_AIRCRAFT_JSON, **overrides):
    data = dict(load_aircraft_definition(path)["payload"])
    data.update({key: value for key, value in overrides.items() if value is not None})
    return Payload(**data)


def systems_from_definition(path=DEFAULT_AIRCRAFT_JSON, **overrides):
    data = dict(load_aircraft_definition(path).get("systems", {}))
    data.update({key: value for key, value in overrides.items() if value is not None})
    return Systems(**data)


def interiors_from_definition(path=DEFAULT_AIRCRAFT_JSON, **overrides):
    data = dict(load_aircraft_definition(path).get("interiors", {}))
    data.update({key: value for key, value in overrides.items() if value is not None})
    return Interiors(**data)


def airplane_geometry(airplane):
    """Return derived adapter geometry from an AeroSandbox airplane."""
    main_wing = next(wing for wing in airplane.wings if wing.name == "Main Wing")
    vtail = next(wing for wing in airplane.wings if wing.name == "VTail")
    fuselage = airplane.fuselages[0]

    planform_area_m2 = main_wing.area()
    span_m = main_wing.span()
    vtail_area_m2 = vtail.area()
    vtail_span_m = vtail.span()
    vtail_dihedral_rad = np.radians(airplane.vtail_dihedral_angle_deg)
    horizontal_tail_area_m2 = vtail_area_m2 * np.cos(vtail_dihedral_rad) ** 2
    vertical_tail_area_m2 = vtail_area_m2 * np.sin(vtail_dihedral_rad) ** 2
    horizontal_tail_span_m = vtail_span_m * np.cos(vtail_dihedral_rad)
    vertical_tail_height_m = 0.5 * vtail_span_m * np.sin(vtail_dihedral_rad)
    fuselage_length_m = fuselage.xsecs[-1].xyz_c[0] - fuselage.xsecs[0].xyz_c[0]
    fuselage_height_m = fuselage.xsecs[0].height
    fuselage_width_m = fuselage.xsecs[0].width
    root_chord = main_wing.xsecs[0].chord
    tip_chord = main_wing.xsecs[-1].chord

    return {
        "planform_area_m2": planform_area_m2,
        "aspect_ratio": span_m**2 / planform_area_m2,
        "taper_ratio": tip_chord / root_chord,
        "sweep_25_rad": airplane.sweep_25_rad,
        "root_thickness_to_chord": airplane.root_thickness_to_chord,
        "span_m": span_m,
        "fuselage_length_m": fuselage_length_m,
        "fuselage_height_m": fuselage_height_m,
        "fuselage_width_m": fuselage_width_m,
        "fuselage_fineness_ratio": 2.0
        * fuselage_length_m
        / (fuselage_height_m + fuselage_width_m),
        "horizontal_tail_area_m2": horizontal_tail_area_m2,
        "horizontal_tail_span_m": horizontal_tail_span_m,
        "vertical_tail_area_m2": vertical_tail_area_m2,
        "vertical_tail_aspect_ratio": vertical_tail_height_m**2
        / vertical_tail_area_m2,
        "vertical_tail_height_m": vertical_tail_height_m,
        "tail_length_m": airplane.tail_length_m,
        "rudder_area_m2": airplane.rudder_area_m2,
        "wing_mounted_control_area_m2": airplane.wing_mounted_control_area_m2,
    }


def build_geometric_asb_airplane(
    *,
    planform_area_m2=80.0,
    fuselage_length_m=30.0,
    fuselage_height_m=3.6,
    fuselage_width_m=3.0,
    vtail_area_m2=20.8,
    vtail_span_m=5.12,
    vtail_dihedral_angle_deg=37.0,
    main_wing_tip_le_x_m=2.066,
    vtail_le_x_m=16.5,
    tail_length_m=16.5,
    rudder_area_m2=2.0,
    wing_mounted_control_area_m2=6.4,
    sweep_25_rad=np.radians(60.0),
    root_thickness_to_chord=0.06,
    main_gear_length_in=42.0,
    nose_gear_length_in=30.0,
    name="Sizing Aircraft",
    aspect_ratio=3.0,
    taper_ratio=0.25,
):
    """Build sizing geometry using AeroSandbox objects."""
    span_m = (planform_area_m2 * aspect_ratio) ** 0.5
    root_chord_m = 2.0 * planform_area_m2 / (span_m * (1.0 + taper_ratio))
    tip_chord_m = taper_ratio * root_chord_m

    main_wing = asb.Wing(
        name="Main Wing",
        symmetric=True,
        xsecs=[
            asb.WingXSec(
                xyz_le=[0.0, 0.0, 0.0],
                chord=root_chord_m,
                airfoil=asb.Airfoil("naca0008"),
            ),
            asb.WingXSec(
                xyz_le=[main_wing_tip_le_x_m, span_m / 2.0, 0.0],
                chord=tip_chord_m,
                airfoil=asb.Airfoil("naca0008"),
            ),
        ],
    )

    vtail_chord_m = vtail_area_m2 / vtail_span_m
    vtail_rotation = np.rotation_matrix_3D(
        angle=np.radians(vtail_dihedral_angle_deg),
        axis="X",
    )
    vtail = asb.Wing(
        name="VTail",
        symmetric=True,
        xsecs=[
            asb.WingXSec(
                xyz_le=vtail_rotation @ np.array([0.0, 0.0, 0.0]),
                chord=vtail_chord_m,
                airfoil=asb.Airfoil("naca0008"),
            ),
            asb.WingXSec(
                xyz_le=vtail_rotation @ np.array([0.0, vtail_span_m / 2.0, 0.0]),
                chord=vtail_chord_m,
                airfoil=asb.Airfoil("naca0008"),
            ),
        ],
    ).translate([vtail_le_x_m, 0.0, 0.0])

    fuselage = asb.Fuselage(
        name="Fuselage",
        xsecs=[
            asb.FuselageXSec(
                xyz_c=[0.0, 0.0, 0.0],
                width=fuselage_width_m,
                height=fuselage_height_m,
            ),
            asb.FuselageXSec(
                xyz_c=[fuselage_length_m, 0.0, 0.0],
                width=fuselage_width_m,
                height=fuselage_height_m,
            ),
        ],
    )

    airplane = asb.Airplane(
        name=name,
        wings=[
            main_wing,
            vtail,
        ],
        fuselages=[
            fuselage,
        ],
    )
    airplane.tail_length_m = tail_length_m
    airplane.rudder_area_m2 = rudder_area_m2
    airplane.wing_mounted_control_area_m2 = wing_mounted_control_area_m2
    airplane.vtail_dihedral_angle_deg = vtail_dihedral_angle_deg
    airplane.sweep_25_rad = sweep_25_rad
    airplane.root_thickness_to_chord = root_thickness_to_chord
    airplane.main_gear_length_in = main_gear_length_in
    airplane.nose_gear_length_in = nose_gear_length_in
    return airplane


def main():
    aircraft = Aircraft.from_json()
    geometry = aircraft.geometry
    weight_inputs = aircraft.to_weight_inputs()
    volume_inputs = aircraft.to_volume_inputs()

    print(" aircraft")
    print(f"Name: {aircraft.airplane.name}")
    print(f"Wings: {len(aircraft.airplane.wings)}")
    print(f"Fuselages: {len(aircraft.airplane.fuselages)}")
    print(f"Planform area: {geometry['planform_area_m2']:.3f} m^2")
    print(f"Aspect ratio: {geometry['aspect_ratio']:.3f}")
    print(f"Fuselage length: {geometry['fuselage_length_m']:.3f} m")
    print(f"Fuselage height: {geometry['fuselage_height_m']:.3f} m")
    print(f"Fuselage width: {geometry['fuselage_width_m']:.3f} m")
    print(f"Fineness ratio: {geometry['fuselage_fineness_ratio']:.3f}")
    print(f"Aircraft mass: {aircraft.mass_kg:.3f} kg")
    print(f"Fuel mass: {aircraft.fuel.mass_kg:.3f} kg")
    print(f"Fuel volume: {aircraft.fuel.volume_m3:.3f} m^3")
    print(f"Propulsion mass: {aircraft.propulsion.mass_kg:.3f} kg")
    print(f"Propulsion volume: {aircraft.propulsion.volume_m3:.3f} m^3")
    print(f"Payload mass: {aircraft.payload.mass_kg:.3f} kg")
    print(f"Payload volume: {aircraft.payload.volume_m3:.3f} m^3")
    print(f"Weight adapter: {type(weight_inputs).__name__}")
    print(f"Volume adapter: {type(volume_inputs).__name__}")


if __name__ == "__main__":
    main()
