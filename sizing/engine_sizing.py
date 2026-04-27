"""Engine sizing rules for the Astromechanic aircraft.

The sizing rule is the thrust-to-weight constraint-analysis master equation
from Zhang et al., "An Improved Method for Initial Sizing of Airbreathing
Hypersonic Aircraft," Aerospace 2023, Sec. 3.2, Eq. 14.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class EngineSizingInputs:
    """Inputs for design-point thrust sizing.

    Use a consistent force unit for `takeoff_gross_weight` and the returned
    thrust. SI is recommended: N for forces, m^2 for area, Pa for dynamic
    pressure, m/s for velocity, and W/N for specific excess power.
    """

    takeoff_gross_weight: object
    planform_area: object
    dynamic_pressure: object
    velocity: object

    installed_full_throttle_thrust_lapse: object = 1.0
    instantaneous_weight_fraction: object = 1.0
    load_factor: object = 1.0
    drag_polar_k1: object = 0.0
    drag_polar_k2: object = 0.0
    zero_lift_drag_coefficient: object = 0.0
    specific_excess_power: object = 0.0


def design_point_thrust_to_weight(inputs: EngineSizingInputs):
    """Return T_dp / W_to from the constraint-analysis master equation."""
    x = inputs
    wing_loading = x.takeoff_gross_weight / x.planform_area
    return design_point_thrust_to_weight_from_wing_loading(
        wing_loading=wing_loading,
        dynamic_pressure=x.dynamic_pressure,
        velocity=x.velocity,
        installed_full_throttle_thrust_lapse=x.installed_full_throttle_thrust_lapse,
        instantaneous_weight_fraction=x.instantaneous_weight_fraction,
        load_factor=x.load_factor,
        drag_polar_k1=x.drag_polar_k1,
        drag_polar_k2=x.drag_polar_k2,
        zero_lift_drag_coefficient=x.zero_lift_drag_coefficient,
        specific_excess_power=x.specific_excess_power,
    )


def design_point_thrust_to_weight_from_wing_loading(
    wing_loading,
    dynamic_pressure,
    velocity,
    installed_full_throttle_thrust_lapse=1.0,
    instantaneous_weight_fraction=1.0,
    load_factor=1.0,
    drag_polar_k1=0.0,
    drag_polar_k2=0.0,
    zero_lift_drag_coefficient=0.0,
    specific_excess_power=0.0,
):
    """Return T_dp / W_to as a function of wing loading W_to / S_plan."""
    beta = instantaneous_weight_fraction
    alpha = installed_full_throttle_thrust_lapse
    wing_loading_term = dynamic_pressure / (beta * wing_loading)
    lift_loading_term = load_factor * beta * wing_loading / dynamic_pressure

    # Source: Zhang et al., Aerospace 2023, Sec. 3.2, Eq. 14. The paper's K2
    # term is retained as `drag_polar_k2`; default is 0.0 for the current model.
    return (
        beta
        / alpha
        * (
            wing_loading_term
            * (
                drag_polar_k1 * lift_loading_term**2
                + drag_polar_k2 * lift_loading_term
                + zero_lift_drag_coefficient
            )
            + specific_excess_power / velocity
        )
    )


def design_point_thrust_to_weight_expanded(inputs: EngineSizingInputs):
    """Return T_dp / W_to using the equation written in W_to and S_plan form."""
    x = inputs
    beta = x.instantaneous_weight_fraction
    alpha = x.installed_full_throttle_thrust_lapse
    wing_loading_term = x.dynamic_pressure * x.planform_area / (
        beta * x.takeoff_gross_weight
    )
    lift_loading_term = (
        x.load_factor
        * beta
        * x.takeoff_gross_weight
        / (x.dynamic_pressure * x.planform_area)
    )

    # Source: Zhang et al., Aerospace 2023, Sec. 3.2, Eq. 14. The paper's K2
    # term is retained as `drag_polar_k2`; default is 0.0 for the current model.
    return (
        beta
        / alpha
        * (
            wing_loading_term
            * (
                x.drag_polar_k1 * lift_loading_term**2
                + x.drag_polar_k2 * lift_loading_term
                + x.zero_lift_drag_coefficient
            )
            + x.specific_excess_power / x.velocity
        )
    )


def design_point_thrust(inputs: EngineSizingInputs):
    """Return design-point thrust in the same force unit as W_to."""
    return design_point_thrust_to_weight(inputs) * inputs.takeoff_gross_weight


def main():
    # Edit run options here.
    inputs = EngineSizingInputs(
        takeoff_gross_weight=100000.0,
        planform_area=50.0,
        dynamic_pressure=5000.0,
        velocity=250.0,
        installed_full_throttle_thrust_lapse=1.0,
        instantaneous_weight_fraction=1.0,
        load_factor=1.0,
        drag_polar_k1=0.05,
        drag_polar_k2=0.0,
        zero_lift_drag_coefficient=0.02,
        specific_excess_power=0.0,
    )
    thrust_to_weight = design_point_thrust_to_weight(inputs)
    thrust = design_point_thrust(inputs)

    print("Astromechanic engine sizing")
    print(f"Design thrust-to-weight: {thrust_to_weight:.6g}")
    print(f"Design-point thrust: {thrust:.6g}")


if __name__ == "__main__":
    main()
