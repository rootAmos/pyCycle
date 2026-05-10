"""Solve wing loading for the minimum required thrust-to-weight ratio."""

from dataclasses import dataclass

import aerosandbox as asb
import aerosandbox.numpy as np


@dataclass(frozen=True)
class ConstraintAnalysisInputs:
    """Vectorized flight conditions for the thrust-to-weight sizing equation.

    Flight speed is specified by Mach and altitude. Dynamic pressure and
    true airspeed are computed with AeroSandbox's atmosphere model.
    """

    names: object
    mach: object
    altitude_m: object
    installed_full_throttle_thrust_lapse: object = 1.0
    instantaneous_weight_fraction: object = 1.0
    load_factor: object = 1.0
    drag_polar_k1: object = 0.0
    drag_polar_k2: object = 0.0
    zero_lift_drag_coefficient: object = 0.0
    climb_rate_m_s: object = 0.0
    acceleration_m_s2: object = 0.0
    gravity_m_s2: object = 9.80665

    @property
    def velocity_m_s(self):
        atmosphere = asb.Atmosphere(altitude=self.altitude_m)
        return self.mach * atmosphere.speed_of_sound()

    @property
    def dynamic_pressure_Pa(self):
        atmosphere = asb.Atmosphere(altitude=self.altitude_m)
        velocity_m_s = self.mach * atmosphere.speed_of_sound()
        return 0.5 * atmosphere.density() * velocity_m_s**2


@dataclass(frozen=True)
class EngineSizingResult:
    wing_loading: float
    thrust_to_weight: float
    active_condition: str
    condition_thrust_to_weight: object


def design_point_thrust_to_weight_from_wing_loading(
    wing_loading_N_m2,
    dynamic_pressure_Pa,
    velocity_m_s,
    installed_full_throttle_thrust_lapse=1.0,
    instantaneous_weight_fraction=1.0,
    load_factor=1.0,
    drag_polar_k1=0.0,
    drag_polar_k2=0.0,
    zero_lift_drag_coefficient=0.0,
    climb_rate_m_s=0.0,
    acceleration_m_s2=0.0,
    gravity_m_s2=9.80665,
):
    """Return design T/W for a condition at a given wing loading W/S."""
    beta = instantaneous_weight_fraction
    alpha = installed_full_throttle_thrust_lapse
    wing_loading_term = dynamic_pressure_Pa / (beta * wing_loading_N_m2)
    lift_loading_term = load_factor * beta * wing_loading_N_m2 / dynamic_pressure_Pa
    excess_power_m_s = climb_rate_m_s + velocity_m_s * acceleration_m_s2 / gravity_m_s2

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
            + excess_power_m_s / velocity_m_s
        )
    )


def _as_condition_column(value):
    try:
        len(value)
    except TypeError:
        return value

    return np.reshape(np.array(value), (-1, 1))


def thrust_to_weight_constraint_residual(
    thrust_to_weight,
    wing_loadings,
    inputs: ConstraintAnalysisInputs,
):
    """Return T/W - required T/W over a range of W/S values for plotting.

    Positive values satisfy the constraint. The result is shaped as
    [condition, wing_loading].
    """
    wing_loading_grid = np.reshape(np.array(wing_loadings), (1, -1))

    required_thrust_to_weight = design_point_thrust_to_weight_from_wing_loading(
        wing_loading_N_m2=wing_loading_grid,
        dynamic_pressure_Pa=_as_condition_column(inputs.dynamic_pressure_Pa),
        velocity_m_s=_as_condition_column(inputs.velocity_m_s),
        installed_full_throttle_thrust_lapse=_as_condition_column(
            inputs.installed_full_throttle_thrust_lapse
        ),
        instantaneous_weight_fraction=_as_condition_column(
            inputs.instantaneous_weight_fraction
        ),
        load_factor=_as_condition_column(inputs.load_factor),
        drag_polar_k1=_as_condition_column(inputs.drag_polar_k1),
        drag_polar_k2=_as_condition_column(inputs.drag_polar_k2),
        zero_lift_drag_coefficient=_as_condition_column(
            inputs.zero_lift_drag_coefficient
        ),
        climb_rate_m_s=_as_condition_column(inputs.climb_rate_m_s),
        acceleration_m_s2=_as_condition_column(inputs.acceleration_m_s2),
        gravity_m_s2=_as_condition_column(inputs.gravity_m_s2),
    )

    return thrust_to_weight - required_thrust_to_weight


def solve_minimum_thrust_to_weight(
    inputs: ConstraintAnalysisInputs,
    wing_loading_min=1000.0,
    wing_loading_max=20000.0,
    wing_loading_guess=None,
):
    """Find the W/S that minimizes the maximum required T/W over conditions."""
    if len(inputs.names) == 0:
        raise ValueError("At least one engine sizing condition is required.")

    if wing_loading_guess is None:
        wing_loading_guess = 0.5 * (wing_loading_min + wing_loading_max)

    opti = asb.Opti()
    wing_loading = opti.variable(init_guess=wing_loading_guess)
    thrust_to_weight = opti.variable(init_guess=0.5)

    opti.subject_to(wing_loading >= wing_loading_min)
    opti.subject_to(wing_loading <= wing_loading_max)

    required_thrust_to_weight = design_point_thrust_to_weight_from_wing_loading(
        wing_loading_N_m2=wing_loading,
        dynamic_pressure_Pa=inputs.dynamic_pressure_Pa,
        velocity_m_s=inputs.velocity_m_s,
        installed_full_throttle_thrust_lapse=(
            inputs.installed_full_throttle_thrust_lapse
        ),
        instantaneous_weight_fraction=inputs.instantaneous_weight_fraction,
        load_factor=inputs.load_factor,
        drag_polar_k1=inputs.drag_polar_k1,
        drag_polar_k2=inputs.drag_polar_k2,
        zero_lift_drag_coefficient=inputs.zero_lift_drag_coefficient,
        climb_rate_m_s=inputs.climb_rate_m_s,
        acceleration_m_s2=inputs.acceleration_m_s2,
        gravity_m_s2=inputs.gravity_m_s2,
    )
    opti.subject_to(thrust_to_weight >= required_thrust_to_weight)

    opti.minimize(thrust_to_weight)
    sol = opti.solve()

    solved_wing_loading = float(sol(wing_loading))
    condition_values = np.array(
        design_point_thrust_to_weight_from_wing_loading(
            wing_loading_N_m2=solved_wing_loading,
            dynamic_pressure_Pa=inputs.dynamic_pressure_Pa,
            velocity_m_s=inputs.velocity_m_s,
            installed_full_throttle_thrust_lapse=(
                inputs.installed_full_throttle_thrust_lapse
            ),
            instantaneous_weight_fraction=inputs.instantaneous_weight_fraction,
            load_factor=inputs.load_factor,
            drag_polar_k1=inputs.drag_polar_k1,
            drag_polar_k2=inputs.drag_polar_k2,
            zero_lift_drag_coefficient=inputs.zero_lift_drag_coefficient,
            climb_rate_m_s=inputs.climb_rate_m_s,
            acceleration_m_s2=inputs.acceleration_m_s2,
            gravity_m_s2=inputs.gravity_m_s2,
        )
    )
    active_condition = inputs.names[int(np.argmax(condition_values))]

    return EngineSizingResult(
        wing_loading=solved_wing_loading,
        thrust_to_weight=float(sol(thrust_to_weight)),
        active_condition=active_condition,
        condition_thrust_to_weight=condition_values,
    )


def main():
    import matplotlib.pyplot as plt

    wing_loading_min = 1000.0
    wing_loading_max = 20000.0

    inputs = ConstraintAnalysisInputs(
        names=("Cruise", "Climb", "Maneuver"),
        mach=np.array([3.0, 0.9, 1.2]),
        altitude_m=np.array([60000.0, 10000.0, 15000.0]) * 0.3048,
        load_factor=np.array([1.0, 1.0, 2.5]),
        drag_polar_k1=np.array([0.050, 0.055, 0.050]),
        zero_lift_drag_coefficient=np.array([0.020, 0.024, 0.026]),
        climb_rate_m_s=np.array([0.0, 12.0, 0.0]),
        acceleration_m_s2=np.array([0.0, 0.0, 0.0]),
    )

    result = solve_minimum_thrust_to_weight(
        inputs=inputs,
        wing_loading_min=wing_loading_min,
        wing_loading_max=wing_loading_max,
    )

    print("Engine sizing")
    print(f"Optimal W/S: {result.wing_loading:.6g}")
    print(f"Minimum design T/W: {result.thrust_to_weight:.6g}")
    print(f"Active condition: {result.active_condition}")
    print("Condition T/W values:")
    for name, thrust_to_weight in zip(inputs.names, result.condition_thrust_to_weight):
        print(f"  {name}: {thrust_to_weight:.6g}")

    wing_loadings = np.linspace(wing_loading_min, wing_loading_max, 300)
    wing_loading_grid = np.reshape(wing_loadings, (1, -1))
    required_curves = design_point_thrust_to_weight_from_wing_loading(
        wing_loading_N_m2=wing_loading_grid,
        dynamic_pressure_Pa=_as_condition_column(inputs.dynamic_pressure_Pa),
        velocity_m_s=_as_condition_column(inputs.velocity_m_s),
        installed_full_throttle_thrust_lapse=_as_condition_column(
            inputs.installed_full_throttle_thrust_lapse
        ),
        instantaneous_weight_fraction=_as_condition_column(
            inputs.instantaneous_weight_fraction
        ),
        load_factor=_as_condition_column(inputs.load_factor),
        drag_polar_k1=_as_condition_column(inputs.drag_polar_k1),
        drag_polar_k2=_as_condition_column(inputs.drag_polar_k2),
        zero_lift_drag_coefficient=_as_condition_column(
            inputs.zero_lift_drag_coefficient
        ),
        climb_rate_m_s=_as_condition_column(inputs.climb_rate_m_s),
        acceleration_m_s2=_as_condition_column(inputs.acceleration_m_s2),
        gravity_m_s2=_as_condition_column(inputs.gravity_m_s2),
    )

    fig, ax = plt.subplots(figsize=(8, 5))
    for name, required_thrust_to_weight in zip(inputs.names, required_curves):
        ax.plot(wing_loadings, required_thrust_to_weight, label=name)

    ax.plot(
        result.wing_loading,
        result.thrust_to_weight,
        marker="*",
        color="red",
        markersize=16,
        linestyle="None",
        label="Minimum feasible T/W",
    )
    ax.axhline(
        result.thrust_to_weight,
        color="red",
        linestyle="--",
        linewidth=1.0,
        alpha=0.6,
    )
    ax.set_xlabel("Wing loading W/S [N/m^2]")
    ax.set_ylabel("Thrust-to-weight ratio T/W [-]")
    ax.set_title("Constraint analysis")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig("_engine_sizing_constraints.png", dpi=200)
    plt.close(fig)
    print("Saved constraint plot: _engine_sizing_constraints.png")


if __name__ == "__main__":
    main()
