"""
Duality Engine Cycle Model — AstroMechanica (Proprietary Concept)
=================================================================

OVERVIEW
--------
This file models the "Duality" engine across three distinct operating modes
using NASA's pyCycle thermodynamic cycle analysis library (built on OpenMDAO).

Each mode has a fundamentally different gas-path topology (different components
and flow routing), which is why THREE separate Python classes are defined rather
than one class with switches. Trying to model all modes in a single class with
a splitter/mixer at near-zero bypass ratio causes singular Jacobians in the
CEA (Chemical Equilibrium with Applications) thermo solver — the flow fractions
become undefined. Separate classes avoid this entirely.

THREE OPERATING MODES
---------------------

  Mode 1 — Ducted Fan  (pure electric propulsion, no combustion)
    Flow path:  FC → Inlet → Fan1 → Fan2 → Duct → Nozzle
    Drive:      Electric motors spin Fan1 and Fan2 (no turbine needed).
    Use:        Low-speed / subsonic operation.

  Mode 2 — Ducted Fan + Afterburner  (electric fans + fuel burn)
    Flow path:  FC → Inlet → Fan1 → Fan2 → Combustor → Nozzle
    Drive:      Same electric fans as Mode 1; fuel is injected AFTER the fans.
    Use:        High-subsonic dash, augmented thrust.

  Mode 3 — RamJet  (no fans, ram compression only)
    Flow path:  FC → Inlet → BypassDuct → Combustor → Nozzle
    Drive:      Fans are stopped; ram pressure from supersonic flight compresses air.
    Use:        Supersonic cruise / acceleration.

DESIGN POINT vs OFF-DESIGN
--------------------------
pyCycle distinguishes two analysis modes:
  - DESIGN (on-design):  The solver finds what physical geometry (areas, map
                         scaling factors) is required to hit a thrust target.
                         This "sizes" the engine.
  - OFF-DESIGN (OD):     Given the fixed geometry from DESIGN, the solver finds
                         what flow conditions (W, N) naturally occur at a
                         different flight condition.

ELECTRIC MOTOR DRIVES
---------------------
Fan1 and Fan2 are electrically driven — there is no turbine and no mechanical
shaft connecting compressors to a turbine (unlike a conventional jet engine).
Fan1 and Fan2 are COUNTER-ROTATING — they spin in opposite physical directions,
which eliminates net exit swirl and allows each stage to be independently
controlled.  Each fan has its own speed variable (N_fan1 for fan1, N_fan2 for
fan2) driven by a separate motor.  In the 1-D thermodynamic model the direction
of rotation does not affect the map (maps are defined in terms of |Nc|), so
both fans use the same FanMap but with independent speed states.
Motor power demand is read from fan.power output — no shaft power balance
equation is needed because there is no turbine to balance against.

SOLVER APPROACH
---------------
Each Cycle uses Newton's method (om.NewtonSolver) with a direct linear solver
(om.DirectSolver — factorizes the full Jacobian). The balance equations drive
specific unknowns (W, N_fan, FAR) to match physical constraints (throat area,
operating line, temperature targets).
"""

import sys
import math
import os

os.environ.setdefault("OPENMDAO_REPORTS", "0")
import openmdao.api as om   # OpenMDAO: the multidisciplinary optimisation framework
import pycycle.api as pyc   # pyCycle: thermodynamic cycle components built on OpenMDAO


# ============================================================================
# Shared solver helper
# ============================================================================

def _add_newton(system, atol=1e-6, rtol=1e-6, maxiter=50):
    """
    Attach a Newton nonlinear solver + direct linear solver to `system`.

    This is called at the end of every Cycle's setup() to configure how the
    coupled nonlinear system of equations is solved.

    Parameters
    ----------
    system   : OpenMDAO System (the Cycle instance)
    atol     : absolute residual tolerance — solver stops when ||R|| < atol
    rtol     : relative residual tolerance — solver stops when ||R||/||R0|| < rtol
    maxiter  : maximum Newton iterations before giving up

    Newton solver settings explained
    ---------------------------------
    iprint=2              → print residual at every iteration (useful for debugging)
    solve_subsystems=True → before the outer Newton step, run each subsystem's own
                            solver to get a consistent starting state.  This is
                            critical for convergence when subsystem residuals are large.
    max_sub_solves=100    → allow up to 100 subsystem solves per Newton iteration
                            (generous limit to avoid premature failure)
    reraise_child_analysiserror=False
                          → if a subsystem throws an AnalysisError (e.g., CEA
                            failed for an unphysical state), catch it and
                            continue rather than crashing the whole run.

    Line search (ArmijoGoldstein)
    -----------------------------
    Plain Newton can overshoot badly for thermodynamic problems — a full Newton
    step might walk into a region where CEA cannot evaluate (negative pressures,
    temperatures above dissociation limits).  The Armijo-Goldstein line search
    starts with the full Newton step and backtracks (multiplying the step size
    by `rho`=0.75 each time) until the residual actually decreases.
    iprint=-1 → suppress line-search iteration printing to keep output clean.

    Linear solver
    -------------
    om.DirectSolver factorizes the full assembled Jacobian matrix using LU
    decomposition.  This is exact (no iterative approximation) and works well
    for the O(100) variable systems typical of these cycle models.
    """
    # --- Nonlinear solver: Newton ---
    newton = system.nonlinear_solver = om.NewtonSolver()
    newton.options['atol']      = atol      # stop if |residual| < 1e-6
    newton.options['rtol']      = rtol      # stop if relative residual < 1e-6
    newton.options['iprint']    = 2         # print residual every iteration
    newton.options['maxiter']   = maxiter   # max 50 Newton steps
    newton.options['solve_subsystems']          = True   # pre-solve each sub-system
    newton.options['max_sub_solves']            = 100    # generous subsystem solve limit
    newton.options['reraise_child_analysiserror'] = False # survive CEA failures

    # --- Line search: Armijo-Goldstein backtracking ---
    newton.linesearch = om.ArmijoGoldsteinLS()
    newton.linesearch.options['rho']    = 0.75  # backtrack factor (step *= 0.75 each try)
    newton.linesearch.options['iprint'] = -1    # suppress line-search printing

    # --- Linear solver: LU factorisation of the full Jacobian ---
    system.linear_solver = om.DirectSolver()


# ============================================================================
# Mode 1: Ducted Fan  (FC → Inlet → Fan1 → Fan2 → Duct → Nozzle)
# ============================================================================

class DualityFanOnly(pyc.Cycle):
    """
    Mode 1: Fan-only (pure electric, no combustion).

    Gas path:  FC → Inlet → Fan1 → Fan2 → Duct(ab) → Nozzle
    Drive:     Counter-rotating electric motors (N_fan1 for fan1, N_fan2 for fan2).
    'ab' element: plain Duct (pressure loss dPqP), NOT a combustor.
               The name 'ab' is kept for slot consistency with Mode 2.

    Used for
    --------
    - DESIGN_mode1  (design=True)  — sizes the engine at SLS, sets all fan areas
    - OD_mode1      (design=False) — subsonic cruise analysis

    BALANCE EQUATIONS
    -----------------
    DESIGN mode (design=True):
      Variable: W  (inlet mass flow rate, lbm/s)
      Equation: perf.Fn == balance.rhs:W
      N_fan1 and N_fan2 are FIXED design inputs (set via set_input_defaults).

    OFF-DESIGN mode (design=False):
      Variable 1: W  → nozz.Throat:stat:area == rhs:W  (DESIGN_mode1 throat area)
      Variable 2: N_fan1 → fan1.map.RlineMap == 2.0  (fan1 independent speed)
      Variable 3: N_fan2 → fan2.map.RlineMap == 2.0  (fan2 independent speed)
      Each fan independently seeks its operating line — correct for independently
      driven counter-rotating stages.
    """

    def setup(self):
        # 'design' is an OpenMDAO option set when instantiating the class
        # (True for DESIGN point, False for all OD points)
        design = self.options['design']

        # --- Thermodynamic library selection ---
        # CEA (Chemical Equilibrium with Applications) is NASA's high-fidelity
        # thermochemistry solver.  It calculates mixture enthalpy, entropy,
        # specific heats, etc. accounting for real-gas effects and dissociation.
        # janaf = Joint Army-Navy-Air Force thermodynamic data tables — the
        # standard dataset of species properties used by CEA.
        self.options['thermo_method'] = 'CEA'
        self.options['thermo_data']   = pyc.species_data.janaf

        # -----------------------------------------------------------------------
        # ADD COMPONENTS (the gas-path elements, in flow order)
        # -----------------------------------------------------------------------

        # FlightConditions: computes freestream total/static state from altitude
        # and Mach number using the International Standard Atmosphere (ISA).
        # Outputs: Fl_O (flow object with Tt, Pt, Ts, Ps, MN, W)
        self.add_subsystem('fc',    pyc.FlightConditions())

        # Inlet: models the intake ram recovery.
        # Key parameter: ram_recovery (Pt_out/Pt_in) — set globally via
        # pyc_add_cycle_param.  Also computes F_ram = inlet momentum drag.
        self.add_subsystem('inlet', pyc.Inlet())

        # Fan1: first-stage electrically-driven fan.
        #   map_data=pyc.FanMap   → use pyCycle's built-in fan map (scaled by s_PR etc.)
        #   map_extrap=True       → allow extrapolation beyond the map table edges
        #                           (important for OD points that wander off-design)
        # promotes_inputs=[('Nmech','N_fan1')]:
        #   The compressor's internal shaft-speed input is called 'Nmech'.
        #   We rename it to 'N_fan1' so fan1 gets its own independent speed
        #   variable — distinct from fan2's N_fan2.  This allows counter-rotating
        #   stages to be driven at different speeds.
        self.add_subsystem('fan1',  pyc.Compressor(map_data=pyc.FanMap, map_extrap=True),
                           promotes_inputs=[('Nmech', 'N_fan1')])

        # Fan2: counter-rotating second stage, independent speed N_fan2.
        # In DESIGN mode N_fan1 and N_fan2 are both fixed inputs (typically equal).
        # In OD mode each has its own balance equation targeting RlineMap==2.0.
        self.add_subsystem('fan2',  pyc.Compressor(map_data=pyc.FanMap, map_extrap=True),
                           promotes_inputs=[('Nmech', 'N_fan2')])

        # ab (afterburner slot): in Mode 1 this is just a plain Duct — it
        # passes flow through with a fractional pressure loss (dPqP).
        # No fuel is added here.  The name 'ab' is kept for consistency with
        # Mode 2 which uses a Combustor in this same slot.
        self.add_subsystem('ab',    pyc.Duct())          # pass-through, no fuel

        # Nozzle: converts total enthalpy to kinetic energy, produces thrust.
        #   nozzType='CV' → convergent nozzle (area decreases to throat, no diverging section)
        #   lossCoef='Cv' → loss specified as velocity coefficient Cv = V_actual/V_ideal
        #                   (Cv < 1.0 accounts for friction/mixing losses in nozzle)
        self.add_subsystem('nozz',  pyc.Nozzle(nozzType='CV', lossCoef='Cv'))

        # Performance: calculates overall engine figures of merit.
        #   num_nozzles=1  → one nozzle contributing to thrust
        #   num_burners=0  → no burners (pure fan mode, no TSFC calculation)
        self.add_subsystem('perf',  pyc.Performance(num_nozzles=1, num_burners=0))

        # -----------------------------------------------------------------------
        # CONNECT FLOW STATIONS (pyc_connect_flow propagates thermodynamic state)
        # -----------------------------------------------------------------------
        # pyc_connect_flow(upstream, downstream, ...)
        # A "flow connection" passes total state (Tt, Pt, composition) and
        # optionally static state (Ts, Ps, area) and/or mass flow (W).
        #
        # connect_w=False: DO NOT connect mass flow from fc to inlet.
        #   Why? The Inlet is the first "real" component.  The mass flow W is
        #   a free variable (solved by the balance), not imposed by FlightConditions.
        #   FlightConditions only sets the stagnation state (Pt, Tt) of the
        #   incoming airstream; W is determined by the balance equation.
        #
        # connect_stat=False: DO NOT connect static state between internal stations.
        #   Why? Static state (Ts, Ps, MN) depends on the local duct area.  Each
        #   component calculates its own exit static state from its geometry.
        #   Propagating static state upstream-to-downstream would over-constrain
        #   the system.  Only total state (Tt, Pt, composition) and W are passed.
        self.pyc_connect_flow('fc.Fl_O',    'inlet.Fl_I',  connect_w=False)
        self.pyc_connect_flow('inlet.Fl_O', 'fan1.Fl_I',   connect_stat=False)
        self.pyc_connect_flow('fan1.Fl_O',  'fan2.Fl_I',   connect_stat=False)
        self.pyc_connect_flow('fan2.Fl_O',  'ab.Fl_I',     connect_stat=False)
        self.pyc_connect_flow('ab.Fl_O',    'nozz.Fl_I',   connect_stat=False)

        # -----------------------------------------------------------------------
        # SCALAR CONNECTIONS (non-flow variables needed by perf and nozzle)
        # -----------------------------------------------------------------------

        # Nozzle needs the ambient static pressure to compute the pressure thrust
        # term.  fc.Fl_O:stat:P is the freestream static pressure (psia).
        self.connect('fc.Fl_O:stat:P',   'nozz.Ps_exhaust')

        # Performance component needs certain pressures for OPR calculation:
        #   Pt2 = total pressure at compressor face (inlet exit = fan1 inlet)
        #   Pt3 = total pressure at combustor inlet (fan2 exit in this mode)
        self.connect('inlet.Fl_O:tot:P', 'perf.Pt2')
        self.connect('fan2.Fl_O:tot:P',  'perf.Pt3')

        # F_ram: momentum drag of ingesting air at flight speed.
        # F_ram = W * V_flight.  Net thrust = Fg - F_ram.
        self.connect('inlet.F_ram', 'perf.ram_drag')

        # Fg: gross thrust produced by the nozzle (momentum + pressure term)
        self.connect('nozz.Fg', 'perf.Fg_0')

        # -----------------------------------------------------------------------
        # BALANCE EQUATIONS
        # -----------------------------------------------------------------------
        # om.BalanceComp is an OpenMDAO component that represents the implicit
        # equation:  LHS == RHS
        # It holds a state variable that Newton drives until LHS == RHS.
        balance = self.add_subsystem('balance', om.BalanceComp())

        if design:
            # ---- DESIGN BALANCE: find W to hit a thrust target ----
            #
            # State variable: W (inlet mass flow, lbm/s)
            #   Initial guess: val=50. lbm/s
            #   eq_units='lbf' → the equality is in pound-force
            #
            # Equation:  perf.Fn  ==  balance.rhs:W
            #   LHS: actual net thrust computed from the cycle (lbf)
            #   RHS: target thrust (set externally, e.g. ~1100 lbf)
            #   → Newton adjusts W until Fn == target
            balance.add_balance('W', units='lbm/s', eq_units='lbf', val=50.)
            self.connect('balance.W', 'inlet.Fl_I:stat:W')  # W drives mass flow
            self.connect('perf.Fn',   'balance.lhs:W')       # LHS = computed thrust

        else:
            # ---- OFF-DESIGN BALANCE 1: W — match frozen throat area ----
            # throat area RHS connected from DESIGN_mode1.nozz.Throat:stat:area
            balance.add_balance('W', val=50., units='lbm/s', eq_units='inch**2')
            self.connect('balance.W',             'inlet.Fl_I:stat:W')
            self.connect('nozz.Throat:stat:area', 'balance.lhs:W')

            # ---- OFF-DESIGN BALANCE 2: inlet area — hold a scheduled diffuser exit MN ----
            balance.add_balance('inlet_area', val=260., units='inch**2',
                                lower=180., upper=400., eq_units=None)
            self.connect('balance.inlet_area',   'inlet.area')
            self.connect('inlet.Fl_O:stat:MN',   'balance.lhs:inlet_area')
            # Fan speeds are scheduled independently in Mode 1, so off-design
            # solves only for mass flow against the mode-specific throat area.

        # -----------------------------------------------------------------------
        # EXECUTION ORDER & SOLVER
        # -----------------------------------------------------------------------
        # set_order ensures components execute in physical flow order, which
        # gives Newton a better-structured residual (each component's outputs
        # feed naturally into the next).  'balance' goes last because it depends
        # on outputs from all flow components.
        self.set_order(['fc', 'inlet', 'fan1', 'fan2', 'ab', 'nozz', 'perf', 'balance'])

        # Attach Newton solver + DirectSolver (see _add_newton docstring)
        _add_newton(self)

        # super().setup() finalises pyCycle's internal connection machinery
        # (CEA property tables, flow port wiring, etc.)
        super().setup()


# ============================================================================
# Mode 2: Fan + Afterburner  (FC → Inlet → Fan1 → Fan2 → Combustor → Nozzle)
# ============================================================================

class DualityFanAB(pyc.Cycle):
    """
    Mode 2: Fan + Afterburner (electric fans + fuel combustion after fans).

    Gas path:  FC → Inlet → Fan1 → Fan2 → Combustor(ab) → Nozzle
    Drive:     Counter-rotating electric fans (N_fan1, N_fan2) + Jet-A afterburner.
    'ab' element: Combustor (Jet-A fuel).

    Used for
    --------
    - DESIGN_mode2 (design=True)  — sizes the Mode 2 nozzle at M=0.8, 30 000 ft.
      Fan hardware (areas, map scalars) comes from DESIGN_mode1.
      Only the nozzle throat area is new from this run.
    - OD_mode2     (design=False) — subsonic dash analysis.

    BALANCE EQUATIONS
    -----------------
    DESIGN mode (design=True):
      Variable 1: W   → perf.Fn == rhs:W  (thrust target)
      Variable 2: FAR → ab.Fl_O:tot:T == rhs:FAR  (T4 target)
      N_fan1 and N_fan2 are fixed inputs.

    OFF-DESIGN mode (design=False):
      Variable 1: W     → nozz.Throat:stat:area == rhs:W  (DESIGN_mode2 area)
      Variable 2: N_fan1 → fan1.map.RlineMap == 2.0
      Variable 3: N_fan2 → fan2.map.RlineMap == 2.0
      Variable 4: FAR   → ab.Fl_O:tot:T == rhs:FAR
    """

    def setup(self):
        design = self.options['design']   # True = size nozzle, False = OD analysis
        self.options['thermo_method'] = 'CEA'
        self.options['thermo_data']   = pyc.species_data.janaf
        FUEL_TYPE = 'Jet-A(g)'

        # -----------------------------------------------------------------------
        # ADD COMPONENTS
        # -----------------------------------------------------------------------
        self.add_subsystem('fc',    pyc.FlightConditions())
        self.add_subsystem('inlet', pyc.Inlet())
        # Counter-rotating fans: each has its own independent speed variable.
        self.add_subsystem('fan1',  pyc.Compressor(map_data=pyc.FanMap, map_extrap=True),
                           promotes_inputs=[('Nmech', 'N_fan1')])
        self.add_subsystem('fan2',  pyc.Compressor(map_data=pyc.FanMap, map_extrap=True),
                           promotes_inputs=[('Nmech', 'N_fan2')])

        # KEY DIFFERENCE from Mode 1: 'ab' is now a Combustor, not a Duct.
        # pyc.Combustor models the fuel injection, mixing, and equilibrium
        # combustion chemistry via CEA.  It takes Fl_I:FAR as the fuel-to-air
        # ratio input and outputs Fl_O with the post-combustion total state.
        self.add_subsystem('ab',    pyc.Combustor(fuel_type=FUEL_TYPE))

        self.add_subsystem('nozz',  pyc.Nozzle(nozzType='CV', lossCoef='Cv'))

        # num_burners=1: tell Performance that one burner exists so it can
        # compute TSFC (thrust-specific fuel consumption = Wfuel / Fn)
        self.add_subsystem('perf',  pyc.Performance(num_nozzles=1, num_burners=1))

        # -----------------------------------------------------------------------
        # CONNECT FLOW STATIONS (same as Mode 1)
        # -----------------------------------------------------------------------
        self.pyc_connect_flow('fc.Fl_O',    'inlet.Fl_I',  connect_w=False)
        self.pyc_connect_flow('inlet.Fl_O', 'fan1.Fl_I',   connect_stat=False)
        self.pyc_connect_flow('fan1.Fl_O',  'fan2.Fl_I',   connect_stat=False)
        self.pyc_connect_flow('fan2.Fl_O',  'ab.Fl_I',     connect_stat=False)
        self.pyc_connect_flow('ab.Fl_O',    'nozz.Fl_I',   connect_stat=False)

        # -----------------------------------------------------------------------
        # SCALAR CONNECTIONS
        # -----------------------------------------------------------------------
        self.connect('fc.Fl_O:stat:P',   'nozz.Ps_exhaust')
        self.connect('inlet.Fl_O:tot:P', 'perf.Pt2')
        self.connect('fan2.Fl_O:tot:P',  'perf.Pt3')

        # NEW vs Mode 1: connect fuel flow from afterburner to Performance
        # Wfuel (lbm/s) is used to compute TSFC.
        self.connect('ab.Wfuel',         'perf.Wfuel_0')

        self.connect('inlet.F_ram', 'perf.ram_drag')
        self.connect('nozz.Fg',     'perf.Fg_0')

        # -----------------------------------------------------------------------
        # BALANCE EQUATIONS
        # -----------------------------------------------------------------------
        balance = self.add_subsystem('balance', om.BalanceComp())

        if design:
            # DESIGN: find W for thrust target; find FAR for T4 target.
            # N_fan1 and N_fan2 are fixed inputs (set_input_defaults in MPDuality).
            balance.add_balance('W', units='lbm/s', eq_units='lbf', val=40.)
            self.connect('balance.W', 'inlet.Fl_I:stat:W')
            self.connect('perf.Fn',   'balance.lhs:W')
            # rhs:W = thrust target, connected externally

            # FAR balance: find fuel-air ratio that hits T4 target.
            # Same logic as OD; in DESIGN this sizes the combustor operating line.
            balance.add_balance('FAR', eq_units='degR', lower=1e-4, val=0.017)
            self.connect('balance.FAR',    'ab.Fl_I:FAR')
            self.connect('ab.Fl_O:tot:T', 'balance.lhs:FAR')
            # rhs:FAR = T4 target, connected externally

        else:
            # OD: W locks throat area; inlet area is scheduled within bounds to
            # hold a realistic diffuser exit Mach; N_fan1/N_fan2 lock each fan's
            # operating line; FAR drives the afterburner to the temperature target.

            # Balance 1: W — nozzle throat area from DESIGN_mode2
            balance.add_balance('W', val=50., units='lbm/s', eq_units='inch**2')
            self.connect('balance.W',             'inlet.Fl_I:stat:W')
            self.connect('nozz.Throat:stat:area', 'balance.lhs:W')
            # rhs:W connected from DESIGN_mode2.nozz.Throat:stat:area in MPDuality

            # Balance 2: inlet area — keep diffuser exit Mach in a plausible range
            balance.add_balance('inlet_area', val=260., units='inch**2',
                                lower=180., upper=450., eq_units=None)
            self.connect('balance.inlet_area', 'inlet.area')
            self.connect('inlet.Fl_O:stat:MN', 'balance.lhs:inlet_area')

            # Balance 3: N_fan1 — fan1 operating line
            # N_fan1 is set directly on the operating point.
            # Fan speeds are prescribed by the motor speed schedule.

            # Balance 4: N_fan2 — fan2 operating line (independent counter-rotating)
            # N_fan2 is set directly on the operating point.

            # Balance 5: FAR — afterburner exit temperature target
            balance.add_balance('FAR', eq_units='degR', lower=1e-4, val=0.017)
            self.connect('balance.FAR',    'ab.Fl_I:FAR')
            self.connect('ab.Fl_O:tot:T', 'balance.lhs:FAR')

        self.set_order(['fc', 'inlet', 'fan1', 'fan2', 'ab', 'nozz', 'perf', 'balance'])
        _add_newton(self)
        super().setup()


# ============================================================================
# Mode 3: RamJet  (FC → Inlet → BypassDuct → Combustor → Nozzle)
# ============================================================================

class DualityRamjet(pyc.Cycle):
    """
    Mode 3: Ramjet (supersonic, no fans, ram compression only).

    Gas path:  FC → Inlet → BypassDuct → Combustor → Nozzle
    Drive:     None — fans are OFF.  The inlet decelerates supersonic flow,
               converting kinetic energy to pressure (ram compression).
               The bypass_duct is a low-loss passage that carries compressed
               air from the inlet to the combustor.

    WHY a different class?
    ----------------------
    At M=2.5 the ram pressure ratio (1 + 0.2·M²)^3.5 ≈ 17, far exceeding the
    fan pressure ratio (~2.0).  The fans would be ingesting reversed pressure
    gradients and stall immediately.  Separating Mode 3 as its own class cleanly
    removes fan1, fan2 from the problem rather than clamping them near-zero.

    Used for
    --------
    - DESIGN_mode3 (design=True)  — sizes the Mode 3 nozzle at M=2.5, 40 000 ft.
      Inlet area comes from DESIGN_mode1 (same physical intake hardware).
    - OD_mode3     (design=False) — supersonic ramjet analysis.

    BALANCE EQUATIONS
    -----------------
    DESIGN mode (design=True):
      Variable 1: W   → perf.Fn == rhs:W   (thrust target)
      Variable 2: FAR → combustor.Fl_O:tot:T == rhs:FAR  (T4 target)
      Internal MN values (inlet.MN, bypass_duct.MN, combustor.MN) are fixed
      inputs that determine duct cross-section areas.

    OFF-DESIGN mode (design=False):
      Variable 1: W   → nozz.Throat:stat:area == rhs:W  (DESIGN_mode3 area)
      Variable 2: FAR → combustor.Fl_O:tot:T == rhs:FAR
    """

    def setup(self):
        design = self.options['design']
        self.options['thermo_method'] = 'CEA'
        self.options['thermo_data']   = pyc.species_data.janaf
        FUEL_TYPE = 'Jet-A(g)'

        # -----------------------------------------------------------------------
        # ADD COMPONENTS
        # -----------------------------------------------------------------------
        self.add_subsystem('fc',    pyc.FlightConditions())

        # Inlet: decelerates flow from M=2.5 to subsonic (normal shock + subsonic
        # diffuser).  ram_recovery < 1 due to shock losses.
        # Note: in a real oblique-shock inlet the recovery would be modelled more
        # carefully (planned extension in docstring above).
        self.add_subsystem('inlet', pyc.Inlet())

        # bypass_duct: carries high-pressure subsonic flow from inlet to combustor.
        # Replaces the fan stages — the fans are physically stopped (or retracted).
        # Small pressure loss modelled via dPqP.
        self.add_subsystem('bypass_duct', pyc.Duct())   # replaces fan path

        # Combustor: injects and burns Jet-A in the high-pressure ram air.
        # This is the primary combustion zone (not an afterburner) — it sees
        # much higher pressure than in turbofan afterburners.
        self.add_subsystem('combustor',   pyc.Combustor(fuel_type=FUEL_TYPE))

        # Nozzle: same convergent-velocity-coefficient nozzle as other modes.
        # At M=2.5 operation the nozzle would ideally be C-D (convergent-divergent)
        # for supersonic exit, but CV is used here as a simplification.
        self.add_subsystem('nozz',        pyc.Nozzle(nozzType='CV', lossCoef='Cv'))

        # Performance: one nozzle, one burner (the combustor)
        self.add_subsystem('perf',        pyc.Performance(num_nozzles=1, num_burners=1))

        # -----------------------------------------------------------------------
        # CONNECT FLOW STATIONS (simple linear chain — no fan stages)
        # -----------------------------------------------------------------------
        self.pyc_connect_flow('fc.Fl_O',          'inlet.Fl_I',       connect_w=False)
        self.pyc_connect_flow('inlet.Fl_O',        'bypass_duct.Fl_I', connect_stat=False)
        self.pyc_connect_flow('bypass_duct.Fl_O',  'combustor.Fl_I',   connect_stat=False)
        self.pyc_connect_flow('combustor.Fl_O',    'nozz.Fl_I',        connect_stat=False)

        # -----------------------------------------------------------------------
        # SCALAR CONNECTIONS
        # -----------------------------------------------------------------------
        self.connect('fc.Fl_O:stat:P',         'nozz.Ps_exhaust')
        self.connect('inlet.Fl_O:tot:P',       'perf.Pt2')

        # Pt3 is defined as "last compressor exit pressure" for OPR calculation.
        # In ramjet mode there is no compressor, so bypass_duct exit pressure
        # (after the inlet) is used as the pseudo-Pt3.
        self.connect('bypass_duct.Fl_O:tot:P', 'perf.Pt3')

        self.connect('combustor.Wfuel',        'perf.Wfuel_0')
        self.connect('inlet.F_ram',            'perf.ram_drag')
        self.connect('nozz.Fg',                'perf.Fg_0')

        # -----------------------------------------------------------------------
        # BALANCE EQUATIONS
        # -----------------------------------------------------------------------
        balance = self.add_subsystem('balance', om.BalanceComp())

        if design:
            # DESIGN: find W to hit thrust target; find FAR to hit T4 target.
            # Internal MN inputs (set via set_input_defaults in MPDuality) determine
            # bypass_duct and combustor cross-section areas for this mode.
            balance.add_balance('W', units='lbm/s', eq_units='lbf', val=100.)
            self.connect('balance.W', 'inlet.Fl_I:stat:W')
            self.connect('perf.Fn',   'balance.lhs:W')
            # rhs:W = thrust target, connected externally

            balance.add_balance('FAR', eq_units='degR', lower=1e-4, val=0.04)
            self.connect('balance.FAR',          'combustor.Fl_I:FAR')
            self.connect('combustor.Fl_O:tot:T', 'balance.lhs:FAR')
            # rhs:FAR = T4 target, connected externally

        else:
            # OD: W locks DESIGN_mode3 throat area; FAR hits temperature target.
            balance.add_balance('W', val=50., units='lbm/s', eq_units='inch**2')
            self.connect('balance.W',             'inlet.Fl_I:stat:W')
            self.connect('nozz.Throat:stat:area', 'balance.lhs:W')
            # rhs:W connected from DESIGN_mode3.nozz.Throat:stat:area in MPDuality

            balance.add_balance('FAR', eq_units='degR', lower=1e-4, val=0.04)
            self.connect('balance.FAR',          'combustor.Fl_I:FAR')
            self.connect('combustor.Fl_O:tot:T', 'balance.lhs:FAR')

        self.set_order(['fc', 'inlet', 'bypass_duct', 'combustor', 'nozz', 'perf', 'balance'])
        _add_newton(self)
        super().setup()


# ============================================================================
# Multi-point model
# ============================================================================

class MPDuality(pyc.MPCycle):
    """
    Multi-point (MP) Duality engine model.

    pyc.MPCycle is a pyCycle container that manages multiple operating points
    (instances of pyc.Cycle) within a single OpenMDAO Problem.  It handles
    the "design → off-design scaling" connections that make all OD points
    geometrically consistent with the DESIGN point.

    OPERATING POINTS
    ----------------
    DESIGN   — Mode 2, supersonic fan + afterburner cruise: sizes the shared
                fan-mode hardware and hot nozzle geometry.
    OD_mode1 — Mode 1 fan-only check.
    OD_mode2 — Mode 2 fixed-geometry check.
    OD_mode3 — Mode 3 ramjet check.

    DESIGN → OD SCALING
    --------------------
    After the DESIGN run, pyCycle passes two types of information to OD points:

    1. Compressor map scaling factors (s_PR, s_Wc, s_eff, s_Nc):
       The fan map is a non-dimensional normalised table.  The scaling factors
       stretch/shift the map so that at design conditions the fan operates at
       PR, Wc, and efficiency equal to the design-point values.
       These scalars are computed during DESIGN and are frozen for all OD runs.
       → fan1.s_PR, fan1.s_Wc, fan1.s_eff, fan1.s_Nc connected to each OD point.

    2. Fixed flow areas (duct areas):
       Each duct/component has a fixed cross-sectional area (determined by the
       Mach number and mass flow at DESIGN).  In OD mode the solver must
       respect these fixed areas when computing static conditions.
       → inlet, fan1, fan2, ab exit areas connected from DESIGN to OD points.

    3. Nozzle throat area:
       The throat area from DESIGN is used as the RHS of the W balance in all
       OD points — the physical throat constraint that locks mass flow.
    """

    def setup(self):

        # ====================================================================
        # DESIGN POINT — Mode 2 fan + afterburner, Concorde-like cruise
        # ====================================================================
        # pyc_add_pnt: registers an operating point within the MPCycle container.
        # ====================================================================
        # DESIGN_mode2 — supersonic fan + afterburner cruise — sizes the shared
        # fan hardware and the Mode 2 nozzle
        # ====================================================================
        # This is the primary design point. It determines:
        #   - All fan duct cross-section areas (inlet, fan1, fan2, ab)
        #   - Fan compressor map scaling factors (s_PR, s_Wc, s_eff, s_Nc)
        #   - Mode 2 nozzle throat area
        self.pyc_add_pnt('DESIGN_mode2', DualityFanAB(design=True))

        # Counter-rotating fan design speeds: both start at 6000 rpm.
        # These are the reference speeds used to compute corrected speed scalars.
        self.set_input_defaults('DESIGN_mode2.N_fan1', 6000., units='rpm')
        self.set_input_defaults('DESIGN_mode2.N_fan2', 6000., units='rpm')
        # Internal Mach numbers at each station (determines duct areas)
        self.set_input_defaults('DESIGN_mode2.inlet.MN', 0.60)
        self.set_input_defaults('DESIGN_mode2.fan1.MN',  0.45)
        self.set_input_defaults('DESIGN_mode2.fan2.MN',  0.40)
        self.set_input_defaults('DESIGN_mode2.ab.MN',    0.35)
        self.set_input_defaults('DESIGN_mode2.ab.dPqP',  0.03)

        # NOTE: DESIGN_mode3 is still run as a standalone Problem in __main__
        # because pyCycle MPCycle only allows one design=True point.

        # -----------------------------------------------------------------------
        # CYCLE-LEVEL CONSTANTS (broadcast to every point that has these components)
        # -----------------------------------------------------------------------
        self.pyc_add_cycle_param('inlet.ram_recovery', 0.99)
        self.pyc_add_cycle_param('nozz.Cv',            0.99)

        # ====================================================================
        # OD_mode1 — fan-only, A220-like cruise check
        # ====================================================================
        self.pyc_add_pnt('OD_mode1', DualityFanOnly(design=False))
        self.set_input_defaults('OD_mode1.fc.MN',  CRUISE_CONDITIONS['mode1']['mach'])
        self.set_input_defaults('OD_mode1.fc.alt', CRUISE_CONDITIONS['mode1']['alt_ft'], units='ft')
        self.set_input_defaults('OD_mode1.ab.dPqP', 0.01)

        # ====================================================================
        # OD_mode2 — turbojet / afterburning mode, Concorde-like cruise check
        # ====================================================================
        self.pyc_add_pnt('OD_mode2', DualityFanAB(design=False))
        self.set_input_defaults('OD_mode2.fc.MN',  CRUISE_CONDITIONS['mode2']['mach'])
        self.set_input_defaults('OD_mode2.fc.alt', CRUISE_CONDITIONS['mode2']['alt_ft'], units='ft')
        self.set_input_defaults('OD_mode2.balance.rhs:FAR', 3200., units='degR')
        self.set_input_defaults('OD_mode2.ab.dPqP', 0.03)

        # ====================================================================
        # OD_mode3 — ramjet, SR-71-like cruise check
        # ====================================================================
        self.pyc_add_pnt('OD_mode3', DualityRamjet(design=False))
        self.set_input_defaults('OD_mode3.fc.MN',  CRUISE_CONDITIONS['mode3']['mach'])
        self.set_input_defaults('OD_mode3.fc.alt', CRUISE_CONDITIONS['mode3']['alt_ft'], units='ft')
        self.set_input_defaults('OD_mode3.balance.rhs:FAR', 3800., units='degR')
        self.set_input_defaults('OD_mode3.bypass_duct.dPqP', 0.01)
        self.set_input_defaults('OD_mode3.combustor.dPqP',   0.03)

        # ====================================================================
        # DESIGN → OD CONNECTIONS
        # ====================================================================
        #
        # Fan map scaling factors for both fan modes come from DESIGN_mode2.
        for pt in ('OD_mode1', 'OD_mode2'):
            for sfx in ('s_PR', 's_Wc', 's_eff', 's_Nc'):
                self.connect(f'DESIGN_mode2.fan1.{sfx}', f'{pt}.fan1.{sfx}')
                self.connect(f'DESIGN_mode2.fan2.{sfx}', f'{pt}.fan2.{sfx}')

        # Fan hardware is shared across both fan modes, but the inlet is
        # mode-specific and is therefore not connected here.
        for pt in ('OD_mode1', 'OD_mode2'):
            self.connect('DESIGN_mode2.fan1.Fl_O:stat:area',  f'{pt}.fan1.area')
            self.connect('DESIGN_mode2.fan2.Fl_O:stat:area',  f'{pt}.fan2.area')
            self.connect('DESIGN_mode2.ab.Fl_O:stat:area',    f'{pt}.ab.area')

        # Mode 2 nozzle area is the primary fan+AB throat constraint. Mode 1
        # gets its own nozzle area via prob.set_val() in __main__.
        self.connect('DESIGN_mode2.nozz.Throat:stat:area', 'OD_mode2.balance.rhs:W')

        # Mode 3 nozzle area is set via prob.set_val() in __main__ after the
        # standalone sizing run has been completed.

        # OD_mode3 inlet / bypass duct / combustor areas are injected from
        # standalone DESIGN_mode3 in __main__.

        super().setup()


# ============================================================================
# Cycle station map (visualisation)
# ============================================================================

def plot_cycle_map(prob):
    """
    Generate a 4×4 grid plot of cycle station properties for all operating points.

    Layout:  4 columns (one per operating point) × 4 rows (T, P, MN, Area)
    For each cell: piecewise-linear curve connecting flow stations.

    This function illustrates the thermodynamic state evolution along the gas path
    for each mode, making it easy to see:
      - Temperature rise across fans and combustors
      - Pressure rise/drop across components
      - Mach number changes through ducts
      - Area changes (engine flowpath shape)

    Parameters
    ----------
    prob : om.Problem — the fully-solved OpenMDAO problem instance.
           Values are extracted using prob.get_val().
    """
    import matplotlib.pyplot as plt
    import numpy as np

    # MODES dict: maps operating point name → list of flow station paths + display names.
    # These station paths follow the pattern: '{component}.Fl_O'
    # The pyCycle convention is that every component has a Fl_O (flow outlet) port.
    MODES = {
        'DESIGN_mode2': {
            'label': 'DESIGN 2  —  Fan + AB  (M 2.02 / 60 kft)',
            'stations': ['fc.Fl_O', 'inlet.Fl_O', 'fan1.Fl_O', 'fan2.Fl_O',
                         'ab.Fl_O', 'nozz.Fl_O'],
            'names':    ['Free-\nstream', 'Inlet\nExit', 'Fan1\nExit',
                         'Fan2\nExit', 'AB\nExit', 'Nozzle\nExit'],
        },
        'OD_mode1': {
            'label': 'OD Mode 1  —  Fan Only  (M 0.78 / 35 kft)',
            'stations': ['fc.Fl_O', 'inlet.Fl_O', 'fan1.Fl_O', 'fan2.Fl_O',
                         'ab.Fl_O', 'nozz.Fl_O'],
            'names':    ['Free-\nstream', 'Inlet\nExit', 'Fan1\nExit',
                         'Fan2\nExit', 'Duct\nExit', 'Nozzle\nExit'],
        },
        'OD_mode2': {
            'label': 'OD Mode 2  —  Fan + AB  (M 2.02 / 60 kft)',
            'stations': ['fc.Fl_O', 'inlet.Fl_O', 'fan1.Fl_O', 'fan2.Fl_O',
                         'ab.Fl_O', 'nozz.Fl_O'],
            'names':    ['Free-\nstream', 'Inlet\nExit', 'Fan1\nExit',
                         'Fan2\nExit', 'AB\nExit', 'Nozzle\nExit'],  # 'AB' not 'Duct'
        },
        'OD_mode3': {
            'label': 'OD Mode 3  —  RamJet  (M 3.2 / 80 kft)',
            # Shorter chain: no fan stages
            'stations': ['fc.Fl_O', 'inlet.Fl_O', 'bypass_duct.Fl_O',
                         'combustor.Fl_O', 'nozz.Fl_O'],
            'names':    ['Free-\nstream', 'Inlet\nExit', 'Bypass\nExit',
                         'Combust.\nExit', 'Nozzle\nExit'],
        },
    }

    ROW_LABELS = ['Temperature  (°C)', 'Pressure  (kPa)', 'Mach Number', 'Area  (m²)']

    # Physical plausibility bounds — if a value falls outside these, the solver
    # likely did not converge at that station.  We NaN those points so they
    # do not corrupt the plot scale.
    PHYS = [
        (-273.15,   3226.85),   # Temperature (°C): -273.15 to 3226.85 °C (3500 K upper bound)
        (0.,   3500.),   # Pressure (kPa):   0 to 3500 kPa (well above any expected stagnation)
        (0.,   10.),     # Mach Number:      0 to 10 (comfortably covers hypersonic)
        (0.,   5.),      # Area (m²):        0 to 5 m² (covers large engine flow areas)
    ]

    # Temperature y-axis scaling strategy:
    # T_PER_COLUMN=True → each column (mode) gets its own T y-axis range.
    # This is because fan-only modes (Tmax ~600 K) and combustion modes (Tmax ~2100 K)
    # differ by 3×.  A shared axis would flatten the fan-mode curves to noise.
    T_PER_COLUMN = True

    # -------------------------------------------------------------------------
    # FIRST PASS: collect all data from the solved OpenMDAO problem
    # -------------------------------------------------------------------------
    all_data = {}
    for pt, cfg in MODES.items():
        Tt, Ts, Pt, Ps, MN, A = [], [], [], [], [], []   # lists of per-station values
        fan_power_MW = {}

        for s in cfg['stations']:
            # Full OpenMDAO variable path:  e.g., 'DESIGN.fan1.Fl_O:tot:T'
            base = f'{pt}.{s}'

            def _get(var, units=None, _base=base):
                """
                Safely extract a scalar from the OpenMDAO problem.
                Returns float('nan') if the variable doesn't exist or the query fails
                (e.g., 'nozz.Fl_O:stat:area' may not exist for all nozzle configs).
                """
                try:
                    kw = {'units': units} if units else {}
                    return prob.get_val(f'{_base}:{var}', **kw)[0]  # [0] = first element
                except Exception:
                    return float('nan')   # non-converged or missing → NaN (safe to plot)

            # Collect total (stagnation) temperature and pressure,
            # static temperature and pressure, Mach number, and cross-section area.
            # Units are converted to plotting units (°C, kPa, m²).
            Tt.append(_get('tot:T',     units='degK') - 273.15)    # total temperature  (°C)
            Ts.append(_get('stat:T',    units='degK') - 273.15)    # static temperature (°C)
            Pt.append(_get('tot:P',     units='kPa'))     # total pressure     (kPa)
            Ps.append(_get('stat:P',    units='kPa'))     # static pressure    (kPa)
            MN.append(_get('stat:MN'))                    # Mach number        (dimensionless)
            A.append( _get('stat:area', units='m**2'))    # cross-section area (m²)

        for fan_name in ('fan1', 'fan2'):
            try:
                fan_power_MW[fan_name] = abs(
                    prob.get_val(f'{pt}.{fan_name}.power', units='W')[0]
                ) / 1.0e6
            except Exception:
                fan_power_MW[fan_name] = float('nan')

        # Store as numpy arrays for efficient masking / ylim computation.
        Tt = np.array(Tt)
        Ts = np.array(Ts)
        Pt = np.array(Pt)
        Ps = np.array(Ps)
        MN = np.array(MN)
        A = np.array(A)

        all_data[pt] = dict(
            Tt=Tt,
            Ts=Ts,
            Pt=Pt,
            Ps=Ps,
            MN=MN,
            A=A,
            fan_power_MW=fan_power_MW,
        )

    # -------------------------------------------------------------------------
    # COMPUTE Y-AXIS LIMITS (with 8% padding so lines don't touch the border)
    # -------------------------------------------------------------------------

    def _valid(arrays, lo, hi):
        """Concatenate arrays and return only finite values within (lo, hi)."""
        vals = np.concatenate([a.ravel() for a in arrays])
        return vals[(vals >= lo) & (vals <= hi) & np.isfinite(vals)]

    def _bounds(vals, lo, hi):
        """
        Compute (ymin, ymax) with 8% padding around the data range.
        Falls back to (lo, hi) if no valid data.
        """
        if vals.size == 0:
            return (lo, hi)
        span = vals.max() - vals.min() or 1.   # avoid zero span (flat line)
        pad  = 0.08 * span                     # 8% padding
        return (max(lo, vals.min() - pad), vals.max() + pad)

    # Global ylims for P, MN, Area (rows 1–3): shared across all 4 columns
    # so that cross-mode comparisons are meaningful.
    global_ylims = [None]   # row 0 (T) will be per-column — placeholder here
    for (lo, hi), keys in zip(PHYS[1:], [('Pt', 'Ps'), ('MN',), ('A',)]):
        # Gather all values for this quantity from all modes
        vals = _valid([all_data[pt][k] for pt in MODES for k in keys], lo, hi)
        global_ylims.append(_bounds(vals, lo, hi))

    # Per-column T ylims (row 0): each mode column has its own temperature scale
    col_T_ylims = {}
    for pt in MODES:
        lo, hi = PHYS[0]
        vals = _valid([all_data[pt]['Tt'], all_data[pt]['Ts']], lo, hi)
        col_T_ylims[pt] = _bounds(vals, lo, hi)

    # -------------------------------------------------------------------------
    # SECOND PASS: draw the 4×4 subplot grid
    # -------------------------------------------------------------------------
    fig, axes = plt.subplots(4, 4, figsize=(22, 14))
    fig.suptitle('Duality Engine — Cycle Station Map', fontsize=15, fontweight='bold', y=1.01)

    # ROW_KEYS: each row specifies which data keys to plot, display labels, and colours.
    #   (key_total, key_static, label_total, label_static, colour_total, colour_static)
    # Total quantities: red (#d62728), blue (#1f77b4) — classic matplotlib pair
    # For MN and Area there is no "static" vs "total" concept, so k_stat=None.
    ROW_KEYS = [
        ('Tt', 'Ts', 'T_total', 'T_static', '#d62728', '#1f77b4'),  # row 0: Temperature
        ('Pt', 'Ps', 'P_total', 'P_static', '#d62728', '#1f77b4'),  # row 1: Pressure
        ('MN', None, 'Mach',    None,        '#2ca02c', None),       # row 2: Mach number
        ('A',  None, 'Area',    None,        '#9467bd', None),       # row 3: Area
    ]

    for col, (pt, cfg) in enumerate(MODES.items()):
        # x positions: 0, 1, 2, ... for each flow station
        xs = np.arange(len(cfg['stations']), dtype=float)
        d  = all_data[pt]   # shorthand for this mode's data dict

        for row, (k_tot, k_stat, lbl_tot, lbl_stat, c_tot, c_stat) in enumerate(ROW_KEYS):
            ax     = axes[row, col]
            lo, hi = PHYS[row]   # physical plausibility bounds for masking

            # Select the appropriate y-axis limits
            if row == 0:
                ylim = col_T_ylims[pt]   # per-column for temperature
            else:
                ylim = global_ylims[row] # global for P, MN, area

            def _mask(y, _lo=lo, _hi=hi):
                """
                Replace out-of-bounds / non-finite values with NaN.
                matplotlib will draw a gap in the line at NaN points,
                cleanly showing missing/non-converged data without crashing.
                """
                m = y.copy().astype(float)
                m[~np.isfinite(m) | (m < _lo) | (m > _hi)] = np.nan
                return m

            # Draw a faint vertical dashed line at each station for readability
            for x in xs:
                ax.axvline(x, color='#bbbbbb', linestyle='--', linewidth=0.9, zorder=1)

            # Plot total quantity (solid line with circle markers)
            y_tot = _mask(d[k_tot])
            ax.plot(xs, y_tot, color=c_tot, marker='o', linewidth=2,
                    markersize=5, label=lbl_tot, zorder=3)

            # Plot static quantity if it exists (dashed line with square markers)
            if k_stat is not None:
                y_stat = _mask(d[k_stat])
                ax.plot(xs, y_stat, color=c_stat, marker='s', linewidth=2,
                        markersize=5, linestyle='--', label=lbl_stat, zorder=3)

            # Axis formatting
            ax.set_xticks(xs)
            ax.set_xticklabels(cfg['names'], fontsize=7)   # station names on x-axis
            ax.set_xlim(-0.4, xs[-1] + 0.4)               # small margin around stations
            ax.set_ylim(*ylim)                             # apply computed y-limits
            ax.grid(axis='y', alpha=0.3, zorder=0)        # faint horizontal grid

            # Labels only on outer edges to avoid clutter
            if col == 0:
                ax.set_ylabel(ROW_LABELS[row], fontsize=9)  # y-label only on leftmost column
            if row == 0:
                ax.set_title(cfg['label'], fontsize=8, fontweight='bold', pad=6)  # title on top row
                for fan_name, station_name in (('fan1', 'fan1.Fl_O'), ('fan2', 'fan2.Fl_O')):
                    power_MW = d['fan_power_MW'].get(fan_name, float('nan'))
                    if np.isfinite(power_MW) and station_name in cfg['stations']:
                        station_idx = cfg['stations'].index(station_name)
                        ax.annotate(
                            f'{power_MW:.2f} MW',
                            xy=(xs[station_idx], y_tot[station_idx]),
                            xytext=(0, 12),
                            textcoords='offset points',
                            ha='center',
                            va='bottom',
                            fontsize=7,
                            color='#333333',
                            bbox=dict(
                                boxstyle='round,pad=0.18',
                                fc='white',
                                ec='#999999',
                                alpha=0.82,
                            ),
                            zorder=5,
                        )
            if col == 0 and k_stat is not None:
                ax.legend(fontsize=7, loc='best')  # legend only on left column, rows with two lines

    plt.tight_layout()
    out = 'duality_cycle_map.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'\nCycle map saved -> {out}')


def write_station_geometry_report(prob, out='duality_station_geometry.txt',
                                  inlet_fixed_width=30.0):
    """
    Write a text report of station order and approximate station dimensions.

    The cycle model solves flow station area, but it does not define a physical
    2-D nacelle layout.  For this report, x is the ordered station index along
    the gas path.  Variable inlet stations are represented as fixed-width
    rectangular slots:

        height = area / inlet_fixed_width

    Other stations still use circular-equivalent dimensions:

        D_eq = sqrt(4 * area / pi)
    """
    import math

    modes = {
        'DESIGN_mode2': {
            'label': 'Design Mode 2 - Fan + AB',
            'stations': [
                ('inlet.Fl_O', 'Inlet Exit'),
                ('fan1.Fl_O', 'Fan1 Exit'),
                ('fan2.Fl_O', 'Fan2 Exit'),
                ('ab.Fl_O', 'Afterburner Exit'),
                ('nozz.Throat', 'Nozzle Throat'),
                ('nozz.Fl_O', 'Nozzle Exit'),
            ],
        },
        'OD_mode1': {
            'label': 'Off-Design Mode 1 - Fan Only',
            'stations': [
                ('inlet.Fl_O', 'Inlet Exit'),
                ('fan1.Fl_O', 'Fan1 Exit'),
                ('fan2.Fl_O', 'Fan2 Exit'),
                ('ab.Fl_O', 'Duct Exit'),
                ('nozz.Throat', 'Nozzle Throat'),
                ('nozz.Fl_O', 'Nozzle Exit'),
            ],
        },
        'OD_mode2': {
            'label': 'Off-Design Mode 2 - Fan + AB',
            'stations': [
                ('inlet.Fl_O', 'Inlet Exit'),
                ('fan1.Fl_O', 'Fan1 Exit'),
                ('fan2.Fl_O', 'Fan2 Exit'),
                ('ab.Fl_O', 'Afterburner Exit'),
                ('nozz.Throat', 'Nozzle Throat'),
                ('nozz.Fl_O', 'Nozzle Exit'),
            ],
        },
        'OD_mode3': {
            'label': 'Off-Design Mode 3 - Ramjet',
            'stations': [
                ('inlet.Fl_O', 'Inlet Exit'),
                ('bypass_duct.Fl_O', 'Bypass Duct Exit'),
                ('combustor.Fl_O', 'Combustor Exit'),
                ('nozz.Throat', 'Nozzle Throat'),
                ('nozz.Fl_O', 'Nozzle Exit'),
            ],
        },
    }

    def _area_in2(pt, station):
        try:
            return _scalar(prob, f'{pt}.{station}:stat:area', units='inch**2')
        except Exception:
            return float('nan')

    def _dims_in(area, station):
        diam = math.sqrt(4.0 * area / math.pi) if math.isfinite(area) and area >= 0.0 else float('nan')
        if station == 'inlet.Fl_O' and math.isfinite(area) and inlet_fixed_width > 0.0:
            return area / inlet_fixed_width, inlet_fixed_width, diam
        return diam, diam, diam

    with open(out, 'w') as f:
        print('Duality Engine Station Geometry', file=f)
        print('Units: x = station index from inlet to exhaust; area in inch^2; dimensions in inches', file=f)
        print('Note: pyCycle solves cross-sectional area, not physical height/width.', file=f)
        print(f'      Inlet rows use fixed width = {inlet_fixed_width:.3f} in and height = A/width.', file=f)
        print('      Non-inlet rows use circular-equivalent diameters: D_eq = sqrt(4A/pi).', file=f)

        for pt, cfg in modes.items():
            print('\n' + cfg['label'] + f' ({pt})', file=f)
            print('-' * 98, file=f)
            print(f'{"x":>6}  {"station":<28}  {"path":<26}  {"area":>12}  {"height":>12}  {"width":>12}  {"eq_diam":>12}', file=f)
            print('-' * 98, file=f)

            for idx, (station, name) in enumerate(cfg['stations']):
                area = _area_in2(pt, station)
                height, width, diam = _dims_in(area, station)
                print(f'{idx:6.1f}  {name:<28}  {station:<26}  {area:12.3f}  {height:12.3f}  {width:12.3f}  {diam:12.3f}', file=f)

    print(f'Station geometry report saved -> {out}')


# ============================================================================
# Results viewer (console summary)
# ============================================================================

def viewer(prob, pt, file=sys.stdout):
    """
    Print a concise performance summary for operating point `pt`.

    Extracts key scalars from the solved OpenMDAO problem and formats them
    in a tabular printout, followed by pyCycle's built-in station/compressor/
    burner/nozzle detail tables.

    Parameters
    ----------
    prob : om.Problem — fully-solved problem.
    pt   : str        — operating point name, e.g. 'DESIGN', 'OD_mode3'.
    file : file-like  — output stream (default stdout).
    """

    # -----------------------------------------------------------------------
    # EXTRACT TOP-LEVEL PERFORMANCE SCALARS
    # -----------------------------------------------------------------------
    try:
        mn  = prob.get_val(pt + '.fc.Fl_O:stat:MN')[0]           # freestream Mach number
        alt = prob.get_val(pt + '.fc.alt', units='ft')[0]         # altitude (ft)
        w   = prob.get_val(pt + '.inlet.Fl_O:stat:W')[0]          # inlet mass flow (lbm/s)
        fn  = prob.get_val(pt + '.perf.Fn')[0]                    # net thrust (lbf)
        fg  = prob.get_val(pt + '.perf.Fg')[0]                    # gross thrust (lbf)
        ram = prob.get_val(pt + '.inlet.F_ram')[0]                 # ram drag (lbf)
        opr = prob.get_val(pt + '.perf.OPR')[0]                   # overall pressure ratio
    except Exception as e:
        # If any key variable is missing (e.g., the point didn't converge),
        # print an error and return early rather than crashing.
        print(f'{pt}: could not read perf — {e}', file=file)
        return

    # -----------------------------------------------------------------------
    # TSFC (Thrust-Specific Fuel Consumption, lbm/hr/lbf)
    # -----------------------------------------------------------------------
    # TSFC = Wfuel / Fn — only meaningful when there is a burner.
    # Fan-only mode has no fuel, so TSFC is N/A.
    # Mode detection: check if the point name contains 'mode2' or 'mode3'.
    has_burner = 'mode2' in pt.lower() or 'mode3' in pt.lower()
    if has_burner:
        try:
            tsfc     = prob.get_val(pt + '.perf.TSFC')[0]
            tsfc_str = f'{tsfc:8.5f}'
        except Exception:
            tsfc_str = '     N/A'   # TSFC unavailable (convergence issue?)
    else:
        tsfc_str = '     N/A'   # fan-only: no combustion, no fuel consumption

    # -----------------------------------------------------------------------
    # MODE IDENTIFICATION AND ELECTRIC POWER DRAW
    # -----------------------------------------------------------------------
    if 'mode3' in pt.lower():
        # Ramjet: no fans → no electric power draw to report
        mode_str = 'Mode 3: RamJet'
        pwr_str  = ''
    else:
        # Fan mode: report electric motor power for each fan
        mode_str = 'Mode 2: Fan + AB' if 'mode2' in pt.lower() else 'Mode 1: Fan Only'

        # pyCycle sign convention: compressor power is NEGATIVE (work done ON the gas).
        # abs() converts to a positive number for display purposes.
        # Units: horsepower (hp) for intuitive interpretation.
        p1 = abs(prob.get_val(pt + '.fan1.power', units='hp')[0])  # Fan1 motor power (hp)
        p2 = abs(prob.get_val(pt + '.fan2.power', units='hp')[0])  # Fan2 motor power (hp)
        pwr_str = f'  Fan1={p1:.0f} hp  Fan2={p2:.0f} hp (electric)'

    # -----------------------------------------------------------------------
    # PRINT SUMMARY TABLE
    # -----------------------------------------------------------------------
    print('\n' + '='*70, file=file, flush=True)
    print(f'  {pt}  --  {mode_str}{pwr_str}', file=file, flush=True)
    print('='*70, file=file, flush=True)
    print('  Mach      Alt       W       Fn      Fg    Fram    OPR      TSFC',
          file=file, flush=True)
    print(f' {mn:7.5f}  {alt:7.1f}  {w:7.3f}  {fn:7.1f}  {fg:7.1f}'
          f'  {ram:7.1f}  {opr:6.3f}  {tsfc_str}', file=file, flush=True)

    # -----------------------------------------------------------------------
    # DETAILED FLOW STATION, COMPRESSOR, BURNER, AND NOZZLE TABLES
    # -----------------------------------------------------------------------
    # pyc.print_flow_station: shows Tt, Pt, Ts, Ps, MN, W at each station.
    # pyc.print_compressor:   shows PR, eff, Wc, Nc, R-line from map.
    # pyc.print_burner:       shows FAR, Wfuel, Tt4 for combustors.
    # pyc.print_nozzle:       shows Cv, Fg, throat conditions.

    # Select station list based on mode (Mode 3 has different component names)
    if 'mode3' in pt.lower():
        fs = ['fc.Fl_O', 'inlet.Fl_O', 'bypass_duct.Fl_O', 'combustor.Fl_O', 'nozz.Fl_O']
    else:
        fs = ['fc.Fl_O', 'inlet.Fl_O', 'fan1.Fl_O', 'fan2.Fl_O', 'ab.Fl_O', 'nozz.Fl_O']
    pyc.print_flow_station(prob, [f'{pt}.{s}' for s in fs], file=file)

    # Compressor tables: fan modes only (no fans in ramjet)
    if 'mode3' not in pt.lower():
        pyc.print_compressor(prob, [f'{pt}.fan1', f'{pt}.fan2'], file=file)

    # Burner tables: afterburner in Mode 2, main combustor in Mode 3
    if 'mode2' in pt.lower():
        pyc.print_burner(prob, [f'{pt}.ab'],        file=file)
    if 'mode3' in pt.lower():
        pyc.print_burner(prob, [f'{pt}.combustor'], file=file)

    # Nozzle table: all modes
    pyc.print_nozzle(prob, [f'{pt}.nozz'], file=file)


# ============================================================================
# Standalone nozzle-sizing helpers (Mode 2 and Mode 3)
# ============================================================================
# pyCycle's MPCycle only allows one design=True point.  DESIGN_mode2 and
# DESIGN_mode3 are therefore run as standalone om.Problems to compute their
# nozzle throat areas.  Those areas are then injected into the main MPDuality
# problem via set_val before run_model().
#
# Cruise thrust targets scaled to a PC-24-class 6-passenger jet.
# Approximation: cruise thrust required scales with weight as W/(L/D).
PC24_SCALED_THRUST = {
    'mode1_fan': 1100.0,       # A220-like subsonic transport cruise, L/D ~= 17
    'mode2_turbojet': 2500.0,  # Concorde-like supersonic turbojet cruise, L/D ~= 7.5
    'mode3_ramjet': 3100.0,    # SR-71-like high-Mach cruise, L/D ~= 6.0
}

FAN_RLINE_TARGET = 2.20

CRUISE_CONDITIONS = {
    'mode1': {  # Subsonic fan-only point
        'alt_ft': 35000.0,
        'mach': 0.78,
        'Pt_psia': 5.169,
        'Tt_degR': 441.78,
    },
    'mode2': {  # Concorde-like supersonic cruise
        'alt_ft': 60000.0,
        'mach': 2.02,
        'Pt_psia': 8.396,
        'Tt_degR': 708.22,
    },
    'mode3': {  # SR-71-like high-Mach cruise
        'alt_ft': 80000.0,
        'mach': 3.20,
        'Pt_psia': 19.800,
        'Tt_degR': 1212.68,
    },
}


def _scalar(prob, name, units=None):
    val = prob.get_val(name, units=units) if units else prob.get_val(name)
    return float(val[0])

def _run_design_mode2():
    """
    Standalone DualityFanAB(design=True) sizing run at M=0.8, 30 000 ft.
    Returns the nozzle throat area sized for the PC-24-scaled Concorde-like
    cruise thrust target at T4=3200 R.
    """
    p = om.Problem()
    p.model = DualityFanAB(design=True)
    p.setup()
    p.set_solver_print(level=-1)

    # Shared cycle params that MPCycle normally broadcasts
    p.set_val('inlet.ram_recovery', 0.99)
    p.set_val('nozz.Cv',           0.99)

    # Fan independent design speeds (counter-rotating, same magnitude)
    p.set_val('N_fan1', 6000., units='rpm')
    p.set_val('N_fan2', 6000., units='rpm')

    # Internal station Mach numbers (determine duct areas)
    p.set_val('inlet.MN', 0.60)
    p.set_val('fan1.MN',  0.45)
    p.set_val('fan2.MN',  0.40)
    p.set_val('ab.MN',    0.35)
    p.set_val('ab.dPqP',  0.03)

    # Flight condition: Concorde-like cruise
    p.set_val('fc.alt', CRUISE_CONDITIONS['mode2']['alt_ft'], units='ft')
    p.set_val('fc.MN',  CRUISE_CONDITIONS['mode2']['mach'])

    # Fan design pressure ratios (same as DESIGN_mode1)
    p.set_val('fan1.PR', 1.50)
    p.set_val('fan2.PR', 1.30)

    # Balance targets
    p.set_val('balance.rhs:W',   PC24_SCALED_THRUST['mode2_turbojet'], units='lbf')
    p.set_val('balance.rhs:FAR', 3200., units='degR')   # T4 target

    # Initial guesses consistent with the chosen cruise condition.
    p['balance.W']   = 40.
    p['balance.FAR'] = 0.025
    p['fc.balance.Pt'] = CRUISE_CONDITIONS['mode2']['Pt_psia']
    p['fc.balance.Tt'] = CRUISE_CONDITIONS['mode2']['Tt_degR']
    p.set_val('inlet.Fl_O:tot:T',  CRUISE_CONDITIONS['mode2']['Tt_degR'], units='degR')
    p.set_val('inlet.Fl_O:tot:P',  8.31,  units='lbf/inch**2')
    p.set_val('fan1.Fl_O:tot:T',   796.,  units='degR')
    p.set_val('fan1.Fl_O:tot:P',  12.47,  units='lbf/inch**2')
    p.set_val('fan2.Fl_O:tot:T',   860.,  units='degR')
    p.set_val('fan2.Fl_O:tot:P',  16.21,  units='lbf/inch**2')
    p.set_val('ab.Fl_O:tot:T',    3200.,  units='degR')
    p.set_val('ab.Fl_O:tot:P',    15.72,  units='lbf/inch**2')

    p.run_model()
    result = {
        'nozz': _scalar(p, 'nozz.Throat:stat:area', units='inch**2'),
        'inlet_area': _scalar(p, 'inlet.Fl_O:stat:area', units='inch**2'),
        'fan1_area': _scalar(p, 'fan1.Fl_O:stat:area', units='inch**2'),
        'fan2_area': _scalar(p, 'fan2.Fl_O:stat:area', units='inch**2'),
        'ab_area': _scalar(p, 'ab.Fl_O:stat:area', units='inch**2'),
        'W': _scalar(p, 'balance.W'),
        'FAR': _scalar(p, 'balance.FAR'),
        'Pt': _scalar(p, 'fc.balance.Pt'),
        'Tt': _scalar(p, 'fc.balance.Tt'),
        'inlet_Tt': _scalar(p, 'inlet.Fl_O:tot:T', units='degR'),
        'inlet_Pt': _scalar(p, 'inlet.Fl_O:tot:P', units='lbf/inch**2'),
        'fan1_Tt': _scalar(p, 'fan1.Fl_O:tot:T', units='degR'),
        'fan1_Pt': _scalar(p, 'fan1.Fl_O:tot:P', units='lbf/inch**2'),
        'fan2_Tt': _scalar(p, 'fan2.Fl_O:tot:T', units='degR'),
        'fan2_Pt': _scalar(p, 'fan2.Fl_O:tot:P', units='lbf/inch**2'),
        'ab_Tt': _scalar(p, 'ab.Fl_O:tot:T', units='degR'),
        'ab_Pt': _scalar(p, 'ab.Fl_O:tot:P', units='lbf/inch**2'),
    }
    for sfx in ('s_PR', 's_Wc', 's_eff', 's_Nc'):
        result[f'fan1_{sfx}'] = _scalar(p, f'fan1.{sfx}')
        result[f'fan2_{sfx}'] = _scalar(p, f'fan2.{sfx}')

    fn = _scalar(p, 'perf.Fn')
    print(f"  DESIGN_mode2 sizing: Fn={fn:.1f} lbf  nozzle throat={result['nozz']:.3f} inch²")
    return result


def _initial_total_conditions(alt_ft, mach):
    """Return rough freestream total pressure [psia] and temperature [degR]."""
    gamma = 1.4
    alt_m = alt_ft * 0.3048
    t0 = 288.15
    p0 = 101325.0
    lapse = -0.0065
    r_air = 287.05287
    g = 9.80665
    if alt_m <= 11000.0:
        ts = t0 + lapse * alt_m
        ps = p0 * (ts / t0) ** (-g / (lapse * r_air))
    else:
        ts = 216.65
        p11 = p0 * (ts / t0) ** (-g / (lapse * r_air))
        ps = p11 * math.exp(-g * (alt_m - 11000.0) / (r_air * ts))
    tt = ts * (1.0 + 0.5 * (gamma - 1.0) * mach**2)
    pt = ps * (tt / ts) ** (gamma / (gamma - 1.0))
    return pt / 6894.757293168, tt * 1.8


def _run_design_mode3(alt_ft=None, mach=None, thrust_lbf=None, t4_degR=3800.0):
    """
    Standalone DualityRamjet(design=True) sizing run at M=2.5, 40 000 ft.
    Returns a dict with nozzle throat area, inlet area, bypass_duct area,
    and combustor area (inch²),
    sized for the PC-24-scaled SR-71-like cruise thrust target at T4=3800 R.
    """
    p = om.Problem()
    p.model = DualityRamjet(design=True)
    p.setup()
    p.set_solver_print(level=-1)

    p.set_val('inlet.ram_recovery', 0.99)
    p.set_val('nozz.Cv',           0.99)

    p.set_val('inlet.MN',       0.20)   # subsonic after normal-shock deceleration
    p.set_val('bypass_duct.MN', 0.15)   # low-velocity passage
    p.set_val('combustor.MN',   0.25)   # combustor entrance Mach
    p.set_val('bypass_duct.dPqP', 0.01)
    p.set_val('combustor.dPqP',   0.03)

    alt_ft = CRUISE_CONDITIONS['mode3']['alt_ft'] if alt_ft is None else alt_ft
    mach = CRUISE_CONDITIONS['mode3']['mach'] if mach is None else mach
    thrust_lbf = PC24_SCALED_THRUST['mode3_ramjet'] if thrust_lbf is None else thrust_lbf
    pt_psia, tt_degR = _initial_total_conditions(alt_ft, mach)

    p.set_val('fc.alt', alt_ft, units='ft')
    p.set_val('fc.MN',  mach)

    p.set_val('balance.rhs:W',   thrust_lbf, units='lbf')
    p.set_val('balance.rhs:FAR', t4_degR, units='degR')

    # Initial guesses consistent with the chosen cruise condition.
    p['balance.W']   = 80.
    p['balance.FAR'] = 0.04
    p['fc.balance.Pt'] = pt_psia
    p['fc.balance.Tt'] = tt_degR
    p.set_val('inlet.Fl_O:tot:T',       tt_degR, units='degR')
    p.set_val('inlet.Fl_O:tot:P',       pt_psia, units='lbf/inch**2')
    p.set_val('bypass_duct.Fl_O:tot:T', tt_degR, units='degR')
    p.set_val('bypass_duct.Fl_O:tot:P', 0.99 * pt_psia, units='lbf/inch**2')
    p.set_val('combustor.Fl_O:tot:T',   t4_degR,  units='degR')
    p.set_val('combustor.Fl_O:tot:P',   0.96 * pt_psia, units='lbf/inch**2')

    p.run_model()
    result = {
        'nozz': _scalar(p, 'nozz.Throat:stat:area', units='inch**2'),
        'inlet_area': _scalar(p, 'inlet.Fl_O:stat:area', units='inch**2'),
        'bypass_duct': _scalar(p, 'bypass_duct.Fl_O:stat:area', units='inch**2'),
        'combustor': _scalar(p, 'combustor.Fl_O:stat:area', units='inch**2'),
        'W': _scalar(p, 'balance.W'),
        'FAR': _scalar(p, 'balance.FAR'),
        'Pt': _scalar(p, 'fc.balance.Pt'),
        'Tt': _scalar(p, 'fc.balance.Tt'),
        'inlet_Tt': _scalar(p, 'inlet.Fl_O:tot:T', units='degR'),
        'inlet_Pt': _scalar(p, 'inlet.Fl_O:tot:P', units='lbf/inch**2'),
        'bypass_Tt': _scalar(p, 'bypass_duct.Fl_O:tot:T', units='degR'),
        'bypass_Pt': _scalar(p, 'bypass_duct.Fl_O:tot:P', units='lbf/inch**2'),
        'combustor_Tt': _scalar(p, 'combustor.Fl_O:tot:T', units='degR'),
        'combustor_Pt': _scalar(p, 'combustor.Fl_O:tot:P', units='lbf/inch**2'),
    }
    fn = _scalar(p, 'perf.Fn')
    print(f"  DESIGN_mode3 sizing: Fn={fn:.1f} lbf  nozzle throat={result['nozz']:.3f} inch²"
          f"  inlet={result['inlet_area']:.1f} in²  bypass={result['bypass_duct']:.1f} in²"
          f"  combustor={result['combustor']:.1f} in²")
    return result


# ============================================================================
# Entry point
# ============================================================================

if __name__ == '__main__':
    import time

    # -----------------------------------------------------------------------
    # STEP 1: standalone sizing for Mode 3
    # -----------------------------------------------------------------------
    # Mode 2 is now the primary design=True point inside MPDuality. Only the
    # ramjet geometry still needs a separate standalone sizing run.
    print('\n--- Standalone sizing: DESIGN_mode3 ---')
    d3 = _run_design_mode3()
    nozz_area_mode3    = d3['nozz']
    inlet_area_mode3   = d3['inlet_area']
    bypass_area_mode3  = d3['bypass_duct']
    combust_area_mode3 = d3['combustor']

    # -----------------------------------------------------------------------
    # STEP 2: build and set up the main multi-point problem
    # -----------------------------------------------------------------------
    prob = om.Problem()
    prob.model = mp = MPDuality()   # set the multi-point model as the root system
    prob.setup()                    # triggers setup() on all subsystems recursively

    # -----------------------------------------------------------------------
    # DESIGN_mode2 — Concorde-like supersonic fan + afterburner cruise.
    # This is the primary sizing point for the multi-point model.
    # -----------------------------------------------------------------------
    prob.set_val('DESIGN_mode2.fc.alt',        CRUISE_CONDITIONS['mode2']['alt_ft'], units='ft')
    prob.set_val('DESIGN_mode2.fc.MN',         CRUISE_CONDITIONS['mode2']['mach'])
    prob.set_val('DESIGN_mode2.balance.rhs:W', PC24_SCALED_THRUST['mode2_turbojet'], units='lbf')
    prob.set_val('DESIGN_mode2.balance.rhs:FAR', 3200., units='degR')
    prob.set_val('DESIGN_mode2.fan1.PR',       1.50)
    prob.set_val('DESIGN_mode2.fan2.PR',       1.30)
    prob['DESIGN_mode2.balance.W']       = 35.0
    prob['DESIGN_mode2.balance.FAR']     = 0.035
    prob['DESIGN_mode2.fc.balance.Pt']   = CRUISE_CONDITIONS['mode2']['Pt_psia']
    prob['DESIGN_mode2.fc.balance.Tt']   = CRUISE_CONDITIONS['mode2']['Tt_degR']
    prob.set_val('DESIGN_mode2.inlet.Fl_O:tot:T', CRUISE_CONDITIONS['mode2']['Tt_degR'], units='degR')
    prob.set_val('DESIGN_mode2.inlet.Fl_O:tot:P', 8.31, units='lbf/inch**2')
    prob.set_val('DESIGN_mode2.fan1.Fl_O:tot:T',  796., units='degR')
    prob.set_val('DESIGN_mode2.fan1.Fl_O:tot:P',  12.47, units='lbf/inch**2')
    prob.set_val('DESIGN_mode2.fan2.Fl_O:tot:T',  860., units='degR')
    prob.set_val('DESIGN_mode2.fan2.Fl_O:tot:P',  16.21, units='lbf/inch**2')
    prob.set_val('DESIGN_mode2.ab.Fl_O:tot:T',    3200., units='degR')
    prob.set_val('DESIGN_mode2.ab.Fl_O:tot:P',    15.72, units='lbf/inch**2')

    # -----------------------------------------------------------------------
    # Inject variable-area nozzle throat areas from standalone sizing runs.
    # These are the RHS of the W balance in each OD point — the physical
    # throat constraint.  Must be set AFTER prob.setup() and BEFORE run_model().
    # -----------------------------------------------------------------------
    prob.set_val('OD_mode1.balance.rhs:W', 118.000, units='inch**2')
    prob.set_val('OD_mode3.balance.rhs:W', nozz_area_mode3,  units='inch**2')
    prob.set_val('OD_mode1.balance.rhs:inlet_area', 0.55)
    prob.set_val('OD_mode2.balance.rhs:inlet_area', 0.60)

    # Inject Mode 3 cross-section areas from the standalone Blackbird-like
    # ramjet sizing run.
    prob.set_val('OD_mode3.inlet.area',          inlet_area_mode3,   units='inch**2')
    prob.set_val('OD_mode3.bypass_duct.area', bypass_area_mode3,  units='inch**2')
    prob.set_val('OD_mode3.combustor.area',   combust_area_mode3, units='inch**2')

    # -----------------------------------------------------------------------
    # OD_mode1 — A220-like cruise fan-only check
    # -----------------------------------------------------------------------
    prob['OD_mode1.balance.W']      = 27.
    prob['OD_mode1.balance.inlet_area'] = 260.
    prob.set_val('OD_mode1.N_fan1', 5135., units='rpm')
    prob.set_val('OD_mode1.N_fan2', 4847., units='rpm')
    prob['OD_mode1.fc.balance.Pt']  = CRUISE_CONDITIONS['mode1']['Pt_psia']
    prob['OD_mode1.fc.balance.Tt']  = CRUISE_CONDITIONS['mode1']['Tt_degR']
    prob.set_val('OD_mode1.inlet.Fl_O:tot:T',  CRUISE_CONDITIONS['mode1']['Tt_degR'], units='degR')
    prob.set_val('OD_mode1.inlet.Fl_O:tot:P',  5.12, units='lbf/inch**2')
    prob.set_val('OD_mode1.fan1.Fl_O:tot:T',   500., units='degR')
    prob.set_val('OD_mode1.fan1.Fl_O:tot:P',   7.68, units='lbf/inch**2')
    prob.set_val('OD_mode1.fan2.Fl_O:tot:T',   539., units='degR')
    prob.set_val('OD_mode1.fan2.Fl_O:tot:P',   9.98, units='lbf/inch**2')
    prob.set_val('OD_mode1.ab.Fl_O:tot:T',     539., units='degR')
    prob.set_val('OD_mode1.ab.Fl_O:tot:P',     9.88, units='lbf/inch**2')

    # -----------------------------------------------------------------------
    # OD_mode2 — Concorde-like cruise turbojet / AB check
    # -----------------------------------------------------------------------
    prob['OD_mode2.balance.W']      = 35.
    prob['OD_mode2.balance.inlet_area'] = 260.
    prob['OD_mode2.balance.FAR']    = 0.035
    prob.set_val('OD_mode2.N_fan1', 6000., units='rpm')
    prob.set_val('OD_mode2.N_fan2', 6000., units='rpm')
    prob['OD_mode2.fc.balance.Pt']  = CRUISE_CONDITIONS['mode2']['Pt_psia']
    prob['OD_mode2.fc.balance.Tt']  = CRUISE_CONDITIONS['mode2']['Tt_degR']
    prob.set_val('OD_mode2.inlet.Fl_O:tot:T',  CRUISE_CONDITIONS['mode2']['Tt_degR'], units='degR')
    prob.set_val('OD_mode2.inlet.Fl_O:tot:P',  8.31, units='lbf/inch**2')
    prob.set_val('OD_mode2.fan1.Fl_O:tot:T',   796., units='degR')
    prob.set_val('OD_mode2.fan1.Fl_O:tot:P',   12.47, units='lbf/inch**2')
    prob.set_val('OD_mode2.fan2.Fl_O:tot:T',   860., units='degR')
    prob.set_val('OD_mode2.fan2.Fl_O:tot:P',   16.21, units='lbf/inch**2')
    prob.set_val('OD_mode2.ab.Fl_O:tot:T',     3200., units='degR')
    prob.set_val('OD_mode2.ab.Fl_O:tot:P',     15.72, units='lbf/inch**2')

    # -----------------------------------------------------------------------
    # OD_mode3 — SR-71-like cruise ramjet check
    # -----------------------------------------------------------------------
    prob['OD_mode3.balance.W']     = d3['W']
    prob['OD_mode3.balance.FAR']   = d3['FAR']
    prob['OD_mode3.fc.balance.Pt'] = d3['Pt']
    prob['OD_mode3.fc.balance.Tt'] = d3['Tt']
    prob.set_val('OD_mode3.inlet.Fl_O:tot:T',       d3['inlet_Tt'], units='degR')
    prob.set_val('OD_mode3.inlet.Fl_O:tot:P',       d3['inlet_Pt'], units='lbf/inch**2')
    prob.set_val('OD_mode3.bypass_duct.Fl_O:tot:T', d3['bypass_Tt'], units='degR')
    prob.set_val('OD_mode3.bypass_duct.Fl_O:tot:P', d3['bypass_Pt'], units='lbf/inch**2')
    prob.set_val('OD_mode3.combustor.Fl_O:tot:T',   d3['combustor_Tt'], units='degR')
    prob.set_val('OD_mode3.combustor.Fl_O:tot:P',   d3['combustor_Pt'], units='lbf/inch**2')

    prob.set_solver_print(level=-1)
    prob.set_solver_print(level=2, depth=1)

    t0 = time.time()
    prob.run_model()
    print(f'\nTotal run time: {time.time()-t0:.1f} s')

    for pt in ['DESIGN_mode2', 'OD_mode1', 'OD_mode2', 'OD_mode3']:
        viewer(prob, pt)

    write_station_geometry_report(prob)
    plot_cycle_map(prob)
