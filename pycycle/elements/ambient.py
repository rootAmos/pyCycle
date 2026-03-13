import openmdao.api as om

from pycycle.elements.US1976 import USatm1976Comp

class DeltaTs(om.ExplicitComponent):
    """Computes temperature based on delta from atmospheric"""

    def setup(self):

        # inputs
        self.add_input('Ts_in', val=500.0, units='degR', desc='Temperature from atmospheric model')
        self.add_input('dTs', val=0.0, units='degR', desc='Delta from standard day temperature')

        self.add_output('Ts', shape=1, units='degR', desc='Temperature with delta')

        self.declare_partials('Ts', ['Ts_in', 'dTs'], val=1.0)

    def compute(self, inputs, outputs):
        outputs['Ts'] = inputs['Ts_in'] + inputs['dTs']

    def compute_partials(self, inputs, partials):
        pass


class Ambient(om.Group):
    """Determines pressure, temperature and density base on altitude from an input standard atmosphere table"""

    def setup(self):
        readAtm = self.add_subsystem('readAtmTable', USatm1976Comp(), promotes=('alt', 'Ps', 'rhos'))

        self.add_subsystem('dTs', DeltaTs(), promotes=('dTs', 'Ts'))
        self.connect('readAtmTable.Ts', 'dTs.Ts_in')

        # self.set_order(['readAtmTable','dTs'])


if __name__ == "__main__":
    """Run Ambient standalone. IVC sets alt and dTs; print Ps, Ts, rhos."""

    prob = om.Problem()
    model = prob.model = om.Group()
    ivc = model.add_subsystem('ivc', om.IndepVarComp(), promotes_outputs=['*'])
    ivc.add_output('alt', 30000.0, units='ft')
    ivc.add_output('dTs', 0.0, units='degR')

    model.add_subsystem('ambient', Ambient())
    model.connect('alt', 'ambient.alt')
    model.connect('dTs', 'ambient.dTs')

    prob.setup()
    prob.run_model()

    print("--- Ambient standalone test ---")
    print("alt (ft) =", prob.get_val('alt')[0])
    print("Ps (psi) =", prob.get_val('ambient.Ps')[0])
    print("Ts (degR) =", prob.get_val('ambient.Ts')[0])
    print("rhos (slug/ft**3) =", prob.get_val('ambient.rhos')[0])
