
from feos.pcsaft import PcSaftRecord,PcSaftParameters,PureRecord,Identifier
from feos.eos import EquationOfState,State,PhaseEquilibrium,PhaseDiagram,Contributions
from feos.si import KELVIN,KILO,MOL,METER
import matplotlib.pyplot as plt


import matplotlib.pyplot as plt
import numpy as np
#params = PcSaftParameters.from_json(['hexane', 'octane'], 'gross2001.json')
identifier = Identifier(
    cas='106-97-8',
    name='butane',
    iupac_name='butane',
    smiles='CCCC',
    inchi='InChI=1/C4H10/c1-3-4-2/h3-4H2,1-2H3',
    formula='C4H10')
identifier = Identifier(
    cas='106-97-8',
    name='butane2',
    iupac_name='butane2',
    smiles='CCCC',
    inchi='InChI=1/C4H10/c1-3-4-2/h3-4H2,1-2H3',
    formula='C4H11')
psr = PcSaftRecord(m=2.3316, sigma=3.7086, epsilon_k=222.88)
psr2 = PcSaftRecord(m=2.2316, sigma=3.8086, epsilon_k=202.88)
butane = PureRecord(identifier, molarweight=58.123, model_record=psr)
butane2 = PureRecord(identifier, molarweight=38.123, model_record=psr2)
#parameters = PcSaftParameters.new_pure(butane)
parameters = PcSaftParameters.new_binary([butane, butane2], binary_record=0.015)

eos = EquationOfState.pcsaft(parameters)
#state_nvt = State(eos,temperature=300.15*KELVIN, density=1.25182*KILO*MOL/METER**3, total_moles=100.0*MOL)
#print(state_nvt.dln_phi_dnj())
# Define thermodynamic conditions
#critical_point = State.critical_point(eos)

# Compute properties
#p = critical_point.pressure()
#t = critical_point.temperature
#print(f'Critical point for butane: T={t}, p={p}.')
m=PhaseDiagram.binary_vle(eos,temperature_or_pressure=300.15*KELVIN).to_dict()

plt.plot(m["density vapor"],m["density liquid"])
plt.show()