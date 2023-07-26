

from feos.si import KELVIN,BAR
from feos.eos import EquationOfState
from feos.pcsaft import Identifier, PcSaftRecord, PureRecord,PcSaftParameters
from feos.eos import PhaseDiagram
import matplotlib.pyplot as plt
psr1 = PcSaftRecord(m=2.3316, sigma=3.7086, epsilon_k=222.88)
psr2 = PcSaftRecord(m=2.3016, sigma=3.5086, epsilon_k=202.88)
butane = PureRecord(Identifier(""), molarweight=58.123, model_record=psr1)
something = PureRecord(Identifier(""), molarweight=68.123, model_record=psr2)
parameters=PcSaftParameters.new_binary([butane,something])
saft = EquationOfState.pcsaft(parameters)
# dia_p = PhaseDiagram.binary_vle(saft, 300*KELVIN)
# dia_t = PhaseDiagram.binary_vle(saft, BAR)
# f, ax = plt.subplots(1,3,figsize=(20,5))
# ax[0].plot(dia_p.liquid.molefracs[:,0], dia_p.liquid.pressure/BAR, '-k')
# ax[0].plot(dia_p.vapor.molefracs[:,0], dia_p.vapor.pressure/BAR, '-k')
# ax[0].set_xlim(0,1)
# ax[0].set_xlabel('$x_1,y_1$')
# ax[0].set_ylabel('$p$ / bar')

# ax[1].plot(dia_t.liquid.molefracs[:,0], dia_t.liquid.temperature/KELVIN, '-g')
# ax[1].plot(dia_t.vapor.molefracs[:,0], dia_t.vapor.temperature/KELVIN, '-g')
# ax[1].set_xlim(0,1)
# ax[1].set_xlabel('$x_1,y_1$')
# ax[1].set_ylabel('$T$ / K')

# ax[2].plot([0,1], [0,1], '--k')
# ax[2].plot(dia_t.liquid.molefracs[:,0], dia_t.vapor.molefracs[:,0], '-g')
# ax[2].plot(dia_p.liquid.molefracs[:,0], dia_p.vapor.molefracs[:,0], '-k')
# ax[2].set_xlim(0,1)
# ax[2].set_ylim(0,1)
# ax[2].set_xlabel('$x_1$')
# ax[2].set_ylabel('$y_1$');
# plt.show()