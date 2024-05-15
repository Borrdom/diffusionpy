from __future__ import division, print_function, absolute_import
import numpy as np
from .diffusion import Diffusion_MS,Gammaij
from .surface import time_dep_surface
from .PCSAFT import lngi,vpure,NETVLE,supersaturation
from .utilities import DasDennis
from .xloil_functions import Diffusion_MS_iter_xloil, Diffusion_MS_xloil
