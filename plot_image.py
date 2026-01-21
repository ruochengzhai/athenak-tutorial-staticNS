# AthPlot: manipulate and plot AthenaK binary data
# Copyright (C) 2025, Jacob Fields <jmf6719@psu.ed>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import os
SCRIPTDIR = os.path.dirname(os.path.realpath(__file__))

import sys
sys.path.append(os.path.abspath(os.path.join(SCRIPTDIR, os.pardir, os.pardir, os.pardir, "plot-tools")))
from athplot.load_ath_bin import BinaryData
from athplot.image import Image

extent = [-20.0, 20.0, -20.0, 20.0]

image = Image("../bin/tov.mhd_w_bcc.00001.bin", extent, "dens", slice_loc=['z', 0.0])

data = image.make_image_data(256, 256)

pcm = plt.imshow(data, cmap='inferno',
                 norm=LogNorm(vmin=1.28e-7, vmax=1.28e-2),
                 interpolation='nearest', origin='lower', extent=image.extent)

plt.xlim(extent[0:2])
plt.ylim(extent[2:])
plt.colorbar(pcm, extend='both', label=r'$\rho$')

plt.show()
