"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021-2026 (c) Sarath Menon, Yury Lysogorskiy, Ralf Drautz
Interdisciplinary Centre for Advanced Materials Simulation (ICAMS),
Ruhr University Bochum, 44801 Bochum, Germany

calphy is published and distributed under the Academic Software Licence v1.0 (ASL).
calphy is distributed in the hope that it will be useful for non-commercial academic
research, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the LICENSE file for details.

The ASL permits academic non-commercial use only. Contact
sarath.menon@ruhr-uni-bochum.de to enquire about commercial use rights.

More information about the program can be found in:
Menon, Sarath, Yury Lysogorskiy, Jutta Rogal, and Ralf Drautz.
"Automated Free Energy Calculation from Atomistic Simulations." Physical Review Materials 5(10), 2021
DOI: 10.1103/PhysRevMaterials.5.103801

For more information contact:
sarath.menon@ruhr-uni-bochum.de
"""
"""Dynamic Clausius-Clapeyron integration (``mode: dcci``).

Traces a solid-liquid coexistence line P_coex(T) from one known coexistence
point in a single nonequilibrium run per direction, following de Koning,
Antonelli and Yip, J. Chem. Phys. 115, 11025 (2001).  Both phases are driven
side by side at the fixed kinetic temperature T0 of the starting point with the
scaled Hamiltonian K + lambda*U; the scaled pressure P_RS is ramped linearly and
lambda follows from the Clausius-Clapeyron condition in its scaled form,

    d lambda / d P_RS = -(v_s - v_l) / (u_s - u_l),

evaluated from block averages of the potential energy u = <pe>/(N lambda) and
volume v = <V>/N per atom of the two cells.  Every point maps back to the real
system through T = T0/lambda and P = P_RS/lambda.  See DCCI_PLAN.md.
"""
import copy
import os

import numpy as np
import yaml

import calphy.helpers as ph


class DynamicCCI:
    """
    Coordinator of a dynamic Clausius-Clapeyron integration.

    Owns one :class:`~calphy.solid.Solid` and one :class:`~calphy.liquid.Liquid`
    sub-job, equilibrated independently at the starting coexistence point and
    then stepped in lockstep through the coupled sweep.

    Parameters
    ----------
    calculation : Calculation
        The validated ``mode: dcci`` calculation.
    simfolder : str
        Folder for the coupled-sweep output; the cells live in ``solid/`` and
        ``liquid/`` below it.
    log_to_screen : bool
        Also log to the screen.
    """

    def __init__(self, calculation=None, simfolder=None, log_to_screen=False):
        self.calc = copy.deepcopy(calculation)
        self.simfolder = simfolder
        self.log_to_screen = log_to_screen
        self.publications = ["10.1063/1.1420486"]

        if self.calc.md.seed is None:
            self.calc.md.seed = int(np.random.SeedSequence().entropy % (2**31 - 2)) + 1
        np.random.seed(self.calc.md.seed)

        with open(os.path.join(simfolder, "input_file.yaml"), "w") as fout:
            yaml.safe_dump({"calculations": [self.calc.model_dump()]}, fout)

        logfile = os.path.join(self.simfolder, "calphy.log")
        self.logger = ph.prepare_log(logfile, screen=log_to_screen)

        self.t0 = float(self.calc._temperature)
        self.p_start = float(self.calc._pressure)
        self.p_stop = float(self.calc._pressure_stop)
        self.logger.info(
            "Dynamic Clausius-Clapeyron integration from the coexistence point "
            "T = %.2f K, P = %.2f bar to P = %.2f bar" % (self.t0, self.p_start, self.p_stop)
        )
        self.logger.info("Master random seed is %d (md.seed)" % self.calc.md.seed)

    def calculate_coexistence_line(self):
        """Run the whole mode: equilibrate both cells, sweep, integrate, report."""
        raise NotImplementedError("mode dcci: the coupled sweep is not implemented yet")
