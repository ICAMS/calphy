.. pytint documentation master file, created by
   sphinx-quickstart on Fri Jan 15 14:33:42 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

calphy
======

``calphy``\ (pronounced *cal-phee*) is a Python library and command line
tool for calculation of free energies. It uses
`LAMMPS <https://www.lammps.org/>`__ as the molecular dynamics driver to
enable calculation of free energies using thermodynamic integration in a
completely automated and efficient manner. The various methods that
calphy can perform are:

-  :math:`F(V_i,T_i)` and :math:`G(P_i,T_i)` for both solid and liquid
   phases at the given thermodynamic conditions using `non-equilibrium
   Hamiltonian
   interpolation <https://linkinghub.elsevier.com/retrieve/pii/S0927025615007089>`__.
-  :math:`F(T_i \to T_f)` and :math:`G(T_i \to T_f)`, temperature
   dependence of Gibbs and Helmholtz free energies using the `reversible
   scaling <https://link.aps.org/doi/10.1103/PhysRevLett.83.3973>`__
   approach.
-  Calculation of solid-solid or solid-liquid phase transition
   temperatures.
-  Calculation of coexistence lines using `dynamic Clausius-Clapeyron
   integration <http://aip.scitation.org/doi/10.1063/1.1420486>`__.
-  Calculation of specific heat :math:`c_P(T)` as a function of
   temperature.
-  Calculation of :math:`F(x, T)` and :math:`G(x, T)` for multicomponent
   systems using `alchemical
   transformations <https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.103801>`__.
-  `Upsampling
   calculations <https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.103801>`__.

Calphy works with all `interatomic potentials implemented in
LAMMPS <https://docs.lammps.org/pairs.html>`__ and can obtain reliable
results with low error bars in as low as 50 ps of simulation time with
less than 3000 atoms. More information about the methods in calphy can
be found in the `associated
publication <https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.5.103801>`__.


Documentation
~~~~~~~~~~~~~

.. toctree::
   gettingstarted 
   documentation

Examples
~~~~~~~~

.. toctree::
   :maxdepth: 2

   examples

API reference
~~~~~~~~~~~~~

.. toctree::
   calphy

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. math::