Input file
----------

``pytint`` uses a ``yaml`` file for specifying the input options. In
this section, the various blocks of the input file is discussed. A
complete sample input file is also provided in the Examples section.

``main`` block
~~~~~~~~~~~~~~

``main`` block consists of the major options that the user has to
provided. A sample block is shown below.

::

    main:
      temperature: [1000, 1400]
      pressure: [0]
      element: 'Cu'
      lattice: [FCC, LQD]
      nsims: 3

-  ``temperature``: List of temperatures at which NEHI calculations have
   to be done. In case of regular use, one calculation will be started
   for each temperature. In case of ``--mode rs``, one NEHI calculation
   will be done for the first temperature. After that, reversible
   scaling calculation will be done to extend the free energy upto the
   last temperature specified.

-  ``pressure``: Pressure for the calculation. Currently only zero
   pressure (Helmholtz free energy) is supported.

-  ``element``: The chemical symbol of the element used.

-  ``lattice``: The lattices for which the free energy calculations have
   to be done. Supported lattices are ``BCC``, ``FCC``, ``HCP``,
   ``DIA``, ``SC`` and ``LQD``.

-  ``nsims``: The number of independent calculations to be carried out
   for finding the error in the estimated free energy.

``md`` block
~~~~~~~~~~~~

::

    md:
      timestep: 0.001
      pair_style: pace
      pair_coeff: "* * Cu.ace Cu"
      mass: 63.546
      tdamp: 0.1
      pdamp: 0.1
      nx: 5
      ny: 5
      nz: 5
      te: 25000
      ts: 50000

-  ``timestep``: Timestep in ps for the md simulations
-  ``pair_style``: Pair style used in LAMMPS. Supported pair styles are
   ``pace``, ``eam`` and its variants, ``sw`` and ``snap``.
-  ``pair_coeff``: Pair coefficient command used in LAMMPS. The
   relative/full path of the input potential file is also specified
   here.
-  ``mass``: Atomic mass of the element.
-  ``tdamp``: Thermostat damping in units of time.
-  ``pdamp``: Barostat damping in units of time.
-  ``nx``: Number of units cells in the 001 direction.
-  ``ny``: Number of units cells in the 010 direction.
-  ``nz``: Number of units cells in the 100 direction.
-  ``te``: Number of time steps for equilibration runs.
-  ``ts``: Number of time steps for switching runs.

``queue`` block
~~~~~~~~~~~~~~~

This block specifies the input parameters for submitting the job on a
cluster

::

    queue:
      scheduler: slurm
      cores: 40
      jobname: cu
      walltime: "23:50:00"
      queuename: shorttime
      memory: 3GB
      modules:
        - anaconda/4
      commands:
        - source .bashrc
        - conda activate py3
      #any other extra options
      #options:
      # - "-j Y"

-  ``scheduler``: The scheduler to be used for calculations. Supported
   options are ``local``, ``slurm`` or ``sge``.
-  ``cores``: Number of cores to be used for the md runs.
-  ``jobname``: Name of the job. Ignored for ``local``.
-  ``walltime``: Walltime for the job. Ignored for ``local``.
-  ``queuename``: Name of the submission queue. Ignored for ``local``.
-  ``memory``: Total memory requested per core. Ignored for ``local``.
-  ``modules``: Name of module that need to be loaded.
-  ``commands``: Extra commands that will be run in the beginning of the
   submission script. If a conda environment is used, the activate
   statements will be here.
-  ``options``: Further special options for the submission script.

