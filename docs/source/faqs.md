# FAQs

**I get a `BrokenPipe` error `.../site-packages/pylammpsmpi/utils/lammps.py", line 57, in _send self._process.stdin.flush() BrokenPipeError: [Errno 32] Broken pipe`**

Calphy uses pylammpsmpi which calls the LAMMPS library. pylammpsmpi works only with the OpenMPI version of LAMMPS. Compiling LAMMPS with OpenMPI would solve this issue.

