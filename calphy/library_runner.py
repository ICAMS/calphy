"""
calphy: a Python library and command line interface for automated free
energy calculations.

Copyright 2021  (c) Sarath Menon^1, Yury Lysogorskiy^2, Ralf Drautz^2
^1: Max Planck Institut fuer Eisenforschung, Dusseldorf, Germany
^2: Ruhr-University Bochum, Bochum, Germany

calphy is published and distributed under the Academic Software License v1.0 (ASL).
See the LICENSE FILE for more details.
"""
"""pylammpsmpi-backed LAMMPS runner (the optional ``library`` execution mode).

``LibraryRunner`` drives a live in-memory LAMMPS session through
:mod:`pylammpsmpi` instead of segmenting scripts for an external binary.
Commands are forwarded immediately, so :meth:`sync` is a no-op -- there is no
segmentation, no restart replay, and no ``.seg*`` files.  Data accessors
currently read the same files the executable backend does; because every read
goes through the runner they can later be swapped for direct library calls
without touching driver code.

pylammpsmpi is an optional dependency (``pip install calphy[library]``) and is
imported lazily inside :class:`LibraryRunner` -- importing this module (or any
other part of calphy) never requires it.
"""
import os

from calphy.runner import BaseRunner, _normalize_cmdargs


class LibraryRunner(BaseRunner):
    """Drives a live ``pylammpsmpi.LammpsLibrary`` session.

    The raw library handle is exposed as :attr:`lmp` so future backend-specific
    accessors (direct thermo/structure reads) can build on it.
    """

    def __init__(self, *, cores, cmdargs, directory):
        super().__init__(directory)
        try:
            from pylammpsmpi import LammpsLibrary
        except ImportError as exc:
            raise ImportError(
                "execution_mode 'library' requires pylammpsmpi, which is an "
                "optional dependency. Install it with `pip install "
                "calphy[library]` (the `lammps` python module must also be "
                "importable, e.g. from the conda-forge lammps package), or "
                "remove execution_mode from the input file to use the default "
                "executable mode."
            ) from exc

        self.cores = cores
        self._log_index = 0
        self._closed = False
        cmdargs = _normalize_cmdargs(cmdargs)
        # Suppress LAMMPS stdout; Python logging handles screen output.
        if "-screen" not in cmdargs:
            cmdargs.extend(["-screen", "none"])
        # Redirect the session log to a known scratch name so rotate_logs can
        # collect per-stage logs (LAMMPS would otherwise write log.lammps).
        if "-log" not in cmdargs:
            cmdargs.extend(["-log", self._log_name(0)])
        self.cmdargs = cmdargs

        self.lmp = LammpsLibrary(
            cores=cores, working_directory=directory, cmdargs=cmdargs
        )

    @staticmethod
    def _log_name(k):
        return "calphy.live.%d.log" % k

    # -- backend contract ----------------------------------------------------- #
    def _dispatch(self, cmd):
        """Forward the validated command to the live session immediately."""
        self.lmp.command(cmd)

    def sync(self):
        """No-op: the session is live and LAMMPS flushes fix output as it runs."""

    def close(self):
        """End the pylammpsmpi session (flushes and closes the LAMMPS log)."""
        if not self._closed:
            self.lmp.close()
            self._closed = True

    def rotate_logs(self, stage_name):
        """Move the log written since the last rotation to ``<stage_name>.log.lammps``.

        On a live session, switching the LAMMPS ``log`` target closes the current
        scratch log so it can be renamed; the drivers however always rotate right
        after :meth:`close`, where the scratch log is already complete and no
        command can (or need) be sent.  The ``log`` command is internal
        bookkeeping, not part of the logical command stream, so it bypasses
        :meth:`command`.
        """
        self._log_index += 1
        if not self._closed:
            self.lmp.command("log %s" % self._log_name(self._log_index))
        prev = os.path.join(self.directory, self._log_name(self._log_index - 1))
        out = os.path.join(self.directory, "%s.log.lammps" % stage_name)
        if os.path.exists(prev):
            os.replace(prev, out)
        else:
            open(out, "w").close()
