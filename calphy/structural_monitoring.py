"""
Structural monitoring functions for detecting phase transformations during simulations.

Contains utilities for:
- Initializing FrameAccumulator with descriptors and detectors
- Pre-training on equilibration trajectories
- Extracting atomic data from LAMMPS
- Processing monitoring frames and extracting statistics
- Generating monitoring plots
"""

import os
import numpy as np
import warnings

# Import structid components
try:
    from structid import FrameAccumulator
    from structid.descriptors import SteinhardtDescriptor
    from structid.detectors import (
        GPRDetector,
        AdaptiveZScoreDetector,
        CumulativeMeanDetector,
    )

    STRUCTID_AVAILABLE = True
except ImportError:
    STRUCTID_AVAILABLE = False


def initialize_accumulator(l_values, cutoff, detector_type="cumulative", threshold=3):
    """
    Initialize a FrameAccumulator with specified descriptor and detector.

    Parameters
    ----------
    l_values : list of int
        Steinhardt parameter orders (e.g., [4, 5, 6, 8, 10])
    cutoff : float
        Neighbor cutoff distance in Angstroms
    detector_type : str, optional
        Type of detector: 'cumulative', 'gpr', 'adaptive' (default: 'cumulative')
    threshold : float, optional
        Detection threshold for cumulative/adaptive detectors (default: 3)

    Returns
    -------
    accumulator : FrameAccumulator or None
        Initialized accumulator, or None if structid not available
    """
    if not STRUCTID_AVAILABLE:
        return None

    descriptor = SteinhardtDescriptor(
        l_values=tuple(l_values),
        cutoff=cutoff,
    )

    if detector_type == "cumulative":
        detector = CumulativeMeanDetector(threshold=threshold)
    elif detector_type == "gpr":
        detector = GPRDetector()
    elif detector_type == "adaptive":
        detector = AdaptiveZScoreDetector(threshold=threshold)
    else:
        raise ValueError(f"Unknown detector type: {detector_type}")

    return FrameAccumulator(descriptor=descriptor, detector=detector)


def setup_pretraining_dump(lmp, simfolder, dump_interval=100, dump_id="_pretrain"):
    """
    Set up LAMMPS dump for equilibration trajectory pre-training.

    Parameters
    ----------
    lmp : lammps object
        LAMMPS instance
    simfolder : str
        Simulation folder path
    dump_interval : int, optional
        Dump every N steps (default: 100)

    Returns
    -------
    pretrain_dump : str
        Path to dump file
    """
    pretrain_dump = os.path.join(simfolder, f"{dump_id}_tmp_pretrain_traj.dump")
    lmp.command(
        f"dump {dump_id} all custom {dump_interval} {pretrain_dump} id type x y z"
    )
    return pretrain_dump, dump_id


def pretrain_accumulator(accumulator, pretrain_dump, logger=None):
    """
    Pre-train accumulator on equilibration trajectory.

    Parameters
    ----------
    accumulator : FrameAccumulator
        The accumulator to pre-train
    pretrain_dump : str
        Path to dump file
    logger : logger object, optional
        Logger for info/warning messages

    Returns
    -------
    success : bool
        True if pre-training succeeded, False otherwise
    """
    if not os.path.exists(pretrain_dump):
        if logger:
            logger.warning(f"Pre-training dump file not found: {pretrain_dump}")
        return False

    try:
        from ase.io import read

        if logger:
            logger.info("Pre-training FrameAccumulator on equilibration trajectory")

        # Read all frames from equilibration
        traj = read(pretrain_dump, format="lammps-dump-text", index=":")

        if logger:
            logger.info(f"Read {len(traj)} frames from equilibration trajectory")

        # Pre-train on all frames
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=Warning, module="sklearn")
            for atoms in traj:
                accumulator.calculate_single_frame(atoms)

        if logger:
            logger.info("FrameAccumulator pre-trained successfully")

        return True

    except Exception as e:
        if logger:
            logger.warning(f"Failed to pre-train FrameAccumulator: {e}")
        return False


def extract_atoms_from_lammps(lmp):
    """
    Extract atomic data from LAMMPS and create ASE Atoms object.

    Parameters
    ----------
    lmp : lammps object
        LAMMPS instance

    Returns
    -------
    atoms : ase.Atoms
        Atoms object with positions, types, and cell
    """
    from ase import Atoms

    # Gather atomic positions and types
    x = lmp.gather_atoms("x", 3)  # positions (natoms × 3)
    atom_types = lmp.gather_atoms("type", 0)  # atom types (natoms,)

    # Get box dimensions
    boxlo, boxhi, xy, yz, xz, periodicity, box_change = lmp.extract_box()

    # Reshape positions to (natoms, 3)
    positions = x.reshape(-1, 3)

    # Create cell from box bounds
    xlo, ylo, zlo = boxlo
    xhi, yhi, zhi = boxhi
    cell = [[xhi - xlo, 0, 0], [xy, yhi - ylo, 0], [xz, yz, zhi - zlo]]

    # Create ASE Atoms object
    atoms = Atoms(numbers=atom_types, positions=positions, cell=cell, pbc=True)

    return atoms


def process_monitoring_frame(accumulator, atoms, logger=None):
    """
    Process a single monitoring frame and extract statistics.

    Parameters
    ----------
    accumulator : FrameAccumulator
        The accumulator to use for analysis
    atoms : ase.Atoms
        Atoms object to analyze
    logger : logger object, optional
        Logger for debug messages

    Returns
    -------
    stats : dict
        Dictionary with keys:
        - 'is_flagged': bool, whether frame is anomalous
        - 'distance': float, descriptor distance
        - 'mean': float, detector mean statistic
        - 'std': float, detector std statistic
    """
    stats = {
        "is_flagged": False,
        "distance": -1.0,
        "mean": -1.0,
        "std": -1.0,
    }

    try:
        # Suppress sklearn convergence warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=Warning, module="sklearn")
            accumulator.calculate_single_frame(atoms)

        # Get results from accumulator
        frames_arr, vectors, distances_arr, flagged = accumulator.get_results()

        # Get the latest distance
        if len(distances_arr) > 0:
            stats["distance"] = distances_arr[-1]

        # Check if current frame is flagged
        current_frame_idx = len(frames_arr) - 1 if len(frames_arr) > 0 else -1
        stats["is_flagged"] = current_frame_idx in flagged

        # Get detector statistics using common interface. Use previous-only distances
        # (exclude current) to match the data the detector saw during detection.
        if hasattr(accumulator.detector, "get_statistics"):
            prev_dists = distances_arr[:-1] if len(distances_arr) > 0 else distances_arr
            try:
                # Try stateless signature first: get_statistics(historical_distances)
                detector_stats = accumulator.detector.get_statistics(prev_dists)
            except TypeError:
                # Fallback: compute mean/std from previous distances directly
                if len(prev_dists) == 0:
                    detector_stats = {"mean": None, "std": None}
                else:
                    finite = prev_dists[np.isfinite(prev_dists)]
                    if len(finite) == 0:
                        detector_stats = {"mean": None, "std": None}
                    else:
                        detector_stats = {
                            "mean": float(np.mean(finite)),
                            "std": float(np.std(finite)),
                        }

            stats["mean"] = (
                detector_stats.get("mean", -1.0)
                if detector_stats.get("mean") is not None
                else -1.0
            )
            stats["std"] = (
                detector_stats.get("std", -1.0)
                if detector_stats.get("std") is not None
                else -1.0
            )

    except Exception as e:
        if logger:
            logger.debug(f"Could not extract detector statistics: {e}")

    return stats


def generate_monitoring_plot(monitor_file, output_file, iteration=1):
    """
    Generate monitoring plot from data file.

    Parameters
    ----------
    monitor_file : str
        Path to monitoring data file
    output_file : str
        Path to output PNG file
    iteration : int, optional
        Iteration number for title (default: 1)

    Returns
    -------
    success : bool
        True if plot generated successfully
    """
    try:
        import matplotlib

        matplotlib.use("Agg")  # Non-interactive backend
        import matplotlib.pyplot as plt

        # Read monitoring data
        data = np.loadtxt(monitor_file)
        if len(data.shape) == 1:  # Single row
            data = data.reshape(1, -1)

        block_idx_arr = data[:, 0]
        step_arr = data[:, 1]
        lambda_arr = data[:, 2]
        temp_arr = data[:, 3]
        is_flagged_arr = data[:, 4]
        distance_arr = data[:, 5]

        # Use temperature as x-axis for plotting (preserve frame order)
        x_arr = temp_arr

        # Use logged mean/std columns if present (they reflect detector-facing statistics)
        if data.shape[1] >= 8:
            mean_arr = data[:, 6]
            std_arr = data[:, 7]
        else:
            mean_arr = np.array([np.mean(distance_arr[: i + 1]) for i in range(len(distance_arr))])
            std_arr = np.array([np.std(distance_arr[: i + 1]) for i in range(len(distance_arr))])

        # Create plot
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(
            x_arr,
            distance_arr,
            "o-",
            label="Descriptor Distance",
            linewidth=2,
            markersize=6,
        )
        ax.plot(x_arr, mean_arr, "k-", label="Detector Mean", linewidth=2)
        ax.plot(
            x_arr,
            mean_arr + 2.0 * std_arr,
            "--",
            label="Mean + 2σ",
            linewidth=1.5,
        )
        ax.plot(
            x_arr,
            mean_arr + 2.5 * std_arr,
            "--",
            label="Mean + 2.5σ",
            linewidth=1.5,
        )
        ax.plot(
            x_arr,
            mean_arr + 3.0 * std_arr,
            "--",
            label="Mean + 3σ",
            linewidth=1.5,
        )

        # Highlight flagged points (anomalies)
        flagged_indices = np.where(is_flagged_arr > 0.5)[0]
        if len(flagged_indices) > 0:
            ax.scatter(
                x_arr[flagged_indices],
                distance_arr[flagged_indices],
                c="red",
                s=200,
                marker="X",
                edgecolors="darkred",
                linewidth=2,
                label="Flagged Anomalies",
                zorder=5,
            )

        ax.set_xlabel("Temperature (K)", fontsize=12)
        ax.set_ylabel("Descriptor Distance", fontsize=12)

        # Detect sweep direction from filename
        sweep_type = "Forward" if "forward" in monitor_file.lower() else "Backward"
        ax.set_title(
            f"Structural Monitoring - {sweep_type} Sweep (Iteration {iteration})",
            fontsize=14,
        )
        ax.legend(loc="best", fontsize=10)
        ax.grid(True, alpha=0.3)

        # Save plot
        plt.tight_layout()
        plt.savefig(output_file, dpi=150, bbox_inches="tight")
        plt.close()

        return True

    except Exception as e:
        return False
