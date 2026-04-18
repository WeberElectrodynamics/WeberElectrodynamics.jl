"""Rolling circular buffer for per-frame simulation snapshots.

Numpy port of ``RollingBuffer`` in ``ext/WeberElectrodynamicsMakieExt.jl``.
All storage is preallocated at construction; ``push_step`` advances a modular
cursor. ``linearize`` returns the last ``tail`` entries in chronological order.
"""

from __future__ import annotations

import numpy as np

__all__ = ["RollingBuffer"]


class RollingBuffer:
    """Fixed-capacity circular buffer of per-frame state.

    Arrays are shaped ``(particle, dim, slot)`` or ``(slot,)``; ``slot`` is the
    fastest-varying index, matching the column-major Julia layout that the
    Makie extension uses.
    """

    def __init__(self, capacity: int, n_particles: int, dims: int):
        if capacity <= 0:
            raise ValueError("capacity must be positive")
        if n_particles < 1 or dims < 1:
            raise ValueError("n_particles and dims must be positive")

        self.capacity = capacity
        self.n = n_particles
        self.dims = dims
        self.count = 0
        self.cursor = 0

        self.t = np.empty(capacity, dtype=np.float64)
        self.positions = np.empty((n_particles, dims, capacity), dtype=np.float64)
        self.total_energy = np.empty(capacity, dtype=np.float64)
        self.particle_q = np.empty((n_particles, dims, capacity), dtype=np.float64)
        self.particle_p = np.empty((n_particles, dims, capacity), dtype=np.float64)

        self.pair_separation: dict[tuple[int, int], np.ndarray] = {}
        self.pair_radial_velocity: dict[tuple[int, int], np.ndarray] = {}
        for i in range(1, n_particles + 1):
            for j in range(i + 1, n_particles + 1):
                self.pair_separation[(i, j)] = np.empty(capacity, dtype=np.float64)
                self.pair_radial_velocity[(i, j)] = np.empty(capacity, dtype=np.float64)

    def reset(self) -> None:
        self.count = 0
        self.cursor = 0

    def push_step(
        self,
        t: float,
        q: np.ndarray,
        p: np.ndarray,
        masses: np.ndarray,
        energy: float,
    ) -> None:
        """Record a single frame at the current cursor and advance."""
        idx = self.cursor
        dims = self.dims

        self.t[idx] = t
        self.total_energy[idx] = energy

        q2 = q.reshape(self.n, dims)
        p2 = p.reshape(self.n, dims)
        self.positions[:, :, idx] = q2
        self.particle_q[:, :, idx] = q2
        self.particle_p[:, :, idx] = p2

        for (i, j), sep in self.pair_separation.items():
            # 1-indexed pair keys match the Julia convention
            dq = q2[j - 1] - q2[i - 1]
            dv = p2[j - 1] / masses[j - 1] - p2[i - 1] / masses[i - 1]
            r = float(np.linalg.norm(dq))
            sep[idx] = r
            self.pair_radial_velocity[(i, j)][idx] = float(np.dot(dq, dv) / r) if r > 0.0 else 0.0

        self.cursor = (idx + 1) % self.capacity
        self.count = min(self.count + 1, self.capacity)

    def linearize_1d(self, arr: np.ndarray, tail: int | None = None) -> np.ndarray:
        """Return last ``tail`` entries of a 1D array in chronological order."""
        if self.count == 0:
            return np.empty(0, dtype=arr.dtype)
        cap, cur, cnt = self.capacity, self.cursor, self.count
        if cnt < cap:
            out = arr[:cnt]
        else:
            out = np.concatenate([arr[cur:], arr[:cur]])
        if tail is not None and tail < out.size:
            out = out[-tail:]
        return out

    def linearize_last_axis(self, arr: np.ndarray, tail: int | None = None) -> np.ndarray:
        """Linearize a ND array whose last axis is the circular-buffer slot."""
        if self.count == 0:
            shape = arr.shape[:-1] + (0,)
            return np.empty(shape, dtype=arr.dtype)
        cap, cur, cnt = self.capacity, self.cursor, self.count
        if cnt < cap:
            out = arr[..., :cnt]
        else:
            out = np.concatenate([arr[..., cur:], arr[..., :cur]], axis=-1)
        if tail is not None and tail < out.shape[-1]:
            out = out[..., -tail:]
        return out

    def latest(self, arr: np.ndarray) -> np.ndarray:
        """Return the most recently pushed entry along the last axis."""
        if self.count == 0:
            raise IndexError("buffer is empty")
        last_idx = (self.cursor - 1) % self.capacity
        return arr[..., last_idx]
