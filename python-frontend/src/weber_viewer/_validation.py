"""Shared validation helpers for public Python entry points."""

from __future__ import annotations

import math
import operator
from typing import Sequence

import numpy as np


def as_int(
    name: str,
    value: object,
    *,
    min_value: int | None = None,
    max_value: int | None = None,
) -> int:
    """Return ``value`` as an integer, rejecting bools and floats."""
    if isinstance(value, bool):
        raise TypeError(f"{name} must be an integer")
    try:
        ivalue = operator.index(value)
    except TypeError as exc:
        raise TypeError(f"{name} must be an integer") from exc
    if min_value is not None and ivalue < min_value:
        raise ValueError(f"{name} must be >= {min_value} (got {ivalue})")
    if max_value is not None and ivalue > max_value:
        raise ValueError(f"{name} must be <= {max_value} (got {ivalue})")
    return ivalue


def as_finite_float(name: str, value: object) -> float:
    """Return ``value`` as a finite float."""
    if isinstance(value, bool):
        raise TypeError(f"{name} must be a finite number")
    try:
        fvalue = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a finite number") from exc
    if not math.isfinite(fvalue):
        raise ValueError(f"{name} must be finite (got {value!r})")
    return fvalue


def as_positive_float(name: str, value: object) -> float:
    """Return ``value`` as a finite positive float."""
    fvalue = as_finite_float(name, value)
    if fvalue <= 0.0:
        raise ValueError(f"{name} must be positive (got {fvalue})")
    return fvalue


def as_flat_float_array(
    name: str,
    values: Sequence[float],
    *,
    length: int,
) -> np.ndarray:
    """Return a flattened finite ``Float64`` array of the expected length."""
    try:
        arr = np.asarray(values, dtype=np.float64).reshape(-1)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a numeric sequence") from exc
    if arr.size != length:
        raise ValueError(f"{name} must have length {length} (got {arr.size})")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} must contain only finite values")
    return arr


def validate_tspan(tspan: Sequence[float]) -> tuple[float, float]:
    """Validate and normalize a two-entry integration interval."""
    if len(tspan) != 2:
        raise ValueError("tspan must have exactly two entries")
    t0 = as_finite_float("tspan[0]", tspan[0])
    t1 = as_finite_float("tspan[1]", tspan[1])
    if t1 <= t0:
        raise ValueError(f"tspan end must be greater than start (got {tspan!r})")
    return (t0, t1)
