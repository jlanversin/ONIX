"""Select and use the numerical backend for depletion solves."""

from contextlib import nullcontext

import numpy as np
from scipy import linalg as scipy_linalg


_BACKEND_NAME = 'cpu'
_DEVICE_ID = 0
_DEVICE_NAME = None
_CUPY = None
_CUPYX_LINALG = None


def _load_cupy_modules():
    global _CUPY, _CUPYX_LINALG

    if _CUPY is None or _CUPYX_LINALG is None:
        import cupy
        import cupyx.scipy.linalg

        _CUPY = cupy
        _CUPYX_LINALG = cupyx.scipy.linalg

    return _CUPY, _CUPYX_LINALG


def _normalize_device_name(raw_name):
    if isinstance(raw_name, bytes):
        return raw_name.decode()
    return str(raw_name)


def set_compute_backend(backend='cpu', device_id=0):
    """Set the numerical backend used by ONIX.

    Parameters
    ----------
    backend: str
        Either ``'cpu'`` or ``'gpu'``.
    device_id: int
        CUDA device index to use when ``backend`` is ``'gpu'``.
    """

    backend_name = '{}'.format(backend).strip().lower()
    if backend_name not in ('cpu', 'gpu'):
        raise ValueError("backend must be either 'cpu' or 'gpu'")

    global _BACKEND_NAME, _DEVICE_ID, _DEVICE_NAME

    if backend_name == 'cpu':
        _BACKEND_NAME = 'cpu'
        _DEVICE_ID = 0
        _DEVICE_NAME = None
        return describe_compute_backend()

    cupy, _ = _load_cupy_modules()

    try:
        device_count = cupy.cuda.runtime.getDeviceCount()
    except Exception as exc:
        raise RuntimeError('CuPy is installed but no CUDA runtime is available') from exc

    if device_id < 0 or device_id >= device_count:
        raise ValueError(
            'Requested CUDA device {} but only {} device(s) are available'.format(
                device_id,
                device_count,
            )
        )

    try:
        with cupy.cuda.Device(device_id):
            # Force a small allocation so backend selection fails early on bad runtimes.
            cupy.zeros(1, dtype=cupy.float64)
            cupy.cuda.Stream.null.synchronize()
            device_props = cupy.cuda.runtime.getDeviceProperties(device_id)
    except Exception as exc:
        raise RuntimeError(
            'Failed to activate the GPU backend on CUDA device {}'.format(device_id)
        ) from exc

    _BACKEND_NAME = 'gpu'
    _DEVICE_ID = device_id
    _DEVICE_NAME = _normalize_device_name(device_props['name'])

    return describe_compute_backend()


def get_compute_backend():
    """Return the active backend name."""

    return _BACKEND_NAME


def describe_compute_backend():
    """Return a human-readable description of the active backend."""

    if _BACKEND_NAME == 'gpu':
        if _DEVICE_NAME is None:
            return 'gpu:{}'.format(_DEVICE_ID)
        return 'gpu:{} ({})'.format(_DEVICE_ID, _DEVICE_NAME)

    return 'cpu'


def gpu_available(device_id=0):
    """Return whether CuPy can see and use the requested CUDA device."""

    try:
        cupy, _ = _load_cupy_modules()
        device_count = cupy.cuda.runtime.getDeviceCount()
        if device_id < 0 or device_id >= device_count:
            return False
        with cupy.cuda.Device(device_id):
            cupy.zeros(1, dtype=cupy.float64)
            cupy.cuda.Stream.null.synchronize()
        return True
    except Exception:
        return False


def backend_device():
    """Return a context manager for the active backend device."""

    if _BACKEND_NAME != 'gpu':
        return nullcontext()

    cupy, _ = _load_cupy_modules()
    return cupy.cuda.Device(_DEVICE_ID)


def get_array_module():
    """Return the active array module."""

    if _BACKEND_NAME == 'gpu':
        cupy, _ = _load_cupy_modules()
        return cupy

    return np


def get_expm_module():
    """Return the module that provides ``expm`` for the active backend."""

    if _BACKEND_NAME == 'gpu':
        _, cupyx_linalg = _load_cupy_modules()
        return cupyx_linalg

    return scipy_linalg


def asarray(values, dtype=None):
    """Convert values to the active backend array type."""

    array_module = get_array_module()
    return array_module.asarray(values, dtype=dtype)


def to_numpy(values):
    """Move backend arrays back to NumPy arrays."""

    if _BACKEND_NAME == 'gpu':
        cupy, _ = _load_cupy_modules()
        if isinstance(values, cupy.ndarray):
            return cupy.asnumpy(values)

    return np.asarray(values)


def synchronize():
    """Block until queued GPU work has completed."""

    if _BACKEND_NAME == 'gpu':
        cupy, _ = _load_cupy_modules()
        cupy.cuda.Stream.null.synchronize()


def solve(matrix, rhs):
    """Solve a dense linear system on the active backend."""

    array_module = get_array_module()

    with backend_device():
        matrix_backend = array_module.asarray(matrix)
        rhs_backend = array_module.asarray(rhs)
        solution = array_module.linalg.solve(matrix_backend, rhs_backend)
        synchronize()

    return to_numpy(solution)
