# Copyright (C) 2014-2026. Ben Pruitt & Nick Conway; Wyss Institute
# See LICENSE for full GPLv2 license.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''
tests._leakcheck
~~~~~~~~~~~~~~~~~

Fast, deterministic memory-leak assertions for the test suite.

The primary signal is :mod:`tracemalloc` (stdlib): it tracks Python-object and
``PyMem_Malloc`` allocations -- i.e. the binding layer that primer3-py owns
(e.g. the scratch buffers in ``p3helpers`` and every Python object the bindings
build) -- so it is deterministic, cross-platform and cheap enough to bound
tightly in CI.

What it does NOT see is raw C ``malloc`` inside ``libprimer3`` (the design core,
the ``output_structure`` buffers). A leak living entirely in the C core will not
trip the tracemalloc check; for that use the one-off valgrind/massif profiling
(the repo ships ``valgrind-python.supp``). A best-effort peak-RSS delta is also
available as an optional secondary guard, measured only *after* a warmup so the
allocator's initial arena ramp is excluded (that ramp is bounded and is not a
leak -- see the profiling notes in the PR).
'''
import gc
import sys
import tracemalloc

try:
    import resource

    def _peak_rss_kib():
        # ru_maxrss is bytes on macOS but KiB on Linux.
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return rss / 1024.0 if sys.platform == 'darwin' else float(rss)
except ImportError:  # pragma: no cover - Windows
    def _peak_rss_kib():
        return 0.0


def measure(fn, iters, warmup):
    '''Run ``fn()`` ``warmup`` then ``iters`` times.

    Returns ``(tracemalloc_growth_kib, peak_rss_growth_kib)`` measured across
    the ``iters`` calls only (the warmup allocations are excluded).
    '''
    for _ in range(warmup):
        fn()
    gc.collect()

    started_here = not tracemalloc.is_tracing()
    if started_here:
        tracemalloc.start()
    start_traced = tracemalloc.get_traced_memory()[0]
    start_rss = _peak_rss_kib()

    for _ in range(iters):
        fn()

    gc.collect()
    end_traced = tracemalloc.get_traced_memory()[0]
    end_rss = _peak_rss_kib()
    if started_here:
        tracemalloc.stop()

    return (end_traced - start_traced) / 1024.0, end_rss - start_rss


def assert_no_leak(
        testcase,
        fn,
        iters,
        warmup=None,
        max_tracemalloc_kib=128,
        max_peak_rss_kib=None,
):
    '''Assert that repeatedly calling ``fn`` does not leak memory.

    Args:
        testcase: The :class:`unittest.TestCase` (for its assert helpers)
        fn: Zero-arg callable exercising the code under test
        iters: Measured iterations after warmup
        warmup: Warmup iterations (default ``max(10, iters // 5)``) run before
            measurement so caches/allocator arenas reach steady state
        max_tracemalloc_kib: Hard bound on retained Python/PyMem growth
        max_peak_rss_kib: Optional bound on peak-RSS growth (raw-malloc paths);
            leave ``None`` for paths whose RSS ramp is slow/noisy (e.g. design)
    '''
    if warmup is None:
        warmup = max(10, iters // 5)
    tm_kib, rss_kib = measure(fn, iters, warmup)
    testcase.assertLess(
        tm_kib,
        max_tracemalloc_kib,
        f'Python/PyMem allocation grew {tm_kib:.1f} KiB over {iters} '
        f'iterations (limit {max_tracemalloc_kib} KiB, after {warmup} warmup '
        f'iterations) -- possible memory leak',
    )
    if max_peak_rss_kib is not None:
        testcase.assertLess(
            rss_kib,
            max_peak_rss_kib,
            f'peak RSS grew {rss_kib:.1f} KiB over {iters} iterations '
            f'(limit {max_peak_rss_kib} KiB, after {warmup} warmup '
            f'iterations) -- possible memory leak',
        )
    return tm_kib, rss_kib
