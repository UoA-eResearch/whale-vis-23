from contextlib import contextmanager
from time import perf_counter


@contextmanager
def timer(name):
    """Context manager to time a block of code"""
    t0 = perf_counter()
    yield
    t1 = perf_counter()
    print(f'[{name}] done in {t1 - t0:.3f} s')
