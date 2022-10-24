"""
A simple decorator to measure the time
of function execution.
"""

from functools import wraps
import time


def timeit(func):
    """
    Returns the time of function execution.
    """
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(
            f'Function {func.__name__}{args} {kwargs} Took {total_time:.1e} seconds')
        return result
    return timeit_wrapper
