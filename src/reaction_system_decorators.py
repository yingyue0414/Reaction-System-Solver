"""
Author MYING
Using Python 3.9.7
"""
import numpy as np

def dydt_scalar_decorator(dydt):
    """
    Decorator for a function that calculates the rate of change of a variable in an ODE system, modifying one of the arguments
    by scaling it with a scalar.

    This decorator is designed for functions that calculate the rate of change in an ODE system and allows you to apply a
    scalar transformation to one of the arguments, specified by its index. The decorator takes the last argument of the
    decorated function as the scalar and the index of the argument to be multiplied by the scalar.

    Parameters:
    - dydt (callable): The original function that calculates the rate of change.
    
    Returns:
    - callable: A decorated function that modifies one of the function's arguments by scaling it with a scalar and then
      calls the original function.

    Raises:
    - TypeError: If the decorated function is called with fewer than two positional arguments, indicating that it cannot be
      decorated as it needs at least two arguments to work properly.

    Example:
    ```python
    @dydt_scalar_decorator
    def my_dydt(t, y, a, b, c):
        # Your dydt calculation logic here
        return result

    # Decorated function usage
    decorated_result = my_dydt(0.0, y0, a=2, b=3, c=4, scalar=0.5)
    ```
    In the example, the `scalar` argument (0.5) is used to scale the argument at index 4 (c) before passing it to the `my_dydt`
    function for rate calculation.
    """
    def wrapper(*args):
        if len(args) < 2:
            raise TypeError("The decorated function requires at least two positional arguments.")
            
        # Extract the scalar argument from the last position
        scalar = args[-1]

        # Extract the argument to be multiplied by scalar
        arg_index_to_multiply = 4
        arg_to_multiply = args[arg_index_to_multiply]

        # Multiply the selected argument by scalar
        args_list = list(args)
        args_list[arg_index_to_multiply] = arg_to_multiply * scalar
        args = tuple(args_list)

        result = dydt(*args[:-1])  # Pass all arguments except the last one (scalar)
        return result
    return wrapper