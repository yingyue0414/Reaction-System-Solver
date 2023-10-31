# Reaction IVP Solver
This solver provides a set of functions for working with chemical reaction systems, including parsing reaction strings, solving reaction kinetics, and modifying reaction rates. It is written in Python 3.9.7 and utilizes libraries such as NumPy, Matplotlib, and SciPy.

Note that the module main serves as a tool to generate the `dydt` function from strings of reactions. The module / method to solve IVP for these `dydt` can be customly defined after such function is generated.

## WARNING
The following content of this `README.md` suggests the end goal of this module, not the current progress. Currently, the module is not wrapped to a
package, but is already able to intepret reaction strings and solve IVP using scipy.

## Table of Contents
- [1. Introduction](#introduction)
- [2. Installation](#installation)
- [3. Usage](#usage)
  - [3.1. Parsing Reaction Strings](#parsing-reaction-strings)
  - [3.2. Reaction Matrices](#reaction-matrices)
  - [3.3. Solving Reaction Kinetics](#solving-reaction-kinetics)
  - [3.4. Modifying the Model wtih Decorators](#decorators)

---

## 1. Introduction <a name="introduction"></a>

The Reaction Kinetics solver is a collection of functions that simplify the handling and analysis of chemical reaction systems. It includes functionalities for parsing reaction strings, generating reactant and product matrices, solving reaction kinetics, and modifying reaction rates.

## 2. Installation <a name="installation"></a>

No installation is required for this solver. Simply include the provided Python script in your project and import the necessary functions as needed.

```python
from reaction_kinetics_solver import parse_reaction_string, parse_stoichiometry_string, \
    extract_species_dictionaries_from_reaction_strings, parse_reaction_strings, \
    calculate_macroscopic_reaction_rates, dydt, dydt_scalar_decorator
```

## 3. Usage <a name="usage"></a>

### 3.1. Parsing Reaction Strings <a name="parsing-reaction-strings"></a>

Use the provided functions to parse reaction strings into their constituent parts - reactants, products, and rate constants.

```python
# Example usage:
input_string = "A + 2B -> C, k1"
reactants, products, rate_constant, direction = parse_reaction_string(input_string)
print("Reactants:", reactants)  # Output: A + 2B
print("Products:", products)    # Output: C
print("Rate Constant:", rate_constant)  # Output: k1
print("Direction:", direction)   # Output: 1
```

### 3.2. Reaction Matrices <a name="reaction-matrices"></a>

Define the reactant and product matrices for a set of reactions and solve the kinetics using the provided functions.

```python
# Sample usage:
t_span = [0, 20]
y_initial = [3, 10, 7]
k = [1, 1, 1, 1]
species_names, reactant_matrix, product_matrix, rate_constant_names = \
    parse_reaction_strings(reaction_strings)

# Solve the kinetics using scipy.integrate.solve_ivp
sol = solve_ivp(dydt, t_span, y_initial, args=(reactant_matrix, product_matrix, k), dense_output=True)
```

### 3.3. Solving Reaction Kinetics <a name="solving-reaction-kinetics"></a>

In the module, I utilize the `solve_ivp` function from SciPy to solve the kinetics of a chemical reaction system. Plot the results using Matplotlib. However, after constructing the dydt function,
you can use any method of solving the system of ODEs.

```python
t = np.linspace(0, 20, 300)
y = sol.sol(t)
plt.plot(t, y.T)
plt.xlabel('Time')
plt.ylabel('Concentration')
plt.legend(np.arange(len(y_initial))
plt.show()
```

### 3.4. Modifying the Model wtih Decorators <a name="decorators"></a>

An example usage of a decorator is as follows: assuming we would like to apply a scalar factor (ex. set some concentration to be const) to reaction rates using the `dydt_scalar_decorator` function. This allows for the modulation of reaction rates, which can be useful for sensitivity analyses or parameter variations.

```python
scalar = np.array([1.5, 1, 4, 1])  # Scalar values for each rate constant
dydtmodified = dydt_scalar_decorator(dydt)

# Solve the kinetics with modified rates
solmodified = solve_ivp(dydtmodified, t_span, y_initial,
                        args=(reactant_matrix, product_matrix, k, scalar), dense_output=True)
```

---

Feel free to explore and adapt the Reaction Kinetics solver to meet your specific needs for chemical reaction analysis. If you have any questions or encounter issues, please refer to the code comments and documentation for additional information.
