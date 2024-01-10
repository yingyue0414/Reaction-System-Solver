# Reaction System Solver
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
  - [3.5. Gillespie](#gillespie)
- [4. Testing] (#testing)

---

## 1. Introduction <a name="introduction"></a>

The Reaction Kinetics solver is a collection of functions that simplify the handling and analysis of chemical reaction systems. It includes functionalities for parsing reaction strings, generating reactant and product matrices, solving reaction kinetics, and modifying reaction rates.

## 2. Installation ande dependency <a name="installation"></a>

No installation is required for this solver. Simply include the provided Python script in your project and import the necessary functions as needed.

### 2.1. Dependency

- Numpy
- Matplotlib
- Scipy (default solver in `odesolver` module)

## 3. Usage <a name="usage"></a>



Here is a concise summary of the provided Python script:


### ReactionStringParser Class Documentation

This script contains a Python class, `ReactionStringParser`, designed for parsing and manipulating chemical reaction strings.

#### Usage
```python
from <dir>.reaction_string_parser import *
```

#### Class Attributes
- `rightward_reaction_symbol`: Rightward reaction symbol (default: "-+>")
- `leftward_reaction_symbol`: Leftward reaction symbol (default: "<-+")
- `reversible_reaction_symbol`: Reversible reaction symbol (default: "<-+>")
- `reaction_rate_separator`: Reaction rate separator (default: "[,;]")
- `species_separator`: Species separator (default: '\+')
- `reaction_rate_value_assigner`: Reaction rate value assigner (default: "=")
- `stoich_species_regex`: Regex pattern for stoichiometry (default: "([\d.]+|\d+\s*\/\s*\d+)?\s*(\w+)")
- `is_rate_constant_required`: Boolean indicating if rate constants are required (default: True)

#### Class Methods
- `parse_reaction_string(reaction_string)`: Parse a reaction string into components.
- `parse_stoichiometry_string(reactants_or_products_string)`: Parse stoichiometry strings.
- `extract_species_dictionaries_from_reaction_strings(reaction_strings)`: Extract dictionaries and rate constants.
- `parse_reaction_strings(reaction_strings, dtype=int, sort_reactions_by=None, sort_species_by=None, VERBOSE_MODE=False)`: Parse and sort reaction strings.
- `sort_by_rate_constants(reactant_matrix, product_matrix, rate_constant_names, sort_order)`: Sort matrices based on rate constants.
- `sort_by_species_names(reactant_matrix, product_matrix, species_names, sort_order)`: Sort matrices based on species names.

#### Example Usage
```python
reaction_strings = ["A + B -> C, kon", "2X -> Y, kf", "C -> A + B, koff", "Y + A -> X + C, ki"]
parser = ReactionStringParser()
species_names, rate_constant_names, reactant_matrix, product_matrix = parser.parse_reaction_strings(reaction_strings)
```

Note: Detailed method descriptions and examples are provided in the docstrings within the script.




### 3.1. Parsing Reaction Strings <a name="parsing-reaction-strings"></a>

This Python script contains a class for parsing and manipulating chemical reaction strings. The `ReactionStringParser` class provides a range of methods to parse reaction strings, extract species dictionaries and rate constants, and sort reactions and species. To use this script, follow the example usages provided below for various methods of the `ReactionStringParser` class.

#### Import the `ReactionStringParser` class

```python
from <dir>.reaction_string_parser import ReactionStringParser
```

#### Initialize the `ReactionStringParser` class

```python
parser = ReactionStringParser()
```

### Class Details

The `ReactionStringParser` class provides the following attributes and methods:

#### Attributes

- `__rightward_reaction_symbol`: Rightward reaction symbol (default: `-+>`)
- `__leftward_reaction_symbol`: Leftward reaction symbol (default: `<-+`)
- `__reversible_reaction_symbol`: Reversible reaction symbol (default: `<-+>`)
- `__reaction_rate_separator`: Reaction rate separator (default: `,;`)
- `__species_separator`: Species separator (default: `+`)
- `__reaction_rate_value_assigner`: Reaction rate value assigner (default: `=`)
- `__stoich_species_regex`: Regex pattern for stoichiometry (default: `([\d.]+|\d+\s*\/\s*\d+)?\s*(\w+)`)
- `DEBUG_MODE`: Debug mode (default: False)

#### Functions

- `parse_reaction_string(reaction_string)`: Parse a reaction string into components.
- `parse_stoichiometry_string(reactants_or_products_string)`: Parse stoichiometry strings.
- `extract_species_dictionaries_from_reaction_strings(reaction_strings)`: Extract dictionaries and rate constants.
- `parse_reaction_strings(reaction_strings, dtype=int, sort_reactions_by=None, sort_species_by=None, VERBOSE_MODE=False)`: Parse and sort reaction strings.
- `sort_by_rate_constants(reactant_matrix, product_matrix, rate_constant_names, sort_order)`: Sort matrices based on rate constants.
- `sort_by_species_names(reactant_matrix, product_matrix, species_names, sort_order)`: Sort matrices based on species names.

Please refer to the docstrings for each method for detailed explanations and examples.

### 3.2. Reaction Matrices <a name="reaction-matrices"></a>

Define the reactant and product matrices for a set of reactions.

```python
# Example:
# A + B -> C
# C -> D + E

# species :                  A, B, C, D, E
reactant_matrix = np.array([[1, 1, 0, 0, 0],  # A + B
                            [0, 0, 1, 0, 0]]) # C
product_matrix  = np.array([[0, 0, 1, 0, 0],  # C
                            [0, 0, 0, 1, 1]]) # D + E
                        
# note that you can calculate the matrix for total yield / loss of species for reactions
# as follows:
delta_matrix = product_matrix - reactant_matrix
```

### 3.3. Solving Reaction Kinetics <a name="solving-reaction-kinetics"></a>

In the module, I utilize the `solve_ivp` function from SciPy to solve the kinetics of a chemical reaction system. Plot the results using Matplotlib.
However this is just to showcase the feasbility of solving ODE through constructing matrices. i.e. after constructing the dydt function,
you can use any method of solving the system of ODEs.

```python
# Solve the kinetics using scipy.integrate.solve_ivp
from src.stringparser import *
from src.odesolver import *

rsp = ReactionStringParser()

# Example usage:
# Generate reaction system
reaction_strings = ["E + S -> ES", "ES -> E + P"]
species_names, rate_constant_names, reactant_matrix, product_matrix =\
        rsp.parse_reaction_strings(reaction_strings, sort_species_by="increasing")

print("Species Names:")
print(species_names)
print("Rate Constant Names:")
print(rate_constant_names)

# Rate constant assuming already non-dimensionalized
rate_constants = [0.01, 1.0]

# Define time span and initial concentration, assuming already non-dimensionalized
t_span = [0.0, 200.0]
y_init = [5.0, 0.0, 0.0, 100.0]

solve_reaction_ode(dydt, t_span, y_init, reactant_matrix = reactant_matrix,
                                             product_matrix = product_matrix,
                                             k = rate_constants,
                                             species_names = species_names)
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


Certainly! Here's the added section about the Gillespie algorithm:



### 3.5. Gillespie Algorithm Simulation <a name="gillespie"></a>

The Gillespie algorithm is a stochastic simulation algorithm used to model the time evolution of chemical reaction systems. It is particularly useful for systems with low reactant counts and when stochastic effects play a significant role. The `gillespie_simulation` function provided in `src.reaction_gillespie` module (src/reaction_gillespie.py) allows you to perform Gillespie simulations for a given chemical reaction system. `convert_to_microscopic_rate_constants` enables automatic conversion between microscopic and macroscopic rate constant (assuming concentration in M and volume in L).

#### Function Signature

```python
def gillespie_simulation(max_time, y_init, reactant_matrix, product_matrix, microscopic_rate_constants, full_update_scheme=True):
    """
    Perform Gillespie simulation for a chemical reaction system.

    Args:
        max_time (float): Maximum simulation time.
        y_init (numpy.ndarray): Initial state of the system (species counts).
        reactant_matrix (numpy.ndarray): Matrix representing reactants in each reaction.
        product_matrix (numpy.ndarray): Matrix representing products in each reaction.
        microscopic_rate_constants (numpy.ndarray): Rate constants for each reaction.
        full_update_scheme (bool): Controls if update every propensity entry in each iteration. If
                set to False, optimize the simulation so that not every propensity will be updated each step.

    Returns:
        tuple: A tuple containing arrays for recorded time points (t_record) and
               corresponding system states (y_record).
    Note:
        This function performs a Gillespie simulation for a chemical reaction system.
        It records the system state and corresponding time points during the simulation.
    """
```



#### Example Usage

```python
max_time = 100.0
y_init = np.array([10, 5, 3])  # Initial state (species counts)

macroscopic_rate_constants = np.array([1.0e6, 5.0])

reactant_matrix = np.array([[1, 1, 0], [0, 0, 2]])

product_matrix = np.array([[0, 0, 1], [1, 2, 1]])

y_record, t_record = gillespie_simulation(max_time, y_init,
                                           reactant_matrix, product_matrix,
                                           microscopic_rate_constants)
print(y_record)
print(t_record)
```

Feel free to use the `gillespie_simulation` function to explore the stochastic behavior of your chemical reaction system. Adjust the input parameters as needed for your specific simulation requirements.

## 4. Testing <a name="introduction"></a>

Unit test is written with `pytest6.2.4`. Run all tests in main directory with command `pytest`. 