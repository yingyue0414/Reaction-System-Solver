"""
Author MYING
Using Python 3.9.7
"""
import re
import numpy as np

#import re
#import numpy as np
#Note import re is required for regular expression of strings!

def parse_reaction_string(reaction_string):
    """Parse a reaction string into reactants, products, and rate constant.

    Args:
        reaction_string (str): The input reaction string.

    Returns:
        tuple: A tuple containing three elements - reactants (str), products (str), and rate constant (str).

    Raises:
        ValueError: If the input reaction string does not have the expected format.
        
    Example usage:
        input_string = "A + 2B -> C, k1"
        reactants, products, rate_constant = parse_reaction_string(input_string)
        print(products)  # Output: C
    """
    # Split the reaction string by '->' to separate reactants and products
    reaction_parts = reaction_string.split('->')
    
    if len(reaction_parts) != 2:
        raise ValueError("Invalid reaction string format: " + reaction_string)
    
    # Extract reactants and the second half of the reaction string
    reactants, second_half_string = reaction_parts[0].strip(), reaction_parts[1].strip()
    
    # Split the second half of the reaction string by ',' to separate products and rate constant
    second_half_parts = second_half_string.split(',')
    
    if len(second_half_parts) != 2:
        raise ValueError("Invalid reaction string format: " + reaction_string)
    
    # Extract products and rate constant
    products, rate_constant = second_half_parts[0].strip(), second_half_parts[1].strip()
    
    return reactants, products, rate_constant

#import re

def parse_stoichiometry_string(reactants_or_products_string):
    """Parse a stoichiometry string into a stoichiometry dictionary.

    Args:
        reactants_or_products_string (str): The input stoichiometry string.

    Returns:
        dict: A dictionary representing the stoichiometry, where keys are elements and values are coefficients.

    Raises:
        ValueError: If the input stoichiometry string has an invalid format.
    """
    # Split the input string into individual reactants by the '+' symbol
    reactants = reactants_or_products_string.split('+')
    
    # Initialize an empty dictionary to store the stoichiometry
    stoichiometry_dict = {}
    
    for reactant in reactants:
        # Remove leading and trailing white spaces
        reactant = reactant.strip()
        
        # Use regular expression to match coefficients and elements in the reactant
        match = re.match(r"([\d.]+)?\s*(\w+)", reactant)
        
        if match:
            coefficient, element = match.groups()
            if coefficient:
                coefficient = float(coefficient)
            else:
                coefficient = 1.0

            # Update the stoichiometry dictionary with the coefficient and element
            if element in stoichiometry_dict:
                stoichiometry_dict[element] += coefficient
            else:
                stoichiometry_dict[element] = coefficient
        else:
            raise ValueError(f"Invalid reactant or product format: {reactant}")

    return stoichiometry_dict


def extract_species_dictionaries_from_reaction_strings(reaction_strings):
    """
    Extract species dictionaries, rate constants, and names from a list of reaction strings.

    Args:
        reaction_strings (list): List of reaction strings.

    Returns:
        tuple: A tuple containing species_names_set (set), rate_constant_names (list),
        reactant_dictionaries (list), and product_dictionaries (list).
    """
    species_names_set = set()
    reactant_dictionaries = []
    product_dictionaries = []
    rate_constant_values = [] # for future use, I think rate constant not from string is better for now?
    rate_constant_names = []

    for reaction_string in reaction_strings:
        # split reaction string into three entries
        reactants, products, rate_constant = parse_reaction_string(reaction_string)

        # convert to dict for reactant and product
        reactant_dict = parse_stoichiometry_string(reactants)
        product_dict = parse_stoichiometry_string(products)

        # add new names to set
        species_names_set.update(reactant_dict.keys())
        species_names_set.update(product_dict.keys())

        # append dictionaries to the respective lists
        reactant_dictionaries.append(reactant_dict)
        product_dictionaries.append(product_dict)
        
        # append rate constant to the list
        rate_constant_names.append(rate_constant)
        
    return species_names_set, rate_constant_names, reactant_dictionaries, product_dictionaries


def parse_reaction_strings(reaction_strings, VERBOSE_MODE=True, dtype=int):
    """
    Parse reaction strings into species names, rate constant names, reactant matrices, and product matrices.
    Avoid using dtype = float as long as there is no decimals in stoichiometry is advised.

    Args:
        reaction_strings (list): List of reaction strings.
        VERBOSE_MODE (bool, optional): If True, print additional information for debugging.
        dtype (type, optional): Data type for matrix values (default: float).

    Returns:
        tuple: A tuple containing species_names (list), rate_constant_names (list),
        reactant_matrix_array (numpy.ndarray), and product_matrix_array (numpy.ndarray).
    """

    rate_constant_values = []  # Placeholder

    # Extract species dictionaries, names, and count
    species_names_set, rate_constant_names, reactant_dictionaries, product_dictionaries = \
        extract_species_dictionaries_from_reaction_strings(reaction_strings)

    species_names = list(species_names_set)
    num_species = len(species_names)
    num_reactions = len(reaction_strings)

    if VERBOSE_MODE:
        print(f"species_names : {species_names}")
        print(f"reactant_dictionaries : {reactant_dictionaries}")
        print(f"product_dictionaries : {product_dictionaries}")

    # Initialize reactant and product matrices
    reactant_matrix_array = np.zeros((num_reactions, num_species), dtype=dtype)
    product_matrix_array = np.zeros((num_reactions, num_species), dtype=dtype)

    # Create a dictionary to map species names to indices
    species_to_index = {species: i for i, species in enumerate(species_names)}

    # Fill in reactant and product matrices
    for i in range(num_reactions):
        for species, stoichiometry in reactant_dictionaries[i].items():
            reactant_matrix_array[i, species_to_index[species]] = stoichiometry
        for species, stoichiometry in product_dictionaries[i].items():
            product_matrix_array[i, species_to_index[species]] = stoichiometry

    return species_names, rate_constant_names, reactant_matrix_array, product_matrix_array