#!/usr/bin/env python3.8.10
"""
This script contains functions for simulating simple
gillespie trajectories.

Usage: from <dir>.reaction_gillespie import *

Docstrings and comments are helped written by GPT3.5 and GPT API

@author: MYING
@email: yying7@jh.edu
"""
import numpy as np
import math

def convert_to_microscopic_rate_constants(macroscopic_rate_constants, reactant_matrix, volume,
                                          avogadro=6.02214e23):
    """
    Note:
        All volumes are assumed to be in liters, and concentrations are assumed to be in mol/L!
    
    Convert macroscopic rate constants to microscopic rate constants for Gillespie algorithm.

    Args:
        macroscopic_rate_constants (numpy.ndarray): Array of macroscopic rate constants.
        reactant_matrix (numpy.ndarray): Matrix representing reactants in each reaction.
        volume (float): Volume of the system (assumed to be in liters).
        avogadro (float, optional): Avogadro's number (default: 6.02214e-23).

    Raises:
        ValueError: If any entry in the reactant matrix is not a mathematical integer.

    Returns:
        numpy.ndarray: Array of microscopic rate constants.
    """
    #### ALL volumes are assumed to be in liters, and concentrations are assumed to be in mol/L!
    # Check if all entries in the reactant matrix are mathematical integers
    if not np.all(np.mod(reactant_matrix, 1) == 0):
        raise ValueError("For Gillespie, all entries in the matrix must be mathematical integers.")

    # Initialize an array for microscopic rate constants
    microscopic_rate_constants = np.zeros(len(macroscopic_rate_constants))

    # Calculate microscopic rate constants
    for reaction_index, reaction in enumerate(reactant_matrix):
        scalar = 1
        power = 1
        for species_index, species_count in enumerate(reaction):
            scalar *= math.factorial(int(species_count))
            power -= species_count
        microscopic_rate_constants[reaction_index] = (
            scalar * macroscopic_rate_constants[reaction_index] * np.power((volume * avogadro), power)
        )

    return microscopic_rate_constants

def calculate_propensity(y, reactant_matrix, microscopic_rate_constants):
    """
    Calculate propensities for Gillespie algorithm.

    Args:
        y (numpy.ndarray): Current state of the system (species counts).
        reactant_matrix (numpy.ndarray): Matrix representing reactants in each reaction.
        microscopic_rate_constants (numpy.ndarray): Rate constants for each reaction.

    Returns:
        numpy.ndarray: Array of propensities for each reaction.

    Example:
        Suppose you have the following input matrices and arrays:

        y = np.array([10, 5, 3])  # Current state (species counts)
        
        reactant_matrix = np.array([[2, 1, 0],  # Example reactant matrix
                                    [0, 1, 1]])

        microscopic_rate_constants = np.array([0.1, 0.05])  # Example rate constants

        propensities = calculate_propensity(y, reactant_matrix, microscopic_rate_constants)
        print(propensities)
        # Output: [0.1 * comb(10, 2) * comb(5, 1), 0.05 * comb(5, 1) * comb(3, 1)]

    Note:
        The function calculates the propensity of each reaction in a Gillespie algorithm.
        Propensity is the product of the microscopic rate constant and combinatorial terms
        based on the reactant matrix and current state (y) of the system.
    """
    if not np.all(np.mod(reactant_matrix, 1) == 0):
        raise ValueError("For gillespie, all entries in the matrix must be mathematically integers.")
        
    propensities = np.zeros(len(reactant_matrix))

    # Loop over each reaction
    for reaction_index, reaction in enumerate(reactant_matrix):
        propensity = microscopic_rate_constants[reaction_index]

        # Multiply by the combinatorial term for each reactant
        for species_index, species_count in enumerate(reaction):
            propensity *= math.comb(y[species_index], species_count)

        propensities[reaction_index] = propensity

    return propensities


def gillespie_simulation(max_time, y_init,
                         reactant_matrix, product_matrix, microscopic_rate_constants):
    """
    Perform Gillespie simulation for a chemical reaction system.

    Args:
        max_time (float): Maximum simulation time.
        y_init (numpy.ndarray): Initial state of the system (species counts).
        reactant_matrix (numpy.ndarray): Matrix representing reactants in each reaction.
        product_matrix (numpy.ndarray): Matrix representing products in each reaction.
        microscopic_rate_constants (numpy.ndarray): Rate constants for each reaction.

    Returns:
        tuple: A tuple containing arrays for recorded time points (t_record) and
               corresponding system states (y_record).

    Example:
        Suppose you have the following input matrices and arrays:

        max_time = 100.0
        y_init = np.array([10, 5, 3])  # Initial state (species counts)
        
        reactant_matrix = np.array([[2, 1, 0],  # Example reactant matrix
                                    [0, 1, 1]])

        product_matrix = np.array([[0, 1, 0],  # Example product matrix
                                   [1, 0, 1]])

        microscopic_rate_constants = np.array([0.1, 0.05])  # Example rate constants

        y_record, t_record = gillespie_simulation(max_time, y_init,
                                                   reactant_matrix, product_matrix,
                                                   microscopic_rate_constants)
        print(y_record)
        print(t_record)
        
    Note:
        This function performs a Gillespie simulation for a chemical reaction system.
        It records the system state and corresponding time points during the simulation.
    """
    time = 0.0
    y = y_init  # Initial copy numbers
    propensities = np.zeros(len(reactant_matrix))  # Propensities array
    delta_y = product_matrix - reactant_matrix  # Yield matrix
    index = np.array(range(0, len(reactant_matrix)))  # np.random.choice must be 1-d array; use indexing instead
    y_record = [np.copy(y)]  # Record copy numbers
    t_record = [time]  # Record time

    while time < max_time:  # Control simulation time scale
        # Calculate propensity
        propensities = calculate_propensity(y, reactant_matrix, microscopic_rate_constants)

        # Calculate r_tot and sojourn time
        r_tot = np.sum(propensities)
        tau = - (1.0 / r_tot) * np.log(np.random.rand())

        # Add to species
        y += delta_y[np.random.choice(index, p=propensities / r_tot)]

        # Progress time
        time += tau

        # Record
        y_record.append(np.copy(y))
        t_record.append(time)

    return y_record, t_record
