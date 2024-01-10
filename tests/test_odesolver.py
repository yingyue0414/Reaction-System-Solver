import unittest
import numpy as np
from src.odesolver import calculate_macroscopic_reaction_rates, dydt, solve_reaction_ode

class TestODESolver(unittest.TestCase):

    def test_calculate_macroscopic_reaction_rates(self):
        y = [1.0, 0.5]
        reactant_matrix = np.array([[1, 0], [0, 1]])
        k = [0.1, 0.2]

        result = calculate_macroscopic_reaction_rates(y, reactant_matrix, k)
        expected_result = np.array([0.1, 0.1])
        np.testing.assert_allclose(result, expected_result)

    def test_dydt(self):
        # Model system:
        # E + S -> ES, k1
        # ES -> E + P, k2
        #
        # y = [[E], [S], [ES], [P]]
        t = 0.1
        y = [0.1, 50.0, 0.1, 1.5]
        reactant_matrix = np.array([[1, 1, 0, 0], [0, 0, 1, 0]])
        product_matrix = np.array([[0, 0, 1, 0], [1, 0, 0, 1]])
        k = [100.0, 1.0]

        result = dydt(t, y, reactant_matrix, product_matrix, k)
        print(result)
        expected_result = np.array([-499.9, -500.0, 499.9, 0.1])
        np.testing.assert_allclose(result, expected_result)

    def test_solve_reaction_ode(self):
        t_span = (0, 10)
        y_initial = [1.0, 0.5]
        reactant_matrix = np.array([[1, 0], [0, 1]])
        product_matrix = np.array([[0, 1], [1, 0]])
        k = [0.1, 0.2]

        with self.subTest(msg="Check if solve_reaction_ode runs without errors"):
            solve_reaction_ode(dydt, t_span, y_initial, reactant_matrix, product_matrix, k, plotting=False)

if __name__ == '__main__':
    unittest.main()
