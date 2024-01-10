import unittest
import numpy as np
from src.gillespie import convert_to_microscopic_rate_constants, calculate_propensity, gillespie_simulation

class TestReactionGillespie(unittest.TestCase):

    def setUp(self):
        # Example parameters for testing
        self.macroscopic_rate_constants = np.array([1.0e6, 5.0])
        self.reactant_matrix = np.array([[1, 1, 0], [0, 0, 2]])
        self.product_matrix = np.array([[0, 0, 1], [1, 2, 1]])
        self.volume = 1.0e-18 # Litre!
        self.y = np.array([10, 5, 3])

    def test_convert_to_microscopic_rate_constants(self):
        microscopic_rate_constants = convert_to_microscopic_rate_constants(
            self.macroscopic_rate_constants, self.reactant_matrix, self.volume
        )
        expected_result = np.array([1.66053928e+00, 1.66053928e-05])
        np.testing.assert_allclose(microscopic_rate_constants, expected_result)

    def test_calculate_propensity(self):
        microscopic_rate_constants = np.array([0.1 * 90, 0.05 * 5])
        propensities = calculate_propensity(self.y, self.reactant_matrix, microscopic_rate_constants)
        expected_result = np.array([450, 0.75])
        np.testing.assert_allclose(propensities, expected_result)

    def test_gillespie_simulation(self):
        max_time = 10.0
        y_init = np.array([10, 5, 3])
        
        record_interval = 10
        full_update_scheme = True

        y_record, t_record = gillespie_simulation(
            max_time, y_init, self.reactant_matrix, self.product_matrix,
            self.macroscopic_rate_constants, record_interval, full_update_scheme
        )

        # You can add more assertions based on the expected behavior of your simulation

if __name__ == '__main__':
    unittest.main()
