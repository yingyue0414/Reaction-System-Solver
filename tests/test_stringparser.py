import unittest
import numpy as np
from numpy.testing import assert_array_equal
from src.stringparser import ReactionStringParser  # Replace 'your_module' with the actual module name

class TestReactionStringParser(unittest.TestCase):

    def setUp(self):
        # Initialize the ReactionStringParser with default values
        self.parser = ReactionStringParser()

    def test_parse_reaction_string(self):

        input_string = "A + 2B -> C, k1"
        expected_result = ('A + 2B', 'C', 'k1', 1)
        result = self.parser.parse_reaction_string(input_string)
        self.assertEqual(result, expected_result)

        input_string_2 = "2X + Y <- , k2"
        expected_result_2 = ('2X + Y', '', 'k2', -1)
        result_2 = self.parser.parse_reaction_string(input_string_2)
        self.assertEqual(result_2, expected_result_2)

        input_string_3 = "P + Q <-> R, k3, k4"
        expected_result_3 = ('P + Q', 'R', ('k3', 'k4'), 0)
        result_3 = self.parser.parse_reaction_string(input_string_3)
        self.assertEqual(result_3, expected_result_3)

    def test_parse_stoichiometry_string(self):

        input_string = "A + 1.42857  B   +1/2C+  A"
        expected_result = {'A': 2.0, 'B': 1.42857, 'C': 0.5}
        result = self.parser.parse_stoichiometry_string(input_string)
        self.assertEqual(result, expected_result)

        input_string_2 = "3D + E + F"
        expected_result_2 = {'D': 3.0, 'E': 1.0, 'F': 1.0}
        result_2 = self.parser.parse_stoichiometry_string(input_string_2)
        self.assertEqual(result_2, expected_result_2)

        input_string_3 = "9XX + 2XXX + X"
        expected_result_3 = {'X': 1.0, 'XX': 9.0, 'XXX': 2.0}
        result_3 = self.parser.parse_stoichiometry_string(input_string_3)
        self.assertEqual(result_3, expected_result_3)

    def test_extract_species_dictionaries_from_reaction_strings(self):
        reaction_strings = ["A + B -> C, kon", "2X -> Y, kf"]
        expected_result = (
            {'A', 'B', 'C', 'X', 'Y'},
            ['kon', 'kf'],
            [{'A': 1, 'B': 1}, {'X': 2}],
            [{'C': 1}, {'Y': 1}]
        )
        result = self.parser.extract_species_dictionaries_from_reaction_strings(reaction_strings)
        self.assertCountEqual(result[0], expected_result[0])  # Compare sets ignore order
        self.assertCountEqual(result[1], expected_result[1])  # Compare lists ignore order
        self.assertCountEqual(result[2], expected_result[2])  # Compare lists of dictionaries
        self.assertCountEqual(result[3], expected_result[3])  # Compare lists of dictionaries

    
    def test_parse_reaction_strings(self):
        reaction_strings = ["A + B -> C, kon", "2X -> Y, kf"]
        expected_result = (['A', 'B', 'C', 'X', 'Y'],
                                ['kon', 'kf'],
                                np.array([[1, 1, 0, 0, 0],
                                [0, 0, 0, 2, 0]]),
                                np.array([[0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 1]]))
        result = self.parser.parse_reaction_strings(reaction_strings)
        assert_array_equal(result[0], expected_result[0]) # ignore order
        assert_array_equal(result[1], expected_result[1]) # ignore order
        assert_array_equal(result[2], expected_result[2]) # compare matrices
        assert_array_equal(result[3], expected_result[3]) # compare matrices

        # complicated case
        reaction_strings_2 = [
            "A + DA -> DA_d",
            "DA_d -> A + DA",
            "DA -> DA + MA",
            "DA_d -> DA_d + MA",
            "A + DR -> DR_d",
            "DR_d -> A + DR",
            "DR -> DR + MR",
            "DR_d -> DR_d + MR",
            "MA ->",
            "MA -> MA + A",
            "MR ->",
            "MR -> MR + R",
            "A ->",
            "R ->",
            "A + R -> C",
            "C -> R"
            ]

        result_2 = self.parser.parse_reaction_strings(reaction_strings_2)
        expected_result_2 = (np.array(['A', 'C', 'DA', 'DA_d', 'DR', 'DR_d', 'MA', 'MR', 'R']),
                             ['k0', 'k1', 'k2', 'k3', 'k4', 'k5', 'k6', 'k7', 'k8', 'k9', 'k10', 'k11', 'k12', 'k13', 'k14', 'k15'],
                             np.array([[1, 0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 0, 0, 0, 0, 0],
                                [0, 0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 1, 0, 0, 0, 0, 0],
                                [1, 0, 0, 0, 1, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 1, 0, 0, 0],
                                [0, 0, 0, 0, 1, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 1, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 1, 0],
                                [0, 0, 0, 0, 0, 0, 0, 1, 0],
                                [1, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 1],
                                [1, 0, 0, 0, 0, 0, 0, 0, 1],
                                [0, 1, 0, 0, 0, 0, 0, 0, 0]]), 
                                np.array([[0, 0, 0, 1, 0, 0, 0, 0, 0],
                                [1, 0, 1, 0, 0, 0, 0, 0, 0],
                                [0, 0, 1, 0, 0, 0, 1, 0, 0],
                                [0, 0, 0, 1, 0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 0, 1, 0, 0, 0],
                                [1, 0, 0, 0, 1, 0, 0, 0, 0],
                                [0, 0, 0, 0, 1, 0, 0, 1, 0],
                                [0, 0, 0, 0, 0, 1, 0, 1, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [1, 0, 0, 0, 0, 0, 1, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 1, 1],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 1, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0, 1]]))
        assert_array_equal(result_2[0], expected_result_2[0]) # ignore order
        assert_array_equal(result_2[1], expected_result_2[1]) # ignore order
        assert_array_equal(result_2[2], expected_result_2[2]) # compare matrices
        assert_array_equal(result_2[3], expected_result_2[3]) # compare matrices

    # Add more test cases for other methods as needed

if __name__ == '__main__':
    unittest.main()
