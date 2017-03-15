# -*- coding: utf-8 -*-

import numpy as np
from .context import sajou as sj
model = sajou.model

import unittest


class BasicTestSuite(unittest.TestCase):
    """Basic test cases."""

    def test_model_creation(self):
        # create model
        test_model = sj.model.Model(name='Test Model', dimensionality="2D")
        self.assertEqual(test_model._name, 'Test Model')


if __name__ == '__main__':
    unittest.main()
