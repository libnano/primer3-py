'''
primer3.tests
~~~~~~~~~~~~~

Unit tests for the primer3-py package.

'''

import unittest

import test_lowlevel
import test_primerdesign


def runTests():
    tl = unittest.TestLoader()
    lowLevelSuite = tl.loadTestsFromTestCase(
                        test_lowlevel.TestLowLevelBindings)
    unittest.TextTestRunner(verbosity=2).run(lowLevelSuite)
    designSuite = tl.loadTestsFromTestCase(
                    test_primerdesign.TestDesignBindings)
    unittest.TextTestRunner(verbosity=2).run(designSuite)