import unittest

import test_lowlevel
import test_primerdesign

tl = unittest.TestLoader()
lowLevelSuite = tl.loadTestsFromTestCase(
                    test_lowlevel.TestLowLevelBindings)
unittest.TextTestRunner(verbosity=2).run(lowLevelSuite)
designSuite = tl.loadTestsFromTestCase(
                test_primerdesign.TestDesignBindings)
unittest.TextTestRunner(verbosity=2).run(designSuite)