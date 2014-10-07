
import sys
import unittest

import test_lowlevel
import test_primerdesign

tl = unittest.TestLoader()
lowLevelSuite = tl.loadTestsFromTestCase(
                    test_lowlevel.TestLowLevelBindings)
res1 = unittest.TextTestRunner(verbosity=2).run(lowLevelSuite)
designSuite = tl.loadTestsFromTestCase(
                test_primerdesign.TestDesignBindings)
res2 = unittest.TextTestRunner(verbosity=2).run(designSuite)


success = res1.wasSuccessful() and res2.wasSuccessful()

sys.exit(int(not success))  # Exit 0 on success, 1 on failure