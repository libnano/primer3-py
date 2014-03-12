import random
import unittest

from primer3 import bindings, wrappers

class TestLowLevelBindings(unittest.TestCase):

    def randArgs(self):
        self.seq1 = ''.join([random.choice('ATGC') for _ in 
                             range(random.randint(20, 60))])
        self.seq2 = ''.join([random.choice('ATGC') for _ in 
                             range(random.randint(20, 60))])
        self.mv_conc = random.uniform(1, 200)
        self.dv_conc = random.uniform(0, 40)
        self.dntp_conc = random.uniform(0, 20)
        self.dna_conc = random.uniform(0, 200)
        self.temp_c = random.randint(10, 70)
        self.max_loop = random.randint(10, 30)

    def test_calcTm(self):
        for x in range(25):
            self.randArgs()
            # The oligotm executable requires mv_conc and dna_conc to be ints
            # (technically longs... see oligotm_main.c vs. oligotm.c 
            # discrepency)
            binding_tm = bindings.calcTm(
                                seq=self.seq1,
                                mv_conc=int(self.mv_conc),
                                dv_conc=self.dv_conc,
                                dntp_conc=self.dntp_conc,
                                dna_conc=int(self.dna_conc))
            wrapper_tm = wrappers.calcTm(
                                seq=self.seq1,
                                mv_conc=int(self.mv_conc),
                                dv_conc=self.dv_conc,
                                dntp_conc=self.dntp_conc,
                                dna_conc=int(self.dna_conc))
            self.assertEqual(int(binding_tm), int(wrapper_tm))

    def test_calcHairpin(self):
        for x in range(25):
            self.randArgs()
            binding_res = bindings.calcHairpin(
                                seq=self.seq1,
                                mv_conc=self.mv_conc,
                                dv_conc=self.dv_conc,
                                dntp_conc=self.dntp_conc,
                                dna_conc=self.dna_conc,
                                temp_c=self.temp_c,
                                max_loop=self.max_loop)
            wrapper_res = wrappers.calcHairpin(
                                seq=self.seq1,
                                mv_conc=self.mv_conc,
                                dv_conc=self.dv_conc,
                                dntp_conc=self.dntp_conc,
                                dna_conc=self.dna_conc,
                                temp_c=self.temp_c,
                                max_loop=self.max_loop)
            if not wrapper_res:
                self.assertTrue(binding_res.no_structure == 1)
            else:
                self.assertEqual(int(binding_res.tm), int(wrapper_res.tm))

    def test_calcHomodimer(self):
        for x in range(25):
            self.randArgs()
            binding_res = bindings.calcHomodimer(
                                seq=self.seq1,
                                mv_conc=self.mv_conc,
                                dv_conc=self.dv_conc,
                                dntp_conc=self.dntp_conc,
                                dna_conc=self.dna_conc,
                                temp_c=self.temp_c,
                                max_loop=self.max_loop)
            wrapper_res = wrappers.calcHomodimer(
                                seq=self.seq1,
                                mv_conc=self.mv_conc,
                                dv_conc=self.dv_conc,
                                dntp_conc=self.dntp_conc,
                                dna_conc=self.dna_conc,
                                temp_c=self.temp_c,
                                max_loop=self.max_loop)
            if not wrapper_res:
                self.assertTrue(binding_res.no_structure == 1)
            else:
                self.assertEqual(int(binding_res.tm), int(wrapper_res.tm))


    def test_calcHeterodimer(self):
        for x in range(25):
            self.randArgs()
            binding_res = bindings.calcHeterodimer(
                                seq1=self.seq1,
                                seq2=self.seq2,
                                mv_conc=self.mv_conc,
                                dv_conc=self.dv_conc,
                                dntp_conc=self.dntp_conc,
                                dna_conc=self.dna_conc,
                                temp_c=self.temp_c,
                                max_loop=self.max_loop)
            wrapper_res = wrappers.calcHeterodimer(
                                seq1=self.seq1,
                                seq2=self.seq2,
                                mv_conc=self.mv_conc,
                                dv_conc=self.dv_conc,
                                dntp_conc=self.dntp_conc,
                                dna_conc=self.dna_conc,
                                temp_c=self.temp_c,
                                max_loop=self.max_loop)
            if not wrapper_res:
                self.assertTrue(binding_res.no_structure == 1)
            else:
                self.assertEqual(int(binding_res.tm), int(wrapper_res.tm))


if __name__ == '__main__':
    unittest.main()
