import random
import unittest

from primer3 import bindings, wrappers

class TestLowLevelBindings(unittest.TestCase):

    def randArgs(self):
        self.seq1 = ''.join([random.choice('ATGC') for _ in
                             range(random.randint(20, 59))])
        self.seq2 = ''.join([random.choice('ATGC') for _ in
                             range(random.randint(20, 59))])
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
                try:
                    self.assertEqual(int(binding_res.tm), int(wrapper_res.tm))
                except AssertionError:
                    print(self.seq1, self.mv_conc, self.dv_conc, self.dntp_conc,
                          self.dna_conc, self.temp_c, self.max_loop)
                    raise

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
                try:
                    self.assertEqual(int(binding_res.tm), int(wrapper_res.tm))
                except AssertionError:
                    print(self.seq1, self.mv_conc, self.dv_conc, self.dntp_conc,
                          self.dna_conc, self.temp_c, self.max_loop)
                    raise


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
                try:
                    self.assertEqual(int(binding_res.tm), int(wrapper_res.tm))
                except AssertionError:
                    print(self.seq1, self.mv_conc, self.dv_conc, self.dntp_conc,
                          self.dna_conc, self.temp_c, self.max_loop)
                    raise

class TestDesignBindings(unittest.TestCase):

    def testBindingRepeat(self):
        for x in range(1):
            bindings.designPrimers(
                {
                    'PRIMER_OPT_SIZE': 20,
                    'PRIMER_PICK_INTERNAL_OLIGO': 1,
                    'PRIMER_INTERNAL_MAX_SELF_END': 8,
                    'PRIMER_MIN_SIZE': 18,
                    'PRIMER_MAX_SIZE': 25,
                    'PRIMER_OPT_TM': 60.0,
                    'PRIMER_MIN_TM': 57.0,
                    'PRIMER_MAX_TM': 63.0,
                    'PRIMER_MIN_GC': 20.0,
                    'PRIMER_MAX_GC': 80.0,
                    'PRIMER_MAX_POLY_X': 100,
                    'PRIMER_INTERNAL_MAX_POLY_X': 100,
                    'PRIMER_SALT_MONOVALENT': 50.0,
                    'PRIMER_DNA_CONC': 50.0,
                    'PRIMER_MAX_NS_ACCEPTED': 0,
                    'PRIMER_MAX_SELF_ANY': 12,
                    'PRIMER_MAX_SELF_END': 8,
                    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
                    'PRIMER_PAIR_MAX_COMPL_END': 8,
                },
                {
                    'PRIMER_PRODUCT_SIZE_RANGE': [75,100,100,125,125,150,150,175,175,200,200,225],
                    'SEQUENCE_ID': 'MH1000',
                    'SEQUENCE_TEMPLATE': 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCTCTAGTAG',
                    'SEQUENCE_INCLUDED_REGION': [36,342]
                }
            )

#     def testHuman(self):
#         binding_res = bindings.designPrimers(
#             {
#                 'PRIMER_OPT_SIZE': 20,
#                 'PRIMER_PICK_INTERNAL_OLIGO': 1,
#                 'PRIMER_INTERNAL_MAX_SELF_END': 8,
#                 'PRIMER_MIN_SIZE': 18,
#                 'PRIMER_MAX_SIZE': 25,
#                 'PRIMER_OPT_TM': 60.0,
#                 'PRIMER_MIN_TM': 57.0,
#                 'PRIMER_MAX_TM': 63.0,
#                 'PRIMER_MIN_GC': 20.0,
#                 'PRIMER_MAX_GC': 80.0,
#                 'PRIMER_MAX_POLY_X': 100,
#                 'PRIMER_INTERNAL_MAX_POLY_X': 100,
#                 'PRIMER_SALT_MONOVALENT': 50.0,
#                 'PRIMER_DNA_CONC': 50.0,
#                 'PRIMER_MAX_NS_ACCEPTED': 0,
#                 'PRIMER_MAX_SELF_ANY': 12,
#                 'PRIMER_MAX_SELF_END': 8,
#                 'PRIMER_PAIR_MAX_COMPL_ANY': 12,
#                 'PRIMER_PAIR_MAX_COMPL_END': 8,
#             },
#             {
#                 'PRIMER_PRODUCT_SIZE_RANGE': [75,100,100,125,125,150,150,175,175,200,200,225],
#                 'SEQUENCE_ID': 'MH1000',
#                 'SEQUENCE_TEMPLATE': 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCTCTAGTAG',
#                 'SEQUENCE_INCLUDED_REGION': [36,342]
#             }
#         )
#         wrapper_res = wrappers.runP3Main(
#             {
#                 'PRIMER_OPT_SIZE': 20,
#                 'PRIMER_PICK_INTERNAL_OLIGO': 1,
#                 'PRIMER_INTERNAL_MAX_SELF_END': 8,
#                 'PRIMER_MIN_SIZE': 18,
#                 'PRIMER_MAX_SIZE': 25,
#                 'PRIMER_OPT_TM': 60.0,
#                 'PRIMER_MIN_TM': 57.0,
#                 'PRIMER_MAX_TM': 63.0,
#                 'PRIMER_MIN_GC': 20.0,
#                 'PRIMER_MAX_GC': 80.0,
#                 'PRIMER_MAX_POLY_X': 100,
#                 'PRIMER_INTERNAL_MAX_POLY_X': 100,
#                 'PRIMER_SALT_MONOVALENT': 50.0,
#                 'PRIMER_DNA_CONC': 50.0,
#                 'PRIMER_MAX_NS_ACCEPTED': 0,
#                 'PRIMER_MAX_SELF_ANY': 12,
#                 'PRIMER_MAX_SELF_END': 8,
#                 'PRIMER_PAIR_MAX_COMPL_ANY': 12,
#                 'PRIMER_PAIR_MAX_COMPL_END': 8,
#                 'PRIMER_PRODUCT_SIZE_RANGE': '75-100 100-125 125-150 150-175 175-200 200-225',
#                 'SEQUENCE_ID': 'MH1000',
#                 'SEQUENCE_TEMPLATE': 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCTCTAGTAG',
#                 'SEQUENCE_INCLUDED_REGION': '36,342'
#             }
#         )
#         for k, v in binding_res.items():
#             print('{}: {}, {}'.format(k, v, wrapper_res.get(k)))
#             try:
#                 self.assertEqual(str(wrapper_res.get(k)), v)
#             except AssertionError:
#                 print('Key: {}, Wrapper output: {}, Binding output: {}'.format(k, v, binding_res.get(k)))
#                 raise

if __name__ == '__main__':
    # unittest.main()
    for x in range(100):
        seq1 = ''.join([random.choice('ATGC') for _ in
                             range(random.randint(20, 59))])
        seq2 = ''.join([random.choice('ATGC') for _ in
                             range(random.randint(20, 59))])
        mv_conc = random.uniform(1, 200)
        dv_conc = random.uniform(0, 40)
        dntp_conc = random.uniform(0, 20)
        dna_conc = random.uniform(0, 200)
        temp_c = random.randint(10, 70)
        max_loop = random.randint(10, 30)
        binding_res = bindings.calcHairpin(
                            seq=seq1,
                            mv_conc=mv_conc,
                            dv_conc=dv_conc,
                            dntp_conc=dntp_conc,
                            dna_conc=dna_conc,
                            temp_c=temp_c,
                            max_loop=max_loop)
        wrapper_res = wrappers.calcHairpin(
                            seq=seq1,
                            mv_conc=mv_conc,
                            dv_conc=dv_conc,
                            dntp_conc=dntp_conc,
                            dna_conc=dna_conc,
                            temp_c=temp_c,
                            max_loop=max_loop)
        if wrapper_res:
            print(int(binding_res.tm), int(wrapper_res.tm))
            if not int(binding_res.tm) == int(wrapper_res.tm):
                print(binding_res, '\n', wrapper_res)
