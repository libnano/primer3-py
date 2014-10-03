'''
primer3_test.py
~~~~~~~~~~~~~~~

Unit tests for the primer3-py package.

'''

from __future__ import print_function

import os
import random
import resource
import unittest

from time import sleep

from primer3 import bindings, wrappers, simulatedBindings


LOCAL_DIR = os.path.dirname(os.path.realpath(__file__))


def _getMemUsage():
    """ Get current process memory usage in bytes """
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


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
        for _ in range(25):
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
            self.assertEqual(int(binding_res.tm), int(wrapper_res.tm))

    def test_calcHomodimer(self):
        for _ in range(25):
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
            self.assertEqual(int(binding_res.tm), int(wrapper_res.tm))


    def test_calcHeterodimer(self):
        for _ in range(25):
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
            self.assertEqual(int(binding_res.tm), int(wrapper_res.tm))

    def test_correctionMethods(self):
        self.randArgs()
        for sc_method in ['schildkraut', 'santalucia', 'owczarzy']:
            for tm_method in ['breslauer', 'santalucia']:
                binding_tm = bindings.calcTm(
                                    seq=self.seq1,
                                    mv_conc=int(self.mv_conc),
                                    dv_conc=self.dv_conc,
                                    dntp_conc=self.dntp_conc,
                                    dna_conc=int(self.dna_conc),
                                    tm_method=tm_method,
                                    salt_corrections_method=sc_method)
                wrapper_tm = wrappers.calcTm(
                                    seq=self.seq1,
                                    mv_conc=int(self.mv_conc),
                                    dv_conc=self.dv_conc,
                                    dntp_conc=self.dntp_conc,
                                    dna_conc=int(self.dna_conc),
                                    tm_method=tm_method,
                                    salt_corrections_method=sc_method)
                self.assertEqual(int(binding_tm), int(wrapper_tm))
        self.assertRaises(ValueError, bindings.calcTm,
                                    seq=self.seq1,
                                    mv_conc=int(self.mv_conc),
                                    dv_conc=self.dv_conc,
                                    dntp_conc=self.dntp_conc,
                                    dna_conc=int(self.dna_conc),
                                    tm_method='not_a_tm_method')

    def test_memoryLeaks(self):
        sm = _getMemUsage()
        for x in range(1000):
            self.randArgs()
            bindings.calcHeterodimer(
                seq1=self.seq1,
                seq2=self.seq2,
                mv_conc=self.mv_conc,
                dv_conc=self.dv_conc,
                dntp_conc=self.dntp_conc,
                dna_conc=self.dna_conc,
                temp_c=self.temp_c,
                max_loop=self.max_loop)
        sleep(0.1)  # Pause for any GC
        em = _getMemUsage()
        print('\n\tMemory usage before 1k runs of calcHeterodimer: ', sm)
        print('\tMemory usage after 1k runs of calcHeterodimer:  ', em)
        print('\t\t\t\t\tDifference: \t', em-sm)
        if em-sm > 500:
            raise AssertionError('Memory usage increase after 1k runs of \n\t'
                                 'calcHeterodimer > 500 bytes -- potential \n\t'
                                 'memory leak (mem increase: {})'.format(em-sm))

class TestDesignBindings(unittest.TestCase):

    def _compareResults(self, binding_res, simulated_binding_res, 
                        verbose=False):
        keys_in_sim = set(simulated_binding_res)
        keys_in_binding = set(binding_res)

        if keys_in_sim - keys_in_binding:
            if verbose:
                print('\n\n\nIn wrapper simulation result but missing'
                      ' from binding:')
                fmt = '{:<30} {:<50}'
                print(fmt.format('Output Key', 'SimBinding Result'))
                print('-'*80)
                for k in sorted(keys_in_sim - keys_in_binding):
                    print(fmt.format(k, repr(simulated_binding_res[k])))

        if keys_in_binding - keys_in_sim:
            if verbose:
                print('\n\n\nIn binding result but missing from wrapper '
                      'simulation:')
                fmt = '{:<30} {:<50}'
                print(fmt.format('Output Key', 'Binding Result'))
                print('-'*80)
                for k in sorted(keys_in_binding - keys_in_sim):
                    print(fmt.format(k, repr(binding_res[k])))

        allowable_relative_difference = 0.005
        disagree = set(k for k in keys_in_binding & keys_in_sim
                       if simulated_binding_res[k] != binding_res[k]
                       and (isinstance(binding_res[k], float)
                       and (binding_res[k] != 0)
                            and abs((binding_res[k] - simulated_binding_res[k])
                                    / binding_res[k]) > 
                                    allowable_relative_difference) )
        if disagree:
            fmt = '{:<30} {:<25} {:<25}'
            disagreements = '\n'.join([fmt.format(k,
                                        repr(simulated_binding_res[k]),
                                        repr(binding_res[k])) for k in
                                        sorted(disagree)])
            if verbose:
                print('\n\n\nResults disagree:')
                print(fmt.format('Output Key', 'SimBinding Result', 
                                 'Binding Result'))
                print('-'*80)
                print(disagreements)
            raise RuntimeError('Results between binary output and binding '
                               'output disagree: \n' + disagreements)

        else:
            if verbose:
                print('\n\n\nAll the results in common ({}) agree to within '
                      '{:.2%}'.format(len(keys_in_binding & keys_in_sim),
                                      allowable_relative_difference))

    def _convertBoulderInput(self, boulder_str):
        ''' Convert a boulder IO-style input dictionary into bindings /
        simulated-bindings-friendly dictionaries.
        '''
        boulder_dicts = wrappers._parseMultiRecordBoulderIO(boulder_str)
        input_dicts = []
        for bd in boulder_dicts:
            converted_input = [simulatedBindings.unwrap(arg) for arg in 
                               bd.items()]
            global_args = dict(filter(lambda arg: "SEQUENCE" not in arg[0], 
                                      converted_input))
            seq_args = dict(filter(lambda arg: "SEQUENCE" in arg[0], 
                                    converted_input))
            input_dicts.append((global_args, seq_args))
        return input_dicts

    def testHuman(self):
        sequence_template = 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCTCTAGTAG'
        quality_list = [random.randint(20,90) for i in range(len(sequence_template))]
        seq_args = {
            'SEQUENCE_ID': 'MH1000',
            'SEQUENCE_TEMPLATE': sequence_template,
            'SEQUENCE_QUALITY': quality_list,
            'SEQUENCE_INCLUDED_REGION': [36,342]
        }
        global_args = {
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
            'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],[150,175],[175,200],[200,225]],
        }
        simulated_binding_res = simulatedBindings.designPrimers(seq_args, global_args)
        binding_res = bindings.designPrimers(seq_args, global_args)
        wrapper_res = wrappers.designPrimers(
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
                'PRIMER_PRODUCT_SIZE_RANGE': '75-100 100-125 125-150 150-175 175-200 200-225',
                'SEQUENCE_ID': 'MH1000',
                'SEQUENCE_TEMPLATE': sequence_template,
                'SEQUENCE_QUALITY': ' '.join(map(str, quality_list)),
                'SEQUENCE_INCLUDED_REGION': '36,342'
            }
        )
        print('\n\n\n{:<30} {:<25} {:<25} {:<25}'.format('Output Key', 'Wrapper Result', 'SimBinding Result', 'Binding Result'))
        print('-'*80)
        for result_field in sorted(binding_res.keys()):
            print('{:<30} {:<25} {:<25} {:<25}'.format(result_field,
                                                       repr(wrapper_res.get(result_field)),
                                                       repr(simulated_binding_res.get(result_field)),
                                                       repr(binding_res[result_field])))


    def testCompareSim(self):
        sequence_template = 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCTCTAGTAG'
        quality_list = [random.randint(20,90) for i in range(len(sequence_template))]
        seq_args = {
            'SEQUENCE_ID': 'MH1000',
            'SEQUENCE_TEMPLATE': sequence_template,
            'SEQUENCE_QUALITY': quality_list,
            'SEQUENCE_INCLUDED_REGION': [36,342]
        }
        global_args = {
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
            'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],[150,175],[175,200],[200,225]],
        }
        simulated_binding_res = simulatedBindings.designPrimers(seq_args, global_args)
        binding_res = bindings.designPrimers(seq_args, global_args)
        self._compareResults(binding_res, simulated_binding_res, verbose=True)

    def test_fileBased(self):
        test_file_roots = [
            'primer_internal',
            'primer_internal1',
            'primer_ok_regions',
            'primer_tm_lc_masking',
            'primer_task',
            'primer_renewed_tasks',
        ]
        for fn_root in test_file_roots:
            base_fp = os.path.join(LOCAL_DIR, 'test_files', fn_root)
            input_fp = base_fp + '_input'
            # output_fp = base_fp + '_output'
            with open(input_fp) as input_fd:
                input_raw = input_fd.read()
            input_dicts = self._convertBoulderInput(input_raw)
            # with open(output_fp) as output_fd:
            #     output_dict = output_fd.read()
            print('Testing {}'.format(fn_root))
            current_global_args = {}
            i = 0
            for global_args, seq_args in input_dicts:
                if i == 2:
                    break
                current_global_args.update(global_args)
                simulated_binding_res = simulatedBindings.designPrimers(
                                            seq_args, current_global_args)
                binding_res = bindings.designPrimers(seq_args, 
                                                     current_global_args)
                self._compareResults(binding_res, simulated_binding_res)



    def test_memoryLeaks(self):
        sm = _getMemUsage()
        for x in range(100):
            # bindings.runP3Design()
            bindings.designPrimers(
                {
                    'SEQUENCE_ID': 'MH1000',
                    'SEQUENCE_TEMPLATE': 'GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCCTACATTTTAGCATCAGTGAGTACAGCATGCTTACTGGAAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTGCAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCATGGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTGAGCAAGTNAAAAAAATGTTTGGAAGTGTTACTTTAGCAATGGCAAGAATGATAGTATGGAATAGATTGGCAGAATGAAGGCAAAATGATTAGACATATTGCATTAAGGTAAAAAATGATAACTGAAGAATTATGTGCCACACTTATTAATAAGAAAGAATATGTGAACCTTGCAGATGTTTCCCTCTAGTAG',
                    'SEQUENCE_INCLUDED_REGION': [36,342]
                },
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
                    'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],[150,175],[175,200],[200,225]],
                })
        sleep(0.1)  # Pause for any GC
        em = _getMemUsage()
        print('\n\tMemory usage before 1k runs of designPrimers: ', sm)
        print('\tMemory usage after 1k runs of designPrimers:  ', em)
        print('\t\t\t\t\tDifference: \t', em-sm)
        if em-sm > 1000:
            raise AssertionError('Memory usage increase after 1k runs of \n\t'
                                 'designPrimers > 1000 bytes -- potential \n\t'
                                 'memory leak (mem increase: {})'.format(em-sm))

if __name__ == '__main__':
    import sys

    tl = unittest.TestLoader()
    lowLevelSuite = tl.loadTestsFromTestCase(TestLowLevelBindings)
    res1 = unittest.TextTestRunner(verbosity=2).run(lowLevelSuite)
    designSuite = tl.loadTestsFromTestCase(TestDesignBindings)
    res2 = unittest.TextTestRunner(verbosity=2).run(designSuite)

    success = res1.wasSuccessful() and res2.wasSuccessful()

    sys.exit(int(not success))  # Exit 0 on success, 1 on failure
