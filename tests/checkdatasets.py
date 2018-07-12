import json
import subprocess
import sys
import glob
import os.path
import os
import sys

import click
dirname = os.path.dirname
extra = [dirname(dirname(os.path.abspath(__file__)))]
sys.path = extra + sys.path

from primer3.thermoanalysis import (
    ThermoAnalysis,
    ThermoResult
)
from primer3 import __version__ as p3version


def precision(x):
    '''Round the precision of floats to make comparisons exact
    '''
    return round(x, 2)

def todict_tr(x):
    return {
        'structure_found': x.structure_found,
        'tm': precision(x.tm),
        'dg': precision(x.dg),
        'dh': precision(x.dh),
        'ds': precision(x.ds)
    }

def todict_ta(x):
    return {
        'mv_conc': x.mv_conc,
        'dv_conc': x.dv_conc,
        'dntp_conc': x.dntp_conc,
        'dna_conc': x.dna_conc,
        'temp_c':  x.temp,
        'max_loop': x.max_loop,
        # 'temp_only': x.temp_only,
        # 'debug': x.thalargs.debug,
        'max_nn_length': x.max_nn_length,
        'tm_method': x.tm_method,
        'salt_correction_method': x.salt_correction_method
    }

def gitHash():
    try:
        with open(os.devnull, 'wb') as devnull:
            hash_str = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], stderr=devnull)
        return hash_str.strip().decode('utf8')
    except:
        return 'NA'
# end def

def createDataSet():
    seq1 = 'GGGGCCCCCCAAATTTTTT'
    seq2 = 'AAAAAATTTGGCCCCAAA'

    ta = ThermoAnalysis()
    tests = [
        ['calcHeterodimer', seq1, seq2],
        ['calcHairpin', seq1],
        ['calcHomodimer', seq1],
        ['calcTm', seq1]
    ]
    cases = []
    for test in tests:
        f = getattr(ta, test[0])
        res = f(*test[1:])
        # print(test[0])
        if isinstance(res, ThermoResult):
            cases.append([test, todict_tr(res)])
        else:
            cases.append([test, res])
    # end for
    out = {
            'thermoanalysis': todict_ta(ta),
            'cases': cases
    }
    # import pprint
    # pprint.pprint(out)
    return out
# end def

def dumpDataSet(x):
    import time
    date_str = time.strftime("%Y-%m-%d-%H-%M-%S", time.gmtime())
    git_hash = gitHash()
    py_version = '%s.%s.%s' % sys.version_info[0:3]
    platform = sys.platform
    filename = 'primer3_precision_%s_%s_%s_%s.json' % (date_str, git_hash, platform, py_version)
    with open(filename, 'w') as fd:
        print('dumping')
        json.dump(x, fd)
    return filename
# end def

def readDataSet():
    files = glob.glob('*.json')
    if len(files) > 0:
        the_file = files[0]
        file_name = os.path.splitext(the_file)[0]
        git_hash, platform, py_version = file_name.split('_')[3:6]
        print("\t- Checking reference file: %s" % the_file)
        print("\t\tgit hash: %s\n\t\tPlatform: %s\n\t\tPython: %s" % (git_hash, platform, py_version))
        with open(the_file, 'r') as fd:
            out = json.load(fd)
        return out
    else:
        raise OSError("No json file found")
# end def

def compareSets(reference, test_set):
    for k1, v in reference.items():
        if k1 not in test_set:
            return False
        test_val = test_set[k1]
        if k1 == 'thermoanalysis':
            for k2, item in v.items():
                if item != test_val[k2]:
                    return False
        elif k1 == 'cases':
            for item_ref, item_test in zip(v, test_val):
                ref_args, ref_result = item_ref
                test_args, test_result = item_test
                if ref_args != test_args:
                    print(ref_args, test_args)
                    return False
                if ref_result != test_result:
                    print(ref_result, test_result)
                    return False

    return True
# end def

@click.command()
@click.option('--generate', default=False, help='Generate a set')
def runner(generate):
    create_dict = createDataSet()
    if generate:
        filename = dumpDataSet(create_dict)
        print("> Dumped data set: %s" % (filename))
    else:
        print("> Comparing reference and current revisions")
        print("\t- Current git hash: %s" % (gitHash()))
        print("\t- Current primer3 version: %s" % (p3version))
        read_dict = readDataSet()
        print("\t- Is the reference set equal to the current revisions set?")
        are_equal = compareSets(read_dict, create_dict)
        print('\t\t%s' % (str(are_equal)))

if __name__ == '__main__':
    runner()
