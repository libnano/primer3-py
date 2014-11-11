import primer3.wrappers as wrappers

p3_args = {}

interval_list_tags = set(['SEQUENCE_INCLUDED_REGION',
                      'SEQUENCE_TARGET',
                      'SEQUENCE_EXCLUDED_REGION',
                      'SEQUENCE_INTERNAL_EXCLUDED_REGION'])

size_range_list_tags = set(['PRIMER_PRODUCT_SIZE_RANGE'])

#  (semicolon separated list of integer "quadruples"; default empty)
semi_quad_tags = ['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST']


def wrapListOfQuads(v):
    ''' Wrap a list of ordered quads, potentially containing None's
    Produces a string 'list of semicolon separated list of integer "quadruples"'
    Used for SEQUENCE_PRIMER_PAIR_OK_REGION_LIST

    >>> wrapListOfQuads([1,2,3,4])
    '1,2,3,4'

    >>> wrapListOfQuads([(1,2,3,4)])
    '1,2,3,4'

    >>> wrapListOfQuads([[1,2,3,4]])
    '1,2,3,4'

    >>> wrapListOfQuads([[1,2,3,4], [5,6,7,8]])
    '1,2,3,4 ; 5,6,7,8'

    >>> wrapListOfQuads(((1,2,3,4), [5,6,7,8]))
    '1,2,3,4 ; 5,6,7,8'

    '''
    int_to_str = lambda i: str(i) if i > -1 else ''
    try:
        rv = ';'.join(
            ','.join(map(int_to_str, quad))
            for quad in v)
    except TypeError:
        rv = ','.join(x and str(x) or '' for x in v)
    return rv


def wrapListWithFormat(v, fmt, sep = ' '):
    try:
        rv = fmt % tuple(v)
    except TypeError:
        rv = sep.join(fmt % tuple(x) for x in v)
    return rv

def wrap(t):
    '''Convert a primer3 input in python-friendly bindings-style form
    to a string form for use by the process wrapper

    >>> wrap(('SEQUENCE_TEMPLATE', 'ATCG'))
    ('SEQUENCE_TEMPLATE', 'ATCG')

    >>> wrap(('SEQUENCE_QUALITY', range(5)))
    ('SEQUENCE_QUALITY', '0 1 2 3 4')

    >>> wrap(('SEQUENCE_EXCLUDED_REGION', (5,7)))
    ('SEQUENCE_EXCLUDED_REGION', '5,7')

    >>> wrap(('SEQUENCE_EXCLUDED_REGION', [(5,7), (11,13)]))
    ('SEQUENCE_EXCLUDED_REGION', '5,7 11,13')

    >>> wrap(('PRIMER_PRODUCT_SIZE_RANGE', (7,11)))
    ('PRIMER_PRODUCT_SIZE_RANGE', '7-11')
    '''
    k,v = t

    if isinstance(v, (list, tuple)):
        if len(v) == 0:
            rv = ""
        else:
            if k in semi_quad_tags:
                rv = wrapListOfQuads(v)
            elif k in interval_list_tags:
                rv = wrapListWithFormat(v, '%d,%d')
            elif k in size_range_list_tags:
                rv = wrapListWithFormat(v, '%d-%d')
            elif isinstance(v[0], (list, tuple)):
                rv = wrapListWithFormat(v, '%d-%d')
            elif isinstance(v[0], int):
                rv = ' '.join(map(str, v))
            else:
                rv = v
    else:
        rv = v
    return k, rv

def unwrap(t):
    '''convert a wrapper result into the intended form of the
    bindings result

    >>> unwrap(('A_FLOAT', '1.23'))
    ('A_FLOAT', 1.23)

    >>> unwrap(('AN_INT', '42'))
    ('AN_INT', 42)

    >>> unwrap(('A_SIZE_RANGE', '2020-2520'))
    ('A_SIZE_RANGE', (2020, 2520))

    >>> unwrap(('AN_INTERVAL', '7,11'))
    ('AN_INTERVAL', (7, 11))

    >>> unwrap(('MULTI_SIZE_RANGE', '1-2 3-5 7-11'))
    ('MULTI_SIZE_RANGE', ((1, 2), (3, 5), (7, 11)))

    >>> unwrap(('MULTIPLE_INTS', '2 3 5 7 11 13'))
    ('MULTIPLE_INTS', (2, 3, 5, 7, 11, 13))

    >>> unwrap(('A_SEQUENCE', 'ATCG'))
    ('A_SEQUENCE', 'ATCG')
    '''
    k,v = t
    rv = v
    condInt = lambda s: int(s) if s != '' else -1
    for lam in [lambda x: int(x),
                lambda x: float(x),
                lambda x: tuple(int(s) for s in x.split(' ')),
                lambda x: tuple(int(s) for s in x.split('-')),
                lambda x: tuple(int(s) for s in x.split(',')),
                lambda x: tuple(tuple(int(ss) for ss in s.split('-')) 
                                for s in x.split(' ')),
                lambda x: tuple(tuple(condInt(ss) for ss in s.split(',')) 
                                for s in x.split(' ')),
                lambda x: tuple(tuple(condInt(ss) for ss in 
                                s.strip().split(',')) for s in x.split(';'))
    ]:
        try:
            rv = lam(v)
        except:
            pass
        else:
            break
    """

    try:
        rv = int(v)
    except ValueError:
        try:
            rv = float(v)
        except ValueError:
            try:
                rv = tuple(int(s) for s in v.split(','))
            except ValueError:
                try:
                    rv = [int(x) for x in v.split()]
                except ValueError:
                    rv = v
"""
    return k,rv


def setGlobals(global_args, misprime_lib, mishyb_lib):
    p3_args.update(dict(wrap(v) for v in global_args.items()))

def setSeqArgs(seq_args):
    p3_args.update(dict(wrap(v) for v in  seq_args.items()))

def convertResult(result):
    return dict(unwrap(v) for v in result.items())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Design bindings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def designPrimers(seq_args, global_args=None, reset_args=True,
                  misprime_lib=None, mishyb_lib=None, input_log=None, 
                  output_log=None, err_log=None):
    ''' Run the Primer3 design process, with the same interface as the bindings,
    using the wrapped subprocess of primer3_core to do the work. 

    If the global args have been previously set (either by a pervious
    `designPrimers` call or by a `setGlobals` call), `designPrimers` may be
    called with seqArgs alone (as a means of optimization).

    Args:
        seq_args (dict)     :   Primer3 sequence/design args as per Primer3 docs

    Kwargs:
        global_args (dict)  :   Primer3 global args as per Primer3 docs
        misprime_lib (dict) :   `Sequence name: sequence` dictionary for
                                mispriming checks.
        mishyb_lib (dict)   :   `Sequence name: sequence` dictionary for
                                mishybridization checks.

    Returns:
        A dictionary of Primer3 results (should be identical to the expected
        BoulderIO output from primer3_main)

    '''
    if reset_args:
        p3_args.clear()
    setGlobals(global_args, misprime_lib, mishyb_lib)
    setSeqArgs(seq_args)
    result = wrappers.designPrimers(p3_args,
                                    input_log=input_log,
                                    output_log=output_log,
                                    err_log=err_log)
    return convertResult(result)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
