import wrappers

p3_args = {}

interval_list_tags = ['SEQUENCE_INCLUDED_REGION',
                      'SEQUENCE_TARGET',
                      'SEQUENCE_EXCLUDED_REGION',
                      'SEQUENCE_INTERNAL_EXCLUDED_REGION']

#  (semicolon separated list of integer "quadruples"; default empty)
semi_quad_tags = ['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST']


def wrapListWithFormat(v, fmt, sep = ' '):
    try:
        rv = fmt % tuple(v)
    except TypeError:
        rv = sep.join(fmt % tuple(x) for x in v)
    return rv

def wrapListOfQuads(v):
    rv = ' ; '.join(
        ','.join(x and str(x) or ''
                 for x in quad)
        for quad in v)
    return rv


def wrap(t):
    k,v = t
    if isinstance(v, list):
        if len(v) == 0:
            rv = ""
        else:
            if k in semi_quad_tags:
                rv = wrapListOfQuads(v)
            elif k in interval_list_tags:
                rv = wrapListWithFormat(v, '%d,%d')
            elif isinstance(v[0], list):
                rv = wrapListWithFormat(v, '%d-%d')
            elif isinstance(v[0], int):
                rv = ' '.join(map(str, v))
            else:
                rv = v
    else:
        rv = v
    return k, rv

def unwrap(t):
    k,v = t
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
                    rv = map(int, v.split())
                except ValueError:
                    rv = v
    return k,rv


def setGlobals(global_args, misprime_lib, mishyb_lib):
    p3_args.update(dict(map(wrap, global_args.iteritems())))

def setSeqArgs(seq_args):
    p3_args.update(dict(map(wrap, seq_args.iteritems())))

def convertResult(result):
    return dict(map(unwrap, result.iteritems()))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Design bindings ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


def designPrimers(seq_args, global_args=None, misprime_lib=None,
                  mishyb_lib=None):
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
    setGlobals(global_args, misprime_lib, mishyb_lib)
    setSeqArgs(seq_args)
    result = wrappers.designPrimers(p3_args)
    return convertResult(result)
