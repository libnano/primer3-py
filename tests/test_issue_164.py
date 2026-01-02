import primer3


class SeqObj:
    """Object with str representation."""

    def __init__(self, seq):
        self._seq = seq

    def __str__(self):
        return self._seq


def test_calc_tm_with_object():
    """Test strict analysis function (calc_tm) with a custom object."""
    obj = SeqObj('GTAAAACGACGGCCAGT')
    tm = primer3.calc_tm(obj)

    assert isinstance(tm, float)


def test_calc_end_stability_with_object():
    """analysis function (calc_end_stability) with custom objects."""
    obj = SeqObj('GTAAAACGACGGCCAGT')

    res = primer3.calc_end_stability(obj, obj)

    assert hasattr(res, 'tm')
    assert isinstance(res.tm, float)
