# Copyright (C) 2023. Ben Pruitt & Nick Conway
# See LICENSE for full GPLv2 license.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''
examples.orthogonalprimers

A simple script to find a complementary oligo set.

NOTE: More advanced desing might use windowed-hamming distance between reference
genome and do more off target checks or look at energy in addition to
melting temperature
'''
import os.path as op
import random
import sys
from typing import (
    Set,
    Tuple,
)

try:
    from primer3 import thermoanalysis
except ModuleNotFoundError:
    PACKAGE_DIR = op.dirname(op.dirname(op.abspath(__file__)))
    sys.path.append(PACKAGE_DIR)
    from primer3 import thermoanalysis


def rev_complement(seq: str) -> str:
    '''Compute the reverse complement of a sequence

    Args:
        seq: 5' to 3' sequence reverse complement

    Returns:
        Reverse complement string
    '''
    lut = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
    }
    return ''.join(lut[x] for x in seq[::-1])


def check_3p_prime_end(seq: str) -> bool:
    '''Check the 3 prime end of an oligo for GC content

    Args:
        seq: sequence to check

    Returns:
        True if GC content at the 3' end is 3 or less else False
    '''
    end3p = seq[-5:]
    return True if end3p.count('G') + end3p.count('C') < 4 else False


def check_gc_content(seq: str, seq_len: int) -> bool:
    '''Check a sequence for total GC content

    Args:
        seq: sequence to check
        seq_len: length sequence string

    Returns:
        True if GC content is less than 50 else False
    '''
    gc_total = seq.count('G') + seq.count('C')
    return True if (gc_total / seq_len) < 0.51 else False


def search_for_30_mers() -> Tuple[int, Set[str]]:
    ''' A simple and inexhaustive search for a set 20-mers that do not
    heterodimerize among themselves or form hairpins or homodimers.

    Returns:
        Tuple of the form::

            <number of candidates found>, <set of oligo strings>
    '''

    # Requirements
    oligo_screen_set_size_limit = 1000  # screen search size
    desired_total_seqs = 200            # break out of loop after 200 are found
    oligo_size = 30                     # we are looking for 30 mers

    # Temperature limits and cutoffs in Celcius
    tm_lim_lo_c = 60
    tm_lim_hi_c = 65
    tm_hairpin_homodimer_cutoff_c = 40
    tm_offtarget_cutoff_c = 40

    # Define analysis parameters and instance
    thermo_params = {
        'mv_conc': 50,      # Monovalent cation concentration in mM
        'dv_conc': 1.5,     # Divalent cation concentration in mM
        'dntp_conc': 0.2,   # dNTP concentration in mM
        'dna_conc': 200,    # DNA concentration in nM
    }
    ta_obj = thermoanalysis.ThermoAnalysis()
    ta_obj.set_thermo_args(**thermo_params)

    # Count the number of qualifying sequences
    found_seqs = 0
    candidate_seq_list = []

    for _ in range(oligo_screen_set_size_limit):
        # Generate a random candidate
        candidate_seq = ''.join([
            random.choice('ATGC') for _ in
            range(oligo_size)
        ])
        if not (
            check_gc_content(candidate_seq, 20) and
            check_3p_prime_end(candidate_seq)
        ):
            continue
        cand_tm_c = ta_obj.calc_tm(candidate_seq)
        cand_hrp_tm = ta_obj.calc_hairpin(candidate_seq).tm
        cand_homo_tm = ta_obj.calc_hairpin(candidate_seq).tm
        if (
            (tm_lim_lo_c < cand_tm_c < tm_lim_hi_c) and
            (cand_hrp_tm > tm_hairpin_homodimer_cutoff_c) and
            (cand_homo_tm > tm_hairpin_homodimer_cutoff_c)
        ):
            candidate_seq_list.append(candidate_seq)
            found_seqs += 1
            if found_seqs > desired_total_seqs:
                break

    # Now check heterodimers and create output greedily
    good_set: Set[str] = set()
    for i, seq1 in enumerate(candidate_seq_list):
        seq1_rev = rev_complement(seq1)

        # NOTE: nothing is done with the bad set but this could be used in a
        # different search pattern
        bad_set: Set[str] = set()

        for _, seq2 in enumerate(candidate_seq_list[i + 1:]):
            tr_obj = ta_obj.calc_heterodimer(seq1, seq2)
            if tr_obj.tm > tm_offtarget_cutoff_c:
                bad_set.add(seq2)
            tr_obj = ta_obj.calc_heterodimer(seq1_rev, seq2)
            if tr_obj.tm > tm_offtarget_cutoff_c:
                bad_set.add(seq2)
        good_set.add(seq1)

    return len(candidate_seq_list), good_set


if __name__ == '__main__':
    number_of_candidates, result_set = search_for_30_mers()
    print(
        f'Found {number_of_candidates} candidates, with total primer '
        f'set size of {len(result_set)}',
    )
    print(f'The {result_set=}')
