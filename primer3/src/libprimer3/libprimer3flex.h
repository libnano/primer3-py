/* libprimer3/libprimer3flex.h

 Copyright (c) 1996,1997,1998,1999,2000,2001,2004,2006,2007,2009,2010,
               2011,2012,2016
 Whitehead Institute for Biomedical Research, Steve Rozen
 (http://purl.com/STEVEROZEN/), and Helen Skaletsky
 All rights reserved.

       This file is part of primer3 software suite.

       This software suite is is free software;
       you can redistribute it and/or modify it under the terms
       of the GNU General Public License as published by the Free
       Software Foundation; either version 2 of the License, or (at
       your option) any later version.

       This software is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

       You should have received a copy of the GNU General Public License
       along with this software (file gpl-2.0.txt in the source
       distribution); if not, write to the Free Software
       Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON A THEORY
 OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Modifications to libprimer.c code to create libprimer3flex.h.
Copyright (C) 2023. Ben Pruitt & Nick Conway;

See libprimer3flex.c for more info on libprimer3flex.h
*/
#ifndef _P3PY_LIBPRIMER3FLEX_H
#define _P3PY_LIBPRIMER3FLEX_H

#include "khash.h"
#include "libprimer3.h"

typedef struct dpal_arg_holder {
  dpal_args *local;
  dpal_args *end;
  dpal_args *local_end;
  dpal_args *local_ambig;
  dpal_args *local_end_ambig;
} dpal_arg_holder;

typedef struct thal_arg_holder {
  thal_args *any;
  thal_args *end1;
  thal_args *end2;
  thal_args *hairpin_th;
  int thermodynamic_alignment_length_error;         /* p3py update: formally static global */
  char* thermodynamic_alignment_length_error_msg;   /* p3py update: formally static global */
} thal_arg_holder;

/* Ininitalize khash.h `primer_pair_map` type `khash_t(primer_pair_map)` */
KHASH_MAP_INIT_INT(primer_pair_map, primer_pair* )

/* primer3-py update: 2023.03.21 moved to pairs and max_j_seen to
    function-local (choose_pair_or_triple) for thread safety
 */
typedef struct pairs_args_t {
     int *max_j_seen; /* The maximum value of j (loop index for forward primers)
                         that has been examined for every reverse primer
                         index (i) */
     khash_t(primer_pair_map) **pairs; /* hash map for storing found primer
                                          pairs by pair index, previously
                                          global */
} pairs_args_t;

/**************************/
/* Function declarations. */
/**************************/
static void
pr_set_default_global_args_1(p3_global_settings*);

static void
pr_set_default_global_args_2(p3_global_settings*);

static void
_adjust_seq_args(
    const p3_global_settings* pa,
    seq_args_t* sa,
    pr_append_str* nonfatal_err,
    pr_append_str* warning
);

static void
_optimize_ok_regions_list(
    const p3_global_settings* pa,
    seq_args_t* sa
);

static int
any_5_prime_ol_extension_has_problem(
    const primer_rec*
);

static int
p3_ol_is_uninitialized(const primer_rec*);

static int
fake_a_sequence(
    seq_args_t* sa,
    const p3_global_settings* pa
);

static int
_pr_data_control(
    const p3_global_settings*,
    const seq_args_t*,
    pr_append_str* glob_err,
    pr_append_str* nonfatal_err,
    pr_append_str* warning
);

static int
_pr_need_pair_template_mispriming(
    const p3_global_settings* pa
);
static int
_pr_need_pair_template_mispriming_thermod(
    const p3_global_settings* pa
);

static int
_pr_need_template_mispriming(
    const p3_global_settings*
);
static int
_pr_need_template_mispriming_thermod(
    const p3_global_settings*
);

static void
_pr_substr(
    const char*,
    int,
    int,
    char*
);

static void
_pr_shift_substr(
    const char*,
    int,
    int,
    int,
    char*
);

static int
_pr_violates_poly_x(
    const char* oligo_seq,
    int max_poly_x
);

static int
_check_and_adjust_intervals(
    seq_args_t* sa,
    int seq_len,
    int first_index,
    pr_append_str*  nonfatal_err,
    pr_append_str* warning
);

static int
_check_and_adjust_overlap_pos(
    seq_args_t* sa,
    int* list,
    int* count,
    const char* tag,
    int seq_len,
    int first_index,
    pr_append_str* nonfatal_err,
    pr_append_str* warning
);

static int
_check_and_adjust_1_interval(
    const char*,
    int num,
    interval_array_t,
    int,
    int first_index,
    pr_append_str* err,
    seq_args_t*,
    pr_append_str* warning,
    int
);

static void
sort_primer_array(
    oligo_array*
);

static void
add_pair(
    const primer_pair*,
    pair_array_t*
);

static double
align(
    const char*,
    const char*,
    const dpal_args *a
);

static double
align_thermod(
    const char*,
    const char*,
    const thal_args *a,
    thal_arg_holder* thal_arg_to_use
);

void
destroy_primer_sec_struct(
    primer_rec* p_rec
);

void
destroy_pair_sec_struct(
    primer_pair* ppair
);

void
save_overwrite_sec_struct(
    char** base,
    char* inp
);

void
recalc_secundary_structures(
    const p3_global_settings* pa,
    const seq_args_t* sa,
    const dpal_arg_holder* dpal_arg_to_use,
    const thal_arg_holder* thal_arg_to_use,
    const thal_arg_holder* thal_oligo_arg_to_use,
    p3retval* retval
);

void
recalc_primer_sec_struct(
    primer_rec*p_rec,
    const int rev,
    const p3_global_settings* pa,
    const seq_args_t* sa,
    const dpal_arg_holder* dpal_arg_to_use,
    const thal_arg_holder* thal_arg_to_use
);

void
recalc_pair_sec_struct(
    primer_pair* ppair,
    const p3_global_settings* pa,
    const seq_args_t* sa,
    const dpal_arg_holder* dpal_arg_to_use,
    const thal_arg_holder* thal_arg_to_use
);

static int
characterize_pair(
    p3retval* p,
    const p3_global_settings*,
    const seq_args_t*,
    int,
    int,
    int,
    primer_pair* ,
    const dpal_arg_holder*,
    const thal_arg_holder*,
    int update_stats
);

static void
choose_pair_or_triple(
    p3retval*,
    const p3_global_settings*,
    const seq_args_t*,
    const dpal_arg_holder*,
    const thal_arg_holder*,
    const thal_arg_holder*,
    pair_array_t*,
    pairs_args_t* pairs_args
);

static int
sequence_quality_is_ok(
    const p3_global_settings*,
    primer_rec*,
    oligo_type,
    const seq_args_t*,
    int,
    int,
    oligo_stats* global_oligo_stats,
    const args_for_one_oligo_or_primer*
);

static int
choose_internal_oligo(
    p3retval*,
    const primer_rec*,
    const primer_rec*,
    int*,
    int*,
    const seq_args_t*,
    const p3_global_settings*,
    const dpal_arg_holder*,
    const thal_arg_holder*,
    const thal_arg_holder*
);

void
compute_position_penalty(
    const p3_global_settings*,
    const seq_args_t*,
    primer_rec*,
    oligo_type
);

dpal_arg_holder*
create_dpal_arg_holder(void);

void
destroy_dpal_arg_holder(dpal_arg_holder* h);

thal_arg_holder*
create_thal_arg_holder(
    const args_for_one_oligo_or_primer *po_args
);

void
destroy_thal_arg_holder(
    thal_arg_holder* h
);

static p3retval*
create_p3retval(void);

static char
dna_to_upper(
    char*,
    int
);

static int
find_stop_codon(
    const char*,
    int,
    int
);

static void
gc_and_n_content(
    int,
    int,
    const char*,
    primer_rec*
);

static int
make_detection_primer_lists(
    p3retval*,
    const p3_global_settings*,
    const seq_args_t*,
    const dpal_arg_holder*,
    const thal_arg_holder*
);

static int
make_complete_primer_lists(
    p3retval* retval,
    const p3_global_settings* pa,
    const seq_args_t* sa,
    const dpal_arg_holder* dpal_arg_to_use,
    const thal_arg_holder* thal_arg_to_use,
    const thal_arg_holder* thal_oligo_arg_to_use
);

static int
add_primers_to_check(
    p3retval* retval,
    const p3_global_settings* pa,
    const seq_args_t* sa,
    const dpal_arg_holder* dpal_arg_to_use,
    const thal_arg_holder* thal_arg_to_use,
    const thal_arg_holder* thal_oligo_arg_to_use
);

static int
pick_sequencing_primer_list(
    p3retval* retval,
    const p3_global_settings* pa,
    const seq_args_t* sa,
    const dpal_arg_holder* dpal_arg_to_use,
    const thal_arg_holder* thal_arg_to_use
);

static int
make_internal_oligo_list(
    p3retval*,
    const p3_global_settings*,
    const seq_args_t*,
    const dpal_arg_holder*,
    const thal_arg_holder*
);

static int
pick_only_best_primer(
    const int,
    const int,
    oligo_array *,
    const p3_global_settings*,
    const seq_args_t*,
    const dpal_arg_holder*,
    const thal_arg_holder*,
    p3retval*
);

static int
pick_primer_range(
    const int,
    const int,
    int*,
    oligo_array *,
    const p3_global_settings*,
    const seq_args_t*,
    const dpal_arg_holder*,
    const thal_arg_holder*,
    p3retval* retval
);

static int
add_one_primer(
    const char*,
    int*,
    oligo_array *,
    const p3_global_settings*,
    const seq_args_t*,
    const dpal_arg_holder*,
    const thal_arg_holder*,
    p3retval*
);

static int
add_one_primer_by_position(
    int,
    int,
    int*,
    oligo_array *,
    const p3_global_settings*,
    const seq_args_t*,
    const dpal_arg_holder*,
    const thal_arg_holder*,
    p3retval*
);

static int
pick_primers_by_position(
    const int,
    const int,
    int*,
    oligo_array *,
    const p3_global_settings*,
    const seq_args_t*,
    const dpal_arg_holder*,
    const thal_arg_holder*,
    p3retval*
);

static double
obj_fn(
    const p3_global_settings*,
    primer_pair*
);

static int
left_oligo_in_pair_overlaps_used_oligo(
    const primer_rec* left,
    const primer_pair* best_pair,
    int min_dist
);

static int
intl_oligo_in_pair_overlaps_used_oligo(
    const primer_rec* left,
    const primer_pair* best_pair,
    int min_dist
);

static int
right_oligo_in_pair_overlaps_used_oligo(
    const primer_rec* right,
    const primer_pair* best_pair,
    int min_dist
);

static int
oligo_overlaps_interval(
    int,
    int,
    const interval_array_t,
    int
);

static void
calc_and_check_oligo_features(
    const p3_global_settings* pa,
    primer_rec*,
    oligo_type,
    const dpal_arg_holder*,
    const thal_arg_holder*,
    const seq_args_t*, oligo_stats*,
    p3retval*,
    const char*
);

static void
pr_append(
    pr_append_str* ,
    const char*
);

static void
pr_append_new_chunk(
    pr_append_str* x,
    const char* s
);

static int
pair_spans_target(
    const primer_pair*,
    const seq_args_t*
);

static void
pr_append_w_sep(
    pr_append_str* ,
    const char*,
    const char*
);

static void*
pr_safe_malloc(size_t x);

static void*
pr_safe_realloc(
    void* p,
    size_t x
);

static int
compare_primer_pair(
    const void*,
    const void*
);

static int
primer_rec_comp(
    const void*,
    const void*
);
static int
print_list_header(
    FILE*,
    oligo_type,
    int,
    int,
    int
);
static int
print_oligo(
    FILE*,
    const seq_args_t*,
    int,
    const primer_rec*,
    oligo_type,
    int,
    int,
    int
);
static char* strstr_nocase(char*, char*);

static double
p_obj_fn(
    const p3_global_settings*,
    primer_rec*,
    int
 );

static void
oligo_compl(
    primer_rec*,
    const args_for_one_oligo_or_primer *po_args,
    oligo_stats*,
    const dpal_arg_holder*,
    const char* oligo_seq,
    const char* revc_oligo_seq
);

static void
oligo_compl_thermod(
    primer_rec*, /* TO DO 2012-06-01 -- update by removing last argument. */
    const args_for_one_oligo_or_primer *po_args,
    oligo_stats*,
    const thal_arg_holder*,
    const char* oligo_seq,
    const char* revc_oligo_seq
);

static void
oligo_hairpin(
    primer_rec*,
    const args_for_one_oligo_or_primer *po_args,
    oligo_stats*,
    const thal_arg_holder*,
    const char* oligo_seq
);

static void
oligo_compute_sequence_and_reverse(
    primer_rec*,
    const seq_args_t*,
    oligo_type,
    int*,
    int*,
    char*,
    char*
);

static void
oligo_repeat_library_mispriming(
    primer_rec*,
    const p3_global_settings*,
    const seq_args_t*,
    oligo_type,
    oligo_stats*,
    const dpal_arg_holder*,
    pr_append_str*
);

static void
oligo_template_mispriming(
    primer_rec*,
    const p3_global_settings*,
    const seq_args_t*,
    oligo_type,
    oligo_stats*,
    const dpal_args *,
    const thal_args *,
    thal_arg_holder* thal_arg_to_use
);

static int
pair_repeat_sim(
    primer_pair* ,
    const p3_global_settings*
);

static void
free_repeat_sim_score(p3retval*);

static void
free_primer_repeat_sim_score(primer_rec*);

/* edited by T. Koressaar for lowercase masking:  */
static int
is_lowercase_masked(
    int position,
    const char* sequence,
    primer_rec* h,
    oligo_stats*
);

static int
primer_must_match(
    const p3_global_settings* pa,
    primer_rec* h,
    oligo_stats* stats,
    const char* input_oligo_seq,
    char* match_three_prime,
    char* match_five_prime
);

static int
compare_nucleotides(
    const char a,
    const char b
);

static int
test_must_match_parameters(
    char* test
);

static void
set_retval_both_stop_codons(
    const seq_args_t* sa,
    p3retval* retval
);

/* Functions to set bitfield parameters for oligos (or primers) */
static void
bf_set_overlaps_target(
    primer_rec*,
    int
);

static int
bf_get_overlaps_target(
    const primer_rec*
);

static void
bf_set_overlaps_excl_region(
    primer_rec*,
    int
);

static int
bf_get_overlaps_excl_region(
    const primer_rec*
);

static void
bf_set_infinite_pos_penalty(
    primer_rec*,
    int
);

static int
bf_get_infinite_pos_penalty(
    const primer_rec*
);

/* Functions to record problems with oligos (or primers) */
static void
initialize_op(
    primer_rec*
);

static void
op_set_completely_written(
    primer_rec*
);

static void
op_set_must_match_err(
    primer_rec*
);

static void
op_set_too_many_ns(
    primer_rec*
);

static void
op_set_overlaps_target(
    primer_rec*
);

static void
op_set_high_gc_content(
    primer_rec*
);

static void
op_set_low_gc_content(
    primer_rec*
);

static void
op_set_high_tm(
    primer_rec*
);

static void
op_set_low_tm(
    primer_rec*
);

static void
op_set_high_bound(
    primer_rec*
);

static void
op_set_low_bound(
    primer_rec*
);

static void
op_set_overlaps_excluded_region(
    primer_rec*
);

static void
op_set_not_in_any_ok_region(
    primer_rec*
);

static void
op_set_high_self_any(
    primer_rec* oligo
);

static void
op_set_high_self_end(
    primer_rec* oligo
);

static void
op_set_high_hairpin_th(
    primer_rec* oligo
);

static void
op_set_no_gc_glamp(
    primer_rec*
);

static void
op_set_too_many_gc_at_end(
    primer_rec*
);

static void
op_set_high_end_stability(
    primer_rec*
);

static void
op_set_high_poly_x(
    primer_rec*
);

static void
op_set_low_sequence_quality(
    primer_rec*
);

static void
op_set_low_end_sequence_quality(
    primer_rec* oligo
);

static void
op_set_high_similarity_to_non_template_seq(
    primer_rec*
);

static void
op_set_high_similarity_to_multiple_template_sites(
    primer_rec*
);

static void
op_set_overlaps_masked_sequence(
    primer_rec*
);

static void
op_set_too_long(
    primer_rec*
);

static void
op_set_too_short(
    primer_rec*
);

static void
op_set_does_not_amplify_orf(
    primer_rec*
);
/*******************************************/
/* End functions to set problems in oligos */
/*******************************************/


#endif
