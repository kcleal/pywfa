/*
 *                             The MIT License
 *
 * Wavefront Alignment Algorithms
 * Copyright (c) 2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 * This file is part of Wavefront Alignment Algorithms.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * PROJECT: Wavefront Alignment Algorithms
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "utils/commons.h"
#include "wavefront_bialign.h"
#include "wavefront_unialign.h"
#include "wavefront_bialigner.h"

#include "wavefront_compute.h"
#include "wavefront_compute_affine.h"
#include "wavefront_compute_affine2p.h"
#include "wavefront_compute_edit.h"
#include "wavefront_compute_linear.h"
#include "wavefront_extend.h"
#include "wavefront_plot.h"
#include "wavefront_debug.h"

/*
 * Config
 */
#define WF_BIALIGN_FALLBACK_MIN_SCORE  250
#define WF_BIALIGN_FALLBACK_MIN_LENGTH 100
#define WF_BIALIGN_RECOVERY_MIN_SCORE  500

/*
 * Debug
 */
void wavefront_bialign_debug(
    wf_bialign_breakpoint_t* const breakpoint,
    const int align_level) {
  // Parameters
  const int breakpoint_h = WAVEFRONT_H(breakpoint->k_forward,breakpoint->offset_forward);
  const int breakpoint_v = WAVEFRONT_V(breakpoint->k_forward,breakpoint->offset_forward);
  // Prinf debug info
  fprintf(stderr,"[WFA::BiAlign][Recursion=%d] ",align_level);
  int i; for (i=0;i<align_level;++i) fprintf(stderr,"   ");
  fprintf(stderr,"Breakpoint at (h,v,score,comp) = (%d,%d,%d,",
      breakpoint_h,breakpoint_v,breakpoint->score);
  switch (breakpoint->component) {
    case affine2p_matrix_M:  fprintf(stderr,"M");  break;
    case affine2p_matrix_I1: fprintf(stderr,"I1"); break;
    case affine2p_matrix_I2: fprintf(stderr,"I2"); break;
    case affine2p_matrix_D1: fprintf(stderr,"D1"); break;
    case affine2p_matrix_D2: fprintf(stderr,"D2"); break;
    default: fprintf(stderr,"?"); break;
  }
  fprintf(stderr,")\n");
}
/*
 * Init
 */
void wavefront_bialign_init(
    wavefront_bialigner_t* const bialigner,
    const distance_metric_t distance_metric,
    alignment_form_t* const form,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end,
    const int align_level,
    const int verbose) {
  // Parameters
  wavefront_aligner_t* const wf_forward = bialigner->wf_forward;
  wavefront_aligner_t* const wf_reverse = bialigner->wf_reverse;
  // Configure WF-compute function
  switch (distance_metric) {
    case indel:
    case edit:
      bialigner->wf_align_compute = &wavefront_compute_edit;
      break;
    case gap_linear:
      bialigner->wf_align_compute = &wavefront_compute_linear;
      break;
    case gap_affine:
      bialigner->wf_align_compute = &wavefront_compute_affine;
      break;
    case gap_affine_2p:
      bialigner->wf_align_compute = &wavefront_compute_affine2p;
      break;
    default:
      fprintf(stderr,"[WFA] Distance function not implemented\n");
      exit(1);
      break;
  }
  // Initialize wavefront-aligner (forward)
  alignment_span_t span_forward =
      (form->pattern_begin_free > 0 || form->text_begin_free > 0) ?
          alignment_endsfree : alignment_end2end;
  alignment_form_t form_forward = {
      .span = span_forward,
      .pattern_begin_free = form->pattern_begin_free,
      .pattern_end_free = 0,
      .text_begin_free = form->text_begin_free,
      .text_end_free = 0,
  };
  wf_forward->alignment_form = form_forward;
  wf_forward->component_begin = component_begin;
  wf_forward->component_end = component_end;
  wavefront_aligner_init(wf_forward,align_level);
  // Initialize wavefront-aligner (reverse)
  alignment_span_t span_reverse =
      (form->pattern_end_free > 0 || form->text_end_free > 0) ?
          alignment_endsfree : alignment_end2end;
  alignment_form_t form_reverse = {
      .span = span_reverse,
      .pattern_begin_free = form->pattern_end_free,
      .pattern_end_free = 0,
      .text_begin_free = form->text_end_free,
      .text_end_free = 0,
  };
  wf_reverse->alignment_form = form_reverse;
  wf_reverse->component_begin = component_end;
  wf_reverse->component_end = component_begin;
  wavefront_aligner_init(wf_reverse,align_level);
  // Plot
  const bool plot_enabled = (wf_forward->plot != NULL);
  if (plot_enabled) {
    wavefront_plot(wf_forward,0,align_level);
    wavefront_plot(wf_reverse,0,align_level);
  }
  // DEBUG
  if (verbose >= 2) {
    wavefront_debug_begin(wf_forward);
    wavefront_debug_begin(wf_reverse);
  }
}
/*
 * Bidirectional Alignment (base cases)
 */
int wavefront_bialign_base(
    wavefront_aligner_t* const wf_aligner,
    alignment_form_t* const form,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end,
    const int align_level) {
  // Parameters
  wavefront_aligner_t* const wf_base = wf_aligner->bialigner->wf_base;
  const int verbose = wf_base->system.verbose;
  // Configure
  wf_base->alignment_form = *form;
  wavefront_unialign_init(wf_base,component_begin,component_end);
  // DEBUG
  if (verbose >= 2) wavefront_debug_begin(wf_base);
  // Wavefront align sequences
  wavefront_unialign(wf_base);
  // DEBUG
  if (verbose >= 2) {
    wavefront_debug_end(wf_base);
    wavefront_debug_check_correct(wf_base);
  }
  // Append CIGAR
  cigar_append_forward(wf_aligner->cigar,wf_base->cigar);
  // Set status and return
  const int align_status = wf_base->align_status.status;
  if (align_status == WF_STATUS_ALG_COMPLETED) {
    return WF_STATUS_OK;
  } else {
    return WF_STATUS_UNATTAINABLE;
  }
}
/*
 * Bidirectional check breakpoints
 */
void wavefront_bialign_breakpoint_indel2indel(
    wavefront_aligner_t* const wf_aligner,
    const bool breakpoint_forward,
    const int score_0,
    const int score_1,
    wavefront_t* const dwf_0,
    wavefront_t* const dwf_1,
    const affine2p_matrix_type component,
    wf_bialign_breakpoint_t* const breakpoint) {
  // Parameters
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int text_length = sequences->text_length;
  const int pattern_length = sequences->pattern_length;
  const int gap_open =
      (component==affine2p_matrix_I1 || component==affine2p_matrix_D1) ?
      wf_aligner->penalties.gap_opening1 : wf_aligner->penalties.gap_opening2;
  // Check wavefronts overlapping
  const int lo_0 = dwf_0->lo;
  const int hi_0 = dwf_0->hi;
  const int lo_1 = WAVEFRONT_K_INVERSE(dwf_1->hi,pattern_length,text_length);
  const int hi_1 = WAVEFRONT_K_INVERSE(dwf_1->lo,pattern_length,text_length);
  if (hi_1 < lo_0 || hi_0 < lo_1) return;
  // Compute overlapping interval
  const int min_hi = MIN(hi_0,hi_1);
  const int max_lo = MAX(lo_0,lo_1);
  int k_0;
  for (k_0=max_lo;k_0<=min_hi;k_0++) {
    const int k_1 = WAVEFRONT_K_INVERSE(k_0,pattern_length,text_length);
    // Fetch offsets
    const wf_offset_t doffset_0 = dwf_0->offsets[k_0];
    const wf_offset_t doffset_1 = dwf_1->offsets[k_1];
    const int dh_0 = WAVEFRONT_H(k_0,doffset_0);
    const int dh_1 = WAVEFRONT_H(k_1,doffset_1);
    // Check breakpoint d2d
    if (dh_0 + dh_1 >= text_length && score_0 + score_1 - gap_open < breakpoint->score) {
      if (breakpoint_forward) {
        // Check out-of-bounds coordinates
        const int v = WAVEFRONT_V(k_0,dh_0);
        const int h = WAVEFRONT_H(k_0,dh_0);
        if (v > pattern_length || h > text_length) continue;
        // Set breakpoint
        breakpoint->score_forward = score_0;
        breakpoint->score_reverse = score_1;
        breakpoint->k_forward = k_0;
        breakpoint->k_reverse = k_1;
        breakpoint->offset_forward = dh_0;
        breakpoint->offset_reverse = dh_1;
      } else {
        // Check out-of-bounds coordinates
        const int v = WAVEFRONT_V(k_1,dh_1);
        const int h = WAVEFRONT_H(k_1,dh_1);
        if (v > pattern_length || h > text_length) continue;
        // Set breakpoint
        breakpoint->score_forward = score_1;
        breakpoint->score_reverse = score_0;
        breakpoint->k_forward = k_1;
        breakpoint->k_reverse = k_0;
        breakpoint->offset_forward = dh_1;
        breakpoint->offset_reverse = dh_0;
      }
      breakpoint->score = score_0 + score_1 - gap_open;
      breakpoint->component = component;
      // wavefront_bialign_debug(breakpoint,-1); // DEBUG
      // No need to keep searching
      return;
    }
  }
}
void wavefront_bialign_breakpoint_m2m(
    wavefront_aligner_t* const wf_aligner,
    const bool breakpoint_forward,
    const int score_0,
    const int score_1,
    wavefront_t* const mwf_0,
    wavefront_t* const mwf_1,
    wf_bialign_breakpoint_t* const breakpoint) {
  // Parameters
  wavefront_sequences_t* const sequences = &wf_aligner->sequences;
  const int text_length = sequences->text_length;
  const int pattern_length = sequences->pattern_length;
  // Check wavefronts overlapping
  const int lo_0 = mwf_0->lo;
  const int hi_0 = mwf_0->hi;
  const int lo_1 = WAVEFRONT_K_INVERSE(mwf_1->hi,pattern_length,text_length);
  const int hi_1 = WAVEFRONT_K_INVERSE(mwf_1->lo,pattern_length,text_length);
  if (hi_1 < lo_0 || hi_0 < lo_1) return;
  // Compute overlapping interval
  const int min_hi = MIN(hi_0,hi_1);
  const int max_lo = MAX(lo_0,lo_1);
  int k_0;
  for (k_0=max_lo;k_0<=min_hi;k_0++) {
    const int k_1 = WAVEFRONT_K_INVERSE(k_0,pattern_length,text_length);
    // Fetch offsets
    const wf_offset_t moffset_0 = mwf_0->offsets[k_0];
    const wf_offset_t moffset_1 = mwf_1->offsets[k_1];
    const int mh_0 = WAVEFRONT_H(k_0,moffset_0);
    const int mh_1 = WAVEFRONT_H(k_1,moffset_1);
    // Check breakpoint m2m
    if (mh_0 + mh_1 >= text_length && score_0 + score_1 < breakpoint->score) {
      if (breakpoint_forward) {
        breakpoint->score_forward = score_0;
        breakpoint->score_reverse = score_1;
        breakpoint->k_forward = k_0;
        breakpoint->k_reverse = k_1;
        breakpoint->offset_forward = moffset_0;
        breakpoint->offset_reverse = moffset_1;
      } else {
        breakpoint->score_forward = score_1;
        breakpoint->score_reverse = score_0;
        breakpoint->k_forward = k_1;
        breakpoint->k_reverse = k_0;
        breakpoint->offset_forward = moffset_1;
        breakpoint->offset_reverse = moffset_0;
      }
      breakpoint->score = score_0 + score_1;
      breakpoint->component = affine2p_matrix_M;
      // wavefront_bialign_debug(breakpoint,-1); // DEBUG
      // No need to keep searching
      return;
    }
  }
}
/*
 * Bidirectional find overlaps
 */
void wavefront_bialign_overlap(
    wavefront_aligner_t* const wf_aligner_0,
    wavefront_aligner_t* const wf_aligner_1,
    const int score_0,
    const int score_1,
    const bool breakpoint_forward,
    wf_bialign_breakpoint_t* const breakpoint) {
  // Parameters
  const int max_score_scope = wf_aligner_0->wf_components.max_score_scope;
  const distance_metric_t distance_metric = wf_aligner_0->penalties.distance_metric;
  const int gap_opening1 = wf_aligner_0->penalties.gap_opening1;
  const int gap_opening2 = wf_aligner_0->penalties.gap_opening2;
  // Fetch wavefronts-0
  const int score_mod_0 = score_0 % max_score_scope;
  wavefront_t* const mwf_0 = wf_aligner_0->wf_components.mwavefronts[score_mod_0];
  if (mwf_0 == NULL) return;
  wavefront_t* d1wf_0 = NULL, *i1wf_0 = NULL;
  if (distance_metric >= gap_affine) {
    d1wf_0 = wf_aligner_0->wf_components.d1wavefronts[score_mod_0];
    i1wf_0 = wf_aligner_0->wf_components.i1wavefronts[score_mod_0];
  }
  wavefront_t* d2wf_0 = NULL, *i2wf_0 = NULL;
  if (distance_metric == gap_affine_2p) {
    d2wf_0 = wf_aligner_0->wf_components.d2wavefronts[score_mod_0];
    i2wf_0 = wf_aligner_0->wf_components.i2wavefronts[score_mod_0];
  }
  // Traverse all scores-1
  int i;
  for (i=0;i<max_score_scope;++i) {
    // Compute score
    const int score_i = score_1 - i;
    if (score_i < 0) break;
    const int score_mod_i = score_i % max_score_scope;
    // Check I2/D2-breakpoints (gap_affine_2p)
    if (distance_metric == gap_affine_2p) {
      if (score_0 + score_i - gap_opening2 >= breakpoint->score) continue;
      // Check breakpoint d2d
      wavefront_t* const d2wf_1 = wf_aligner_1->wf_components.d2wavefronts[score_mod_i];
      if (d2wf_0 != NULL && d2wf_1 != NULL) {
        wavefront_bialign_breakpoint_indel2indel(
            wf_aligner_0,breakpoint_forward,score_0,score_i,
            d2wf_0,d2wf_1,affine2p_matrix_D2,breakpoint);
      }
      // Check breakpoint i2i
      wavefront_t* const i2wf_1 = wf_aligner_1->wf_components.i2wavefronts[score_mod_i];
      if (i2wf_0 != NULL && i2wf_1 != NULL) {
        wavefront_bialign_breakpoint_indel2indel(
            wf_aligner_0,breakpoint_forward,score_0,score_i,
            i2wf_0,i2wf_1,affine2p_matrix_I2,breakpoint);
      }
    }
    // Check I1/D1-breakpoints (gap_affine)
    if (distance_metric >= gap_affine) {
      if (score_0 + score_i - gap_opening1 >= breakpoint->score) continue;
      // Check breakpoint d2d
      wavefront_t* const d1wf_1 = wf_aligner_1->wf_components.d1wavefronts[score_mod_i];
      if (d1wf_0 != NULL && d1wf_1 != NULL) {
        wavefront_bialign_breakpoint_indel2indel(
            wf_aligner_0,breakpoint_forward,score_0,score_i,
            d1wf_0,d1wf_1,affine2p_matrix_D1,breakpoint);
      }
      // Check breakpoint i2i
      wavefront_t* const i1wf_1 = wf_aligner_1->wf_components.i1wavefronts[score_mod_i];
      if (i1wf_0 != NULL && i1wf_1 != NULL) {
        wavefront_bialign_breakpoint_indel2indel(
            wf_aligner_0,breakpoint_forward,score_0,score_i,
            i1wf_0,i1wf_1,affine2p_matrix_I1,breakpoint);
      }
    }
    // Check M-breakpoints (indel, edit, gap-linear)
    if (score_0 + score_i >= breakpoint->score) continue;
    wavefront_t* const mwf_1 = wf_aligner_1->wf_components.mwavefronts[score_mod_i];
    if (mwf_1 != NULL) {
      wavefront_bialign_breakpoint_m2m(
          wf_aligner_0,breakpoint_forward,
          score_0,score_i,mwf_0,mwf_1,breakpoint);
    }
  }
}
/*
 * Bidirectional breakpoint detection
 */
int wavefront_bialign_overlap_gopen_adjust(
    wavefront_aligner_t* const wf_aligner,
    const distance_metric_t distance_metric) {
  switch (distance_metric) {
    case gap_affine:
      return wf_aligner->penalties.gap_opening1;
    case gap_affine_2p:
      return MAX(wf_aligner->penalties.gap_opening1,wf_aligner->penalties.gap_opening2);
    case indel:
    case edit:
    case gap_linear:
    default:
      return 0;
  }
}
int wavefront_bialign_find_breakpoint(
    wavefront_bialigner_t* const bialigner,
    const distance_metric_t distance_metric,
    alignment_form_t* const form,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end,
    wf_bialign_breakpoint_t* const breakpoint,
    const int align_level) {
  // Parameters
  wavefront_aligner_t* const wf_forward = bialigner->wf_forward;
  wavefront_aligner_t* const wf_reverse = bialigner->wf_reverse;
  alignment_system_t* const system = &wf_forward->system;
  const bool plot_enabled = (wf_forward->plot != NULL);
  const int verbose = system->verbose;
  // Init bialignment
  wavefront_bialign_init(bialigner,distance_metric,form,component_begin,component_end,align_level,verbose);
  // Sequences
  wavefront_sequences_t* const sequences = &wf_forward->sequences;
  const int text_length = sequences->text_length;
  const int pattern_length = sequences->pattern_length;
  // Operators
  void (*wf_align_compute)(wavefront_aligner_t* const,const int) = bialigner->wf_align_compute;
  // Parameters
  const int max_alignment_steps = wf_forward->system.max_alignment_steps;
  const int max_antidiagonal = DPMATRIX_ANTIDIAGONAL(pattern_length,text_length) - 1; // Note: Even removing -1
  int score_forward = 0, score_reverse = 0, forward_max_ak = 0, reverse_max_ak = 0;
  bool reachability_quit;
  // Prepare and perform first bialignment step
  breakpoint->score = INT_MAX;
  reachability_quit = wavefront_extend_end2end_max(wf_forward,score_forward,&forward_max_ak);
  if (reachability_quit) return wf_forward->align_status.status;
  reachability_quit = wavefront_extend_end2end_max(wf_reverse,score_reverse,&reverse_max_ak);
  if (reachability_quit) return wf_reverse->align_status.status;
  // Compute wavefronts of increasing score until both wavefronts overlap
  int max_ak = 0;
  bool last_wf_forward = false;
  while (true) {
    // Check close-to-collision
    if (forward_max_ak + reverse_max_ak >= max_antidiagonal) break;
    /*
     * Compute next wavefront (Forward)
     */
    ++score_forward;
    (*wf_align_compute)(wf_forward,score_forward);
    if (plot_enabled) wavefront_plot(wf_forward,score_forward,align_level); // Plot
    // Extend
    reachability_quit = wavefront_extend_end2end_max(wf_forward,score_forward,&max_ak);
    if (forward_max_ak < max_ak) forward_max_ak = max_ak;
    last_wf_forward = true;
    // Check end-reached and close-to-collision
    if (reachability_quit) return wf_forward->align_status.status;
    if (forward_max_ak + reverse_max_ak >= max_antidiagonal) break;
    /*
     * Compute next wavefront (Reverse)
     */
    ++score_reverse;
    (*wf_align_compute)(wf_reverse,score_reverse);
    if (plot_enabled) wavefront_plot(wf_reverse,score_reverse,align_level); // Plot
    // Extend
    reachability_quit = wavefront_extend_end2end_max(wf_reverse,score_reverse,&max_ak);
    if (reverse_max_ak < max_ak) reverse_max_ak = max_ak;
    last_wf_forward = false;
    // Check end-reached and max-steps-reached
    if (reachability_quit) return wf_reverse->align_status.status;
    if (score_reverse + score_forward >= max_alignment_steps) return WF_STATUS_MAX_STEPS_REACHED;
    // DEBUG
    if (verbose >= 3 && score_forward % system->probe_interval_global == 0) {
      wavefront_unialign_print_status(stderr,wf_forward,score_forward);
    }
  }
  // Advance until overlap is found
  const int max_score_scope = wf_forward->wf_components.max_score_scope;
  const int gap_opening = wavefront_bialign_overlap_gopen_adjust(wf_forward,distance_metric);
  while (true) {
    if (last_wf_forward) {
      // Check overlapping wavefronts
      const int min_score_reverse = (score_reverse > max_score_scope-1) ? score_reverse - (max_score_scope-1) : 0;
      if (score_forward + min_score_reverse - gap_opening >= breakpoint->score) break; // Done!
      wavefront_bialign_overlap(wf_forward,wf_reverse,score_forward,score_reverse,true,breakpoint);
      /*
       * Compute next wavefront (Reverse)
       */
      ++score_reverse;
      (*wf_align_compute)(wf_reverse,score_reverse);
      if (plot_enabled) wavefront_plot(wf_reverse,score_reverse,align_level); // Plot
      // Extend & check end-reached
      reachability_quit = wavefront_extend_end2end(wf_reverse,score_reverse);
      if (reachability_quit) return wf_reverse->align_status.status;
    }
    // Check overlapping wavefronts
    const int min_score_forward = (score_forward > max_score_scope-1) ? score_forward - (max_score_scope-1) : 0;
    if (min_score_forward + score_reverse - gap_opening >= breakpoint->score) break; // Done!
    wavefront_bialign_overlap(wf_reverse,wf_forward,score_reverse,score_forward,false,breakpoint);
    /*
     * Compute next wavefront (Forward)
     */
    ++score_forward;
    (*wf_align_compute)(wf_forward,score_forward);
    if (plot_enabled) wavefront_plot(wf_forward,score_forward,align_level); // Plot
    // Extend & check end-reached/max-steps-reached
    reachability_quit = wavefront_extend_end2end(wf_forward,score_forward);
    if (reachability_quit) return wf_forward->align_status.status;
    if (score_reverse + score_forward >= max_alignment_steps) return WF_STATUS_MAX_STEPS_REACHED;
    // Enable always
    last_wf_forward = true;
  }
  // Return OK
  return WF_STATUS_OK;
}
int wavefront_bialign_find_breakpoint_exception(
    wavefront_aligner_t* const wf_aligner,
    alignment_form_t* const form,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end,
    const int align_level,
    const int align_status) {
  // Breakpoint was not found, check end reached
  if (align_status == WF_STATUS_END_REACHED) {
    wavefront_aligner_t* const wf_forward = wf_aligner->bialigner->wf_forward;
    wavefront_aligner_t* const wf_reverse = wf_aligner->bialigner->wf_reverse;
    // Retrieve score when end was reached
    int score_reached;
    if (wf_forward->align_status.status == WF_STATUS_END_REACHED) {
      score_reached = wf_forward->align_status.score;
    } else {
      score_reached = wf_reverse->align_status.score;
    }
    // Fallback if possible
    if (score_reached <= WF_BIALIGN_RECOVERY_MIN_SCORE) {
      return wavefront_bialign_base(wf_aligner,form,component_begin,component_end,align_level);
    } else {
      return WF_STATUS_END_UNREACHABLE; // To no avail
    }
  } else { // Other unrecoverable conditions
    return align_status;
  }
}
/*
 * Bidirectional Alignment
 */
void wavefront_bialign_init_half_0(
    alignment_form_t* const global_form,
    alignment_form_t* const half_form) {
  // Align half_0
  const alignment_span_t span_0 =
      (global_form->pattern_begin_free > 0 ||
       global_form->text_begin_free > 0) ?
           alignment_endsfree : alignment_end2end;
  half_form->span = span_0;
  half_form->extension = false;
  half_form->pattern_begin_free = global_form->pattern_begin_free;
  half_form->pattern_end_free = 0;
  half_form->text_begin_free = global_form->text_begin_free;
  half_form->text_end_free = 0;
}
void wavefront_bialign_init_half_1(
    alignment_form_t* const global_form,
    alignment_form_t* const half_form) {
  // Align half_0
  const alignment_span_t span_1 =
      (global_form->pattern_begin_free > 0 ||
       global_form->text_begin_free > 0) ?
           alignment_endsfree : alignment_end2end;
  half_form->span = span_1;
  half_form->extension = false;
  half_form->pattern_begin_free = 0;
  half_form->pattern_end_free = global_form->pattern_end_free;
  half_form->text_begin_free = 0;
  half_form->text_end_free = global_form->text_end_free;
}
int wavefront_bialign_alignment(
    wavefront_aligner_t* const wf_aligner,
    alignment_form_t* const form,
    const affine2p_matrix_type component_begin,
    const affine2p_matrix_type component_end,
    const int score_remaining,
    const int align_level) {
  // Parameters
  wavefront_sequences_t* const sequences = &wf_aligner->bialigner->wf_forward->sequences;
  const int pattern_begin = sequences->pattern_begin;
  const int pattern_end = sequences->pattern_begin + sequences->pattern_length;
  const int text_begin = sequences->text_begin;
  const int text_end = sequences->text_begin + sequences->text_length;
  const int pattern_length = pattern_end - pattern_begin;
  const int text_length = text_end - text_begin;
  // Trivial cases
  if (text_length == 0) {
    cigar_append_deletion(wf_aligner->cigar,pattern_length);
    return WF_STATUS_OK;
  } else if (pattern_length == 0) {
    cigar_append_insertion(wf_aligner->cigar,text_length);
    return WF_STATUS_OK;
  } else if (score_remaining <= WF_BIALIGN_FALLBACK_MIN_SCORE) {
    // Fall back to regular WFA
    return wavefront_bialign_base(wf_aligner,form,
        component_begin,component_end,align_level);
  }
  // Find breakpoint in the alignment
  wf_bialign_breakpoint_t breakpoint;
  int align_status = wavefront_bialign_find_breakpoint(
      wf_aligner->bialigner,wf_aligner->penalties.distance_metric,
      form,component_begin,component_end,&breakpoint,align_level);
  // DEBUG
  if (wf_aligner->system.verbose >= 2) {
    wf_aligner->bialigner->wf_forward->align_status.status = align_status;
    wf_aligner->bialigner->wf_reverse->align_status.status = align_status;
    wavefront_debug_end(wf_aligner->bialigner->wf_forward);
    wavefront_debug_end(wf_aligner->bialigner->wf_reverse);
  }
  // Check status
  if (align_status != WF_STATUS_OK) {
    return wavefront_bialign_find_breakpoint_exception(
        wf_aligner,form,component_begin,component_end,align_level,align_status);
  }
  // Breakpoint found
  const int breakpoint_h = WAVEFRONT_H(breakpoint.k_forward,breakpoint.offset_forward);
  const int breakpoint_v = WAVEFRONT_V(breakpoint.k_forward,breakpoint.offset_forward);
  // DEBUG
  if (wf_aligner->system.verbose >= 3) wavefront_bialign_debug(&breakpoint,align_level);
  // Align half_0
  alignment_form_t form_0;
  wavefront_bialigner_set_sequences_bounds(wf_aligner->bialigner,
      pattern_begin,pattern_begin+breakpoint_v,
      text_begin,text_begin+breakpoint_h);
  wavefront_bialign_init_half_0(form,&form_0);
  align_status = wavefront_bialign_alignment(wf_aligner,
      &form_0,component_begin,breakpoint.component,
      breakpoint.score_forward,align_level+1);
  if (align_status != WF_STATUS_OK) return align_status;
  // Align half_1
  alignment_form_t form_1;
  wavefront_bialigner_set_sequences_bounds(wf_aligner->bialigner,
      pattern_begin+breakpoint_v,pattern_end,
      text_begin+breakpoint_h,text_end);
  wavefront_bialign_init_half_1(form,&form_1);
  align_status = wavefront_bialign_alignment(wf_aligner,
      &form_1,breakpoint.component,component_end,
      breakpoint.score_reverse,align_level+1);
  if (align_status != WF_STATUS_OK) return align_status;
  // Set score (Strictly speaking, only needed at level-0)
  if (align_level == 0) {
    cigar_t* const cigar = wf_aligner->cigar;
    cigar->score = wavefront_compute_classic_score(wf_aligner,pattern_length,text_length,breakpoint.score);
    cigar->end_v = pattern_length;
    cigar->end_h = text_length;
  }
  return WF_STATUS_OK; // All good
}
/*
 * Bidirectional Score-only
 */
int wavefront_bialign_compute_score(
    wavefront_aligner_t* const wf_aligner) {
  // Parameters
  wavefront_aligner_t* const wf_forward = wf_aligner->bialigner->wf_forward;
  wavefront_aligner_t* const wf_reverse = wf_aligner->bialigner->wf_reverse;
  wavefront_sequences_t* const sequences = &wf_forward->sequences;
  const int text_length = sequences->text_length;
  const int pattern_length = sequences->pattern_length;
  // Clear cigar
  cigar_clear(wf_aligner->cigar);
  // Find breakpoint in the alignment
  wf_bialign_breakpoint_t breakpoint;
  const int align_status = wavefront_bialign_find_breakpoint(wf_aligner->bialigner,
      wf_aligner->penalties.distance_metric,&wf_aligner->alignment_form,
      affine_matrix_M,affine_matrix_M,&breakpoint,0);
  // DEBUG
  if (wf_aligner->system.verbose >= 2) {
    wavefront_debug_end(wf_forward);
    wavefront_debug_end(wf_reverse);
  }
  // Check status
  cigar_t* const cigar = wf_aligner->cigar;
  if (align_status == WF_STATUS_OK || align_status == WF_STATUS_END_REACHED) {
    if (align_status == WF_STATUS_END_REACHED) {
      breakpoint.score = (wf_forward->align_status.status == WF_STATUS_END_REACHED) ?
          wf_forward->align_status.score : wf_reverse->align_status.score;
    }
    // Set status & score
    cigar->score = wavefront_compute_classic_score(wf_aligner,pattern_length,text_length,breakpoint.score);
    cigar->end_v = pattern_length;
    cigar->end_h = text_length;
    // Return OK
    return WF_STATUS_OK;
  } else {
    // Other cases
    return align_status;
  }
}
/*
 * Bidirectional dispatcher
 */
void wavefront_bialign(
    wavefront_aligner_t* const wf_aligner) {
  // Select scope
  int align_status;
  if (wf_aligner->alignment_scope == compute_score) {
    align_status = wavefront_bialign_compute_score(wf_aligner);
  } else {
    // Resize CIGAR
    wavefront_sequences_t* const sequences = &wf_aligner->bialigner->wf_forward->sequences;
    const int text_length = sequences->text_length;
    const int pattern_length = sequences->pattern_length;
    cigar_resize(wf_aligner->cigar,2*(pattern_length+text_length)); // Resize & clear
    // Bidirectional alignment
    const bool min_length = MAX(pattern_length,text_length) <= WF_BIALIGN_FALLBACK_MIN_LENGTH;
    align_status = wavefront_bialign_alignment(wf_aligner,
        &wf_aligner->alignment_form,
        affine_matrix_M,affine_matrix_M,
        min_length ? 0 : INT_MAX,0);
  }
  // Check status
  if (align_status == WF_STATUS_OK) {
    wf_aligner->align_status.status = WF_STATUS_ALG_COMPLETED;
  } else if (align_status == WF_STATUS_MAX_STEPS_REACHED || align_status == WF_STATUS_OOM) {
    wf_aligner->align_status.status = align_status;
  } else { // Other cases
    wf_aligner->align_status.status = WF_STATUS_UNATTAINABLE;
  }
}
