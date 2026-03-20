#ifndef MM2_H
#define MM2_H

#include "mahit.h"
#include "aln2tr.h"

std::vector<tr_paf_rec_t> run_minimap2_all_alignments(
    const char *ref_fn,
    const char *query_fn,
    const char *preset,
    int best_n);

std::vector<tr_paf_rec_t> run_minimap2_in_memory(
    const char *ref_name, const char *ref_seq,
    const char *query_name, const char *query_seq,
    int best_n);

tr_p_aln_info process_reads_parallel(const std::string &filepath, int best_n, int num_threads);

tr_p_aln_info process_reads_serial(const std::string &filepath, int best_n);

#endif