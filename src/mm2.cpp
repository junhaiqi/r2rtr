#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include "mm2.h"
#include <omp.h>
// #include "ma"

extern "C"
{
#include "minimap.h"
#include "bseq.h"
}

#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

// Run minimap2 with relaxed parameters to collect multiple candidate alignments
std::vector<tr_paf_rec_t> run_minimap2_all_alignments(
    const char *ref_fn,
    const char *query_fn,
    const char *preset,
    int best_n)
{
    std::vector<tr_paf_rec_t> results;

    // ================================
    // 0. Validate input parameters
    // ================================
    std::string preset_str(preset);
    if (preset_str != "ava-ont" && preset_str != "ava-pb")
    {
        std::cerr << "ERROR: Invalid preset '" << preset << "'. Only 'ava-ont' or 'ava-pb' are supported.\n";
        return results;
    }

    // ================================
    // 1. Initialize minimap2 options
    // ================================
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;

    // Apply the specified preset (ava-ont or ava-pb)
    mm_set_opt(preset, &iopt, &mopt);

    // Relax filtering parameters to retain more candidate alignments
    mopt.best_n = best_n;      // maximum number of secondary alignments (-N) set via parameter
    mopt.pri_ratio = 0.5;      // minimal score ratio for secondary alignments (-p)
    mopt.min_chain_score = 10; // minimal chaining score (-m)
    mopt.mid_occ = 0;          // disable high-frequency minimizer filtering (-f)

    // Enable secondary alignment reporting
    mopt.flag &= ~MM_F_NO_PRINT_2ND;

    // ================================
    // 2. Build reference index
    // ================================
    mm_idx_reader_t *reader = mm_idx_reader_open(ref_fn, &iopt, 0);
    if (!reader)
    {
        std::cerr << "ERROR: failed to open reference\n";
        return results;
    }

    mm_idx_t *idx = mm_idx_reader_read(reader, 1);
    mm_idx_reader_close(reader);

    if (!idx)
    {
        std::cerr << "ERROR: failed to build index\n";
        return results;
    }

    // ================================
    // 3. Open query sequences
    // ================================
    mm_bseq_file_t *fp = mm_bseq_open(query_fn);
    if (!fp)
    {
        std::cerr << "ERROR: failed to open query\n";
        mm_idx_destroy(idx);
        return results;
    }

    int n_seq = 0;
    mm_bseq1_t *seqs = nullptr;

    // Initialize thread buffer (Required by newer minimap2 API to prevent segfaults)
    mm_tbuf_t *tbuf = mm_tbuf_init();

    // ================================
    // 4. Perform mapping
    // ================================
    // Added '0' as the 3rd argument (with_qual = false) to match the new API signature
    while ((seqs = mm_bseq_read(fp, 1, 0, &n_seq)) != NULL)
    {
        for (int i = 0; i < n_seq; ++i)
        {
            int n_reg = 0;

            // Run alignment for one query sequence
            mm_reg1_t *regs = mm_map(
                idx,
                seqs[i].l_seq,
                seqs[i].seq,
                &n_reg,
                tbuf, // Passed the thread buffer here (5th argument)
                &mopt,
                NULL);

            // ================================
            // 5. Store alignment results
            // ================================
            for (int j = 0; j < n_reg; ++j)
            {
                mm_reg1_t &r = regs[j];

                tr_paf_rec_t rec;

                // Query information
                rec.ql = seqs[i].l_seq;
                rec.qs = r.qs;
                rec.qe = r.qe;

                // Target (reference) information
                rec.tl = idx->seq[r.rid].len;
                rec.ts = r.rs;
                rec.te = r.re;

                results.push_back(rec);
            }

            // Free alignment results
            free(regs);

            // Free sequence memory
            free(seqs[i].seq);
            free(seqs[i].name);
        }

        free(seqs);
    }

    // ================================
    // 6. Clean up
    // ================================
    // Destroy the thread buffer to prevent memory leaks
    mm_tbuf_destroy(tbuf);

    mm_bseq_close(fp);
    mm_idx_destroy(idx);

    return results;
}

std::vector<tr_paf_rec_t> run_minimap2_in_memory(
    const char *ref_name, const char *ref_seq,
    const char *query_name, const char *query_seq,
    int best_n)
{
    std::vector<tr_paf_rec_t> results;

    mm_idxopt_t iopt;
    mm_mapopt_t mopt;

   
    memset(&iopt, 0, sizeof(mm_idxopt_t));
    memset(&mopt, 0, sizeof(mm_mapopt_t));

    // preset = "";
    mm_set_opt(NULL, &iopt, &mopt);

    mopt.flag |= MM_F_ALL_CHAINS; 
    mopt.flag |= MM_F_NO_PRINT_2ND;
    mopt.flag &= ~MM_F_NO_PRINT_2ND; 

   
    mopt.seed = 42;
    mopt.best_n = (best_n > 0) ? best_n : 5;
    mopt.pri_ratio = 0.0f;;
    mopt.max_occ = 1000000;    
    mopt.mid_occ = 10000;
    mopt.mask_level = 1.0f;

    /***************************************/
    mopt.bw = 100;
    mopt.max_gap = 100; 
    mopt.max_gap_ref = -1; 
    mopt.min_cnt = 2;         
    mopt.min_chain_score = 20; 
    /***************************************/

    /***************************************/
    // iopt.k = 9;      
    // iopt.w = 5;       
    // mopt.zdrop = 999999;            
    // mopt.zdrop_inv = 999999;      
    // mopt.flag |= MM_F_CIGAR; 
    // mopt.min_chain_score = 15;
    // mopt.min_cnt = 1;      
    // mopt.bw = 500;
    // mopt.max_gap = 500; 
    /***************************************/
  

    const char *seq_arr[1] = {ref_seq};
    const char *name_arr[1] = {ref_name};

    mm_idx_t *idx = mm_idx_str(iopt.w, iopt.k, iopt.flag & MM_I_HPC, iopt.bucket_bits, 1, seq_arr, name_arr);
    if (!idx) return results;

    mm_mapopt_update(&mopt, idx);

    // if (mopt.mid_occ <= 0) {
    //     mopt.mid_occ = 1000; 
    // }
    mopt.mid_occ = 1000000;  

    mm_tbuf_t *tbuf = mm_tbuf_init();
    int n_reg = 0;
    int query_len = (int)strlen(query_seq);

    mm_reg1_t *regs = mm_map(idx, query_len, query_seq, &n_reg, tbuf, &mopt, query_name);

    for (int j = 0; j < n_reg; ++j) {
        mm_reg1_t &r = regs[j];
        tr_paf_rec_t rec;
        rec.ql = query_len;
        rec.qs = r.qs;
        rec.qe = r.qe;
        rec.tl = (int)idx->seq[r.rid].len;
        rec.ts = r.rs;
        rec.te = r.re;
        results.push_back(rec);
        if (r.p) free(r.p);
    }

    free(regs);
    mm_tbuf_destroy(tbuf);
    mm_idx_destroy(idx);

    return results;
}


// Structure to hold sequence region info
struct SeqRegion {
    std::string name;
    std::string seq;
};

// Extract prefix, suffix, and middle regions from a sequence
static std::vector<SeqRegion> extract_regions(const std::string& q_name, const std::string& q_seq, size_t region_len = 600) {
    std::vector<SeqRegion> regions;
    size_t q_len = q_seq.length();

    if (q_len == 0) return regions;

    // Prefix: beginning of sequence
    size_t prefix_len = std::min(region_len, q_len);
    regions.push_back({"prefix_" + q_name, q_seq.substr(0, prefix_len)});

    // Suffix: end of sequence (only if sequence is long enough)
    if (q_len > region_len) {
        regions.push_back({"suffix_" + q_name, q_seq.substr(q_len - region_len, region_len)});
    }

    // Middle: center of sequence (only if sequence is long enough)
    if (q_len > region_len * 2) {
        size_t start_pos = (q_len - region_len) / 2;
        regions.push_back({"mid_" + q_name, q_seq.substr(start_pos, region_len)});
    }

    return regions;
}

tr_p_aln_info process_reads_serial(
    const std::string &filepath,
    int best_n)
{
    tr_p_aln_info final_results;
    gzFile fp = gzopen(filepath.c_str(), "r");
    if (!fp)
    {
        std::cerr << "ERROR: Cannot open file " << filepath << "\n";
        return final_results;
    }

    kseq_t *seq = kseq_init(fp);

    while (kseq_read(seq) >= 0)
    {
        if (seq->seq.l == 0) continue;

        std::string q_name = seq->name.s;
        std::string q_seq = seq->seq.s;

        // Extract prefix, suffix, and middle regions
        auto regions = extract_regions(q_name, q_seq, 500);

        // Strategy: prioritize prefix, only try alternatives if prefix fails
        std::string best_r_name;
        std::vector<tr_paf_rec_t> best_alns;
        bool found_valid = false;

        // First, try prefix (most reliable for tandem repeats)
        for (const auto& region : regions) {
            if (region.name.find("prefix_") == 0) {
                auto alns = run_minimap2_in_memory(
                    region.name.c_str(), region.seq.c_str(),
                    q_name.c_str(), q_seq.c_str(),
                    best_n);

                if (alns.size() > 1) {
                    best_r_name = region.name;
                    best_alns = std::move(alns);
                    found_valid = true;
                    break;  // Prefix succeeded, use it
                }
            }
        }

        // If prefix failed, try suffix and middle
        if (!found_valid) {
            size_t best_aln_count = 0;
            for (const auto& region : regions) {
                if (region.name.find("prefix_") == 0) continue;  // Skip prefix, already tried

                auto alns = run_minimap2_in_memory(
                    region.name.c_str(), region.seq.c_str(),
                    q_name.c_str(), q_seq.c_str(),
                    best_n);

                if (alns.size() > best_aln_count) {
                    best_aln_count = alns.size();
                    best_r_name = region.name;
                    best_alns = std::move(alns);
                    found_valid = true;
                }
            }
        }

        // Store the alignment result if valid
        if (found_valid && best_alns.size() > 1)
        {
            final_results[q_name][best_r_name] = std::move(best_alns);
        }
    }

    kseq_destroy(seq);
    gzclose(fp);

    return final_results;
}

// ================================
// Parallel Processing Function
// ================================
tr_p_aln_info process_reads_parallel(const std::string &filepath, int best_n, int num_threads)
{
    tr_p_aln_info final_results;
    std::vector<SeqRecord> all_seqs;

    gzFile fp = gzopen(filepath.c_str(), "r");
    if (!fp) {
        std::cerr << "ERROR: Cannot open file " << filepath << "\n";
        return final_results;
    }
    kseq_t *seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        all_seqs.push_back({seq->name.s, seq->seq.s});
    }
    kseq_destroy(seq);
    gzclose(fp);

    omp_set_num_threads(num_threads);

#pragma omp parallel
    {
        tr_p_aln_info local_threads_map;
#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < all_seqs.size(); ++i)
        {
            const std::string &q_name = all_seqs[i].name;
            const std::string &q_seq = all_seqs[i].seq;

            // Extract prefix, suffix, and middle regions
            auto regions = extract_regions(q_name, q_seq, 500);

            // Strategy: prioritize prefix, only try alternatives if prefix fails
            std::string best_r_name;
            std::vector<tr_paf_rec_t> best_alns;
            bool found_valid = false;

            // First, try prefix (most reliable for tandem repeats)
            for (const auto& region : regions) {
                if (region.name.find("prefix_") == 0) {
                    auto alns = run_minimap2_in_memory(
                        region.name.c_str(), region.seq.c_str(),
                        q_name.c_str(), q_seq.c_str(),
                        best_n);

                    if (alns.size() > 1) {
                        best_r_name = region.name;
                        best_alns = std::move(alns);
                        found_valid = true;
                        break;  // Prefix succeeded, use it
                    }
                }
            }

            // If prefix failed, try suffix and middle
            if (!found_valid) {
                size_t best_aln_count = 0;
                for (const auto& region : regions) {
                    if (region.name.find("prefix_") == 0) continue;  // Skip prefix, already tried

                    auto alns = run_minimap2_in_memory(
                        region.name.c_str(), region.seq.c_str(),
                        q_name.c_str(), q_seq.c_str(),
                        best_n);

                    if (alns.size() > best_aln_count) {
                        best_aln_count = alns.size();
                        best_r_name = region.name;
                        best_alns = std::move(alns);
                        found_valid = true;
                    }
                }
            }

            // Store the alignment result if valid
            if (found_valid && best_alns.size() > 1)
            {
                local_threads_map[q_name][best_r_name] = std::move(best_alns);
            }
        }
#pragma omp critical(merge_final)
        {
            for (auto &kv : local_threads_map)
            {
                final_results[kv.first] = std::move(kv.second);
            }
        }
    }

    return final_results;
}
