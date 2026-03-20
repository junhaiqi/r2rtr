#include "log.h"
#include "paf.h"
#include "ketopt.h"
#include "mahit.h"
#include "aln2tr.h"
#include "mm2.h"
#include <stdio.h>
#include <string>
#include <iostream>

#define VERSION_MAJOR 1
#define VERSION_MINOR 1
#define VERSION_PATCH 0

int main(int argc, char *argv[])
{
    ketopt_t o = KETOPT_INIT;
    int32_t c;

    char *paf_path = "";
    int min_sp_num = 10;
    int min_tr_l = 10;
    int td_n = 1;
    int batch_size = 100;
    bool hc = 0;
    int len_cutoff = 0;
    int best_n = 100000;

    while ((c = ketopt(&o, argc, argv, 1, "f:n:l:t:b:c:r:N:", 0)) >= 0)
    {
        if (c == 'f')
        {
            // if (o.arg == 0)
            // {
            //     fprintf(stderr, "Error: -f requires an argument\n");
            //     return 1;
            // }
            paf_path = o.arg;
        }

        else if (c == 'r')
        {
            len_cutoff = atoi(o.arg);
        }

        else if (c == 'n')
        {
            min_sp_num = atoi(o.arg);
        }

        else if (c == 'l')
        {
            min_tr_l = atoi(o.arg);
        }

        else if (c == 'c')
        {
            hc = atoi(o.arg);
        }

        else if (c == 't')
        {
            td_n = atoi(o.arg);
        }

        else if (c == 'b')
        {
            batch_size = atoi(o.arg);
        }

        else if (c == 'N')
        {
            best_n = atoi(o.arg);
        }
    }

    if (argc - o.ind < 1)
    {
        fprintf(stderr, "Version %d.%d.%d\n", VERSION_MAJOR, VERSION_MINOR, VERSION_PATCH);
        // fprintf(stderr, "Usage: %s [Options:] <in.paf/paf.gz> -f <read.fa/fq/fq.gz> \n", argv[0]);
        fprintf(stderr, "Usage: %s [Options:] <read.fa/fq/fq.gz> \n", argv[0]);
        fprintf(stderr, "Options:\n");
        fprintf(stderr, "  -f STR     Specify the PAF file path from fasta/fastq "
                        "type [default = None]\n");
        fprintf(stderr, "  -r INT     Specify the minimum read length to infer tandem repeat [default = 0]\n");
        fprintf(stderr, "  -n INT     Specify the minimum number of supports for the tandem repeat length [default = 10]\n");
        fprintf(stderr, "  -l INT     Specify the inferred minimum tandem repeat length [default = 10]\n");
        fprintf(stderr, "  -c BOOL    Specify the accurate reads (e.g., HiFi, ONT-R10) [default = 0, current version's recommendation is 0]\n");
        fprintf(stderr, "  -t INT     Specify the thread number to use minimap2 or extract repeats using multi-threading [default = 1]\n");
        fprintf(stderr, "  -b INT     Specify the batch size to extract repeats using multi-threading [default = 100]\n");
        fprintf(stderr, "Minimap2 Options:\n");
        // fprintf(stderr, "  -X STR     Specify the minimap2 preset - ava-pb/ava-ont [default = ava-ont]\n");
        fprintf(stderr, "  -N INT     Specify a maximum of INT secondary alignments [default = 100000]\n");
        // fprintf(stderr, "  -o STR     Specify the output path [required parameters]\n");
        return 1;
    }

    if (argv[o.ind] == "")
    {
        fprintf(stderr, "Error: requires a paf file path\n");
        return 1;
    }

    // if (fa_path == "")
    // {
    //     fprintf(stderr, "Error: requires a read file path\n");
    //     return 1;
    // }

    start_main_timer();
    tr_p_aln_info tr_aln_lib;
    std::string message;
    if (paf_path != "")
    {
        message = std::string("Reading PAF file ") + paf_path + " and finding read pairs with multiple alignments...";
        log(message);
        f_tr_aln_info(paf_path, tr_aln_lib, len_cutoff);
        message = "A total of " + std::to_string(tr_aln_lib.size()) + " read pairs with multiple alignments.";
        log(message);
    }

    else
    {   
        message = std::string("Reading Fasta file ") + argv[o.ind] + " and establish self-alignment information...";
        log(message);
        tr_aln_lib = process_reads_parallel(argv[o.ind], best_n, td_n);
        message = "Completed.";
        log(message);

    }
    
    trs_h r_trs;
    score_h s_lib;
    message = "Inferring the tandem repeat length by read pairs with multiple alignments...";
    log(message);
    from_alns_to_trs(tr_aln_lib, r_trs, min_sp_num, min_tr_l, s_lib);
    message = "End of inference.";
    log(message);
    message = "Extracting tandem repeat units from reads...";
    log(message);

    if (td_n == 1)
        get_trs(argv[o.ind], r_trs, tr_aln_lib, s_lib, hc);
    else
        para_get_trs(argv[o.ind], r_trs, tr_aln_lib, s_lib, batch_size, td_n, hc);

    message = "End of extraction.";
    log(message);
    print_resource_usage("INFO");
    return 0;
}
