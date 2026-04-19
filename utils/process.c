/* The MIT License

   Copyright (c) 2025 Genome Research Ltd.

   Author: Petr Danecek <pd3@sanger.ac.uk>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.

 */

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <unistd.h>     // for isatty
#include <inttypes.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/bgzf.h>
#include <htslib/kseq.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcfutils.h>
#include <assert.h>
#include <errno.h>
#include <sys/time.h>
#include "mpileup2/mpileup.h"
#include "bcftools.h"
#include "regidx.h"


typedef struct
{
    char *sites_fname, *aln_fname, *fasta_fname;
    mpileup_t *mplp;
}
args_t;

const char *about(void)
{
    return "Prototype program for partial phasing from reads\n";
}
static const char *usage_text(void)
{
    return
        "\n"
        "About: Input preparation for partial phasing of variants\n"
        "Usage: bcftools +process [OPTIONS]\n"
        "\n"
        "Options:\n"
        "   -a, --aln FILE                BAM/CRAM file\n"
        "   -f, --fasta-ref FILE          Reference file in fasta format\n"
        "   -s, --sites FILE              A tab-delimited file name of sites to assess (chr,pos,ref,alt)\n"
        "\n";
}

static void destroy_data(args_t *args)
{
    if ( args->mplp ) mpileup_destroy(args->mplp);
    free(args);
}
static void init_data(args_t *args)
{
    args->mplp = mpileup_alloc();
    mpileup_set(args->mplp, MAX_DP_PER_SAMPLE, 250);
    mpileup_set(args->mplp, MIN_MQ, 0);
    mpileup_set(args->mplp, MAX_BQ, 60);
    mpileup_set(args->mplp, DELTA_BQ, 30);
    mpileup_set(args->mplp, MIN_REALN_FRAC, 0.05);
    mpileup_set(args->mplp, MIN_REALN_DP, 2);
    mpileup_set(args->mplp, MAX_REALN_DP, 250);
    mpileup_set(args->mplp, MAX_REALN_LEN, 500);
    mpileup_set(args->mplp, SKIP_ANY_SET, BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP);

    int ret = mpileup_set(args->mplp, REGIONS_FNAME, args->sites_fname);
    if ( ret ) error("Error: could not initialize site list %s\n",args->sites_fname);

    mpileup_set(args->mplp, LEGACY_MODE, 1);
    if ( mpileup_set(args->mplp, FASTA_REF, args->fasta_fname)!=0 ) error("Error: could not read the reference %s\n",args->fasta_fname);
    if ( mpileup_set(args->mplp, BAM, args->aln_fname)!=0 ) error("Error: could not reat %s\n",args->aln_fname);
    if ( mpileup_init(args->mplp)!=0 ) error("Error: could not initialize mpileup2\n");
}
static void process_data(args_t *args)
{
    int i,j,ret;
    int nsmpl  = mpileup_get_val(args->mplp,int,N_SAMPLES);
    int *n_plp = mpileup_get_val(args->mplp,int*,N_READS);

    // process the entire bam
    printf("#chr\tpos\tismpl\tbase\tBQ\tMQ\tread\n");
    while ( (ret=mpileup_next(args->mplp))==1 )
    {
        char *chr = mpileup_get_val(args->mplp,char*,CHROM);
        hts_pos_t pos = mpileup_get_val(args->mplp,hts_pos_t,POS);

        bam_pileup1_t **plp = mpileup_get_val(args->mplp,bam_pileup1_t**,LEGACY_PILEUP);
        //char *ref = mpileup_get(args->mplp,REF,&len);

        for (i=0; i<nsmpl; i++)
        {
            for (j=0; j<n_plp[i]; j++)
            {
                const bam_pileup1_t *plp1 = plp[i] + j;
                int qpos = plp1->qpos;
                int bi = bam_seqi(bam_get_seq(plp1->b), qpos);
                char bc = seq_nt16_str[bi];
                printf("%s\t%"PRIhts_pos"\t%d\t%c\t%d\t%d\t%s\n",
                    chr,pos,i,bc, bam_get_qual(plp1->b)[qpos],
                    plp1->b->core.qual, bam_get_qname(plp1->b));
            }
        }
    }
}

int run(int argc, char **argv)
{
    args_t *args = (args_t*) calloc(1,sizeof(args_t));
    static struct option loptions[] =
    {
        {"fasta-ref",required_argument,NULL,'f'},
        {"aln",required_argument,NULL,'a'},
        {"sites",required_argument,NULL,'s'},
        {NULL,0,NULL,0}
    };
    int c;
    while ((c = getopt_long(argc, argv, "o:s:a:f:",loptions,NULL)) >= 0)
    {
        switch (c)
        {
            case 'f': args->fasta_fname = optarg; break;
            case 's': args->sites_fname = optarg; break;
            case 'a': args->aln_fname = optarg; break;
            case 'h':
            case '?':
            default: error("%s", usage_text()); break;
        }
    }
    if ( !args->fasta_fname || !args->sites_fname || !args->aln_fname ) error("%s", usage_text());
    init_data(args);
    process_data(args);
    destroy_data(args);

    return 0;
}
