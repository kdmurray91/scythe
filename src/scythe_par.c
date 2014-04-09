/*
  Scythe - A Simple Bayesian Adapter Contamination Trimmer
*/

#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <zlib.h>
#include <omp.h>

#include "scythe.h"
#include "kseq.h"

__KSEQ_BASIC(static, gzFile)
__KS_GETC(gzread, BUFFER_SIZE)
__KS_GETUNTIL(gzread, BUFFER_SIZE)
__KSEQ_READ(static)

static const float default_prior = 0.3;

#ifndef PROGRAM_NAME
#define PROGRAM_NAME "scythe"
#endif

#ifndef AUTHORS
#define AUTHORS "Vince Buffalo, UC Davis\nEmail: <vsbuffaloAAAAA@ucdavis.edu> (poly-A tail removed)"
#endif

#ifndef VERSION
#define VERSION 0.992
#endif

/* Options drawn from GNU's coreutils/src/system.h */
/* These options are defined so as to avoid conflicting with option
   values used by commands */
enum {
  GETOPT_HELP_CHAR = (CHAR_MIN - 2),
  GETOPT_VERSION_CHAR = (CHAR_MIN - 3)
};
#define GETOPT_HELP_OPTION_DECL \
  "help", no_argument, NULL, GETOPT_HELP_CHAR
#define GETOPT_VERSION_OPTION_DECL \
  "version", no_argument, NULL, GETOPT_VERSION_CHAR
#define case_GETOPT_HELP_CHAR                   \
  case GETOPT_HELP_CHAR:                        \
    usage(EXIT_SUCCESS);                        \
    break;
#define case_GETOPT_VERSION_CHAR(Program_name, Version, Authors)        \
  case GETOPT_VERSION_CHAR:                                             \
  fprintf(stdout, "%s version %0.3f\nCopyright (c) 2011 The Regents "   \
          "of University of California, Davis Campus.\n"                \
          "%s is free software and comes with ABSOLUTELY NO WARRANTY.\n" \
          "Distributed under the MIT License.\n\nWritten by %s\n",      \
          Program_name, (double) Version, Program_name, Authors);       \
  exit(EXIT_SUCCESS);                                                   \
  break;
/* end code drawn from system.h */

static void cpy_kstr(kstring_t *dst, const kstring_t *src)
{
	if (src->l == 0) return;
	if (src->l + 1 > dst->m) {
		dst->m = src->l + 1;
		kroundup32(dst->m);
		dst->s = realloc(dst->s, dst->m);
	}
	dst->l = src->l;
	memcpy(dst->s, src->s, src->l + 1);
}

static void cpy_kseq(kseq_t *dst, const kseq_t *src)
{
	cpy_kstr(&dst->name, &src->name);
	cpy_kstr(&dst->seq,  &src->seq);
	cpy_kstr(&dst->qual, &src->qual);
	cpy_kstr(&dst->comment, &src->comment);
}

static struct option long_options[] = {
  /* Options with an argument */
  {"adapter-file", required_argument, 0, 'a'},
  {"prior", required_argument, 0, 'p'},
  {"partial", required_argument, 0, 'P'},
  {"quality-type", required_argument, 0, 'q'},
  {"matches-file", required_argument, 0, 'm'},
  {"output-file", required_argument, 0, 'c'},
  {"min-match", required_argument, 0, 'n'},
  {"min-keep", required_argument, 0, 'M'},
  {"quiet", no_argument, 0, 'Q'},
  {"tag", no_argument, 0, 't'},
  {GETOPT_HELP_OPTION_DECL},
  {GETOPT_VERSION_OPTION_DECL},
  {NULL, 0, NULL, 0}
};

void usage(int status) {
  fputs("\nUsage: scythe -a adapter_file.fasta sequence_file.fastq\n\
Trim 3'-end adapter contaminants off sequence files. If no output file\n\
is specified, scythe will use stdout.\n\
\n\
Options:\n", stdout);
  printf("\
  -p, --prior		prior (default: %0.3f)\n", default_prior);
  fputs("\
  -q, --quality-type	quality type, either illumina, solexa, or sanger (default: sanger)\n\
  -m, --matches-file	matches file (default: no output)\n\
  -o, --output-file	output trimmed sequences file (default: stdout)\n\
  -t, --tag		add a tag to the header indicating Scythe cut a sequence (default: off)\n", stdout);
  fputs("\
  -n, --min-match	smallest contaminant to consider (default: 5)\n\
  -M, --min-keep	filter sequnces less than or equal to this length (default: 35)\n\
  --quiet		don't output statistics about trimming to stdout (default: off)\n\
  --help		display this help and exit\n\
  --version		output version information and exit\n", stdout);
  fputs("\n\
  Information on quality schemes:\n\
  phred			PHRED quality scores (e.g. from Roche 454). ASCII with no offset, range: [4, 60].\n\
  sanger		Sanger are PHRED ASCII qualities with an offset of 33, range: [0, 93]. From \n\
			NCBI SRA, or Illumina pipeline 1.8+.\n\
  solexa		Solexa (also very early Illumina - pipeline < 1.3). ASCII offset of\n", stdout);
  fputs("\
	 		64, range: [-5, 62]. Uses a different quality-to-probabilities conversion than other\n\
			schemes.\n\
  illumina		Illumina output from pipeline versions between 1.3 and 1.7. ASCII offset of 64,\n\
			range: [0, 62]\n", stdout);
  exit(status);
}

int main(int argc, char *argv[]) {
  int min=5, min_keep=35, index;
  int debug=0, verbose=1;
  int contaminated=0, total=0;
  quality_type qual_type=SANGER;
  float prior=default_prior;
  adapter_array *aa;
  gzFile adapter_fp=NULL, fp;
  FILE *output_fp=stdout, *matches_fp=NULL;
  int optc;
  int add_tag = 0;
  extern char *optarg;
  int need_usage = 0;

  while (1) {
    int option_index = 0;
    optc = getopt_long(argc, argv, "dtfp:a:o:q:m:o:n:M:", long_options, &option_index);

    if (optc == -1)
       break;
    switch (optc) {
      if (long_options[option_index].flag != 0)
        break;
      case 'a':
        adapter_fp = gzopen(optarg, "r");
        if (!adapter_fp) {
          fprintf(stderr, "Could not open adapter file '%s'.\n", optarg);
          return EXIT_FAILURE;
        }
        break;
      case 'd':
        debug = 1;
        break;
      case 't':
        add_tag = 1;
        break;
      case 'Q':
        verbose = 0;
        break;
      case 'o':
        output_fp = fopen(optarg, "w");
        if (!output_fp) {
          fprintf(stderr, "Could not open output file '%s'.\n", optarg);
          return EXIT_FAILURE;
        }
        break;
      case 'm':
        matches_fp = fopen(optarg, "w");
        if (!matches_fp) {
          fprintf(stderr, "Could not open matches file '%s'.\n", optarg);
          return EXIT_FAILURE;
        }
        break;
      case 'n':
        min = atoi(optarg);
        break;
      case 'M':
        min_keep = atoi(optarg);
        break;
      case 'q':
        if (strcmp(optarg, "illumina") == 0)
          qual_type = ILLUMINA;
        else if (strcmp(optarg, "solexa") == 0)
          qual_type = SOLEXA;
        else if (strcmp(optarg, "sanger") == 0)
          qual_type = SANGER;
        else {
          fprintf(stderr, "Unknown quality type '%s'.\n", optarg);
          usage(EXIT_FAILURE);
        }          
        break;
      case 'p':
        /* This truncation is acceptable... priors need not be doubles */
        prior = (float) atof(optarg);
        if (prior > 1 || prior < 0) {
          fprintf(stderr, "Prior must be between 0 and 1\n");
          usage(EXIT_FAILURE);
        }
        break;
        
      case_GETOPT_HELP_CHAR;
      case_GETOPT_VERSION_CHAR(PROGRAM_NAME, VERSION, AUTHORS);
      case '?':
        break;
      default:
        usage(EXIT_FAILURE);
      }
  }
  
  if (debug) {
    matches_fp = stdout;
    output_fp = stdout;
  }
  
  if ((index = optind) == argc) {
    fprintf(stderr, "No FASTQ file specified.\n");
    usage(EXIT_FAILURE);
  }
 
  /* load all adapter sequences into memory */
  if (!adapter_fp) {
    fprintf(stderr, "No adapter file specified.\n");
    usage(EXIT_FAILURE);
  }
  aa = load_adapters(adapter_fp);
  gzclose(adapter_fp);

  fp = strcmp(argv[index], "-") ? gzopen(argv[index], "r") : gzdopen(fileno(stdin), "r");

  if (!fp) {
    fprintf(stderr, "FASTQ file '%s' not found.\n", argv[index]);
    return EXIT_FAILURE;
  }


  /* Loop through entire sequence file. Write trimmed sequences to
     file (or stdout), and record matches in a match file if specifed.
  */
  kseq_t *seq;
  seq = kseq_init(fp);
  #pragma omp parallel default(none) shared(fp, matches_fp, output_fp, stderr, min_keep, need_usage, index, qual_type, prior, min, aa, contaminated, add_tag, total, seq)
  {
    int l, shift;
    match *best_match;
    float *qprobs;
    int our_cont=0, our_total=0;
    kseq_t *seq_cpy = calloc(1, sizeof(*seq_cpy));
    do {
      #pragma omp critical
      {
        l = kseq_read(seq);
        cpy_kseq(seq_cpy, seq);
      }
      if (l < 0) {
        break;
      }
      shift = -1;
      if (!seq_cpy->qual.s) {
        #pragma omp critical
        {
          fputs("Sequence file missing or has malformed quality line.\n", stderr);
          need_usage = 1;
        }
        if (need_usage) break;
      }

      qprobs = qual_to_probs(seq_cpy->qual.s, qual_type);
      best_match = find_best_match(aa, seq_cpy->seq.s, qprobs, prior, 0.25, min);

      if (best_match && best_match->ps->is_contam) {
        our_cont++;
        shift = best_match->shift;
        #pragma omp critical
        {
          if (matches_fp) print_match(seq_cpy, best_match, matches_fp, aa, qual_type);
        }
        /* TODO */
        /* aa->adapters[best_match->adapter_index].occurrences[best_match->n-1]++; */
      }

      #pragma omp critical
      {
        write_fastq(output_fp, seq_cpy, add_tag, shift, min_keep);
        total++;
      }
      if (best_match) destroy_match(best_match);
      free(qprobs);
    } while (1);
    #pragma omp critical
    {
      total += our_total;
      contaminated += our_cont;
    }

    kseq_destroy(seq_cpy);
  }
  kseq_destroy(seq);
  
  if (need_usage) {
    usage(EXIT_FAILURE);
  }
  if (verbose) 
    print_summary(aa, prior, total-contaminated, contaminated, total);

  destroy_adapters(aa, MAX_ADAPTERS);
  gzclose(fp);
  return 0;
}
