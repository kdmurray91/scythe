# Scythe - A very simple adapter trimmer (version 0.93 BETA)

Contact: Vince Buffalo <vsbuffaloAAAAA@gmail.com> (with the poly-A tail removed)

Copyright (c) 2011 The Regents of University of California, Davis Campus.

## About

Scythe uses a Naive Bayesian approach to classify contaminant
substrings in sequence reads. It considers quality information, which
can make it robust in picking out 3'-end adapters, which often include
poor quality bases.

Most next generation sequencing reads have deteriorating quality
towards the 3'-end. It's common for a quality-based trimmer to be
employed before mapping, assemblies, and analysis to remove these poor
quality bases. However, quality-based trimming could remove bases that
are helpful in identifying (and removing) 3'-end adapter
contaminants. Thus, it is recommended you run Scythe *before*
quality-based trimming, as part of a read quality control pipeline.

The Bayesian approach Scythe uses compares two likelihood models: the
probability of seeing the matches in a sequence given contamination,
and not given contamination. Given that the read is contaminated, the
probability of seeing a certain number of matches and mistmatches is a
function of the quality of the sequence. Given the read is not
contaminated (and is thus assumed to be random sequence), the
probability of seeing a certain number of matches and mismatches is
chance. The posterior is calculated across both these likelihood
models, and the class (contaminated or not contaminated) with the
maximum posterior probability is the class selected.

## Requirements

Scythe can be compiled using GCC or Clang; compilation during
development used the latter. Scythe relies on Heng Li's kseq.h, which
is bundled with the source.

Scythe requires Zlib, which can be obtained at <http://www.zlib.net/>.

## Building and Installing Scythe

To build Scythe, cd into the code/ directory and enter:

    make build

Then, copy or move "scythe" to a directory in your $PATH.

## Usage

Scythe can be run minimally with:

    scythe -a adapter_file.fasta -o trimmed_sequences.fasta sequences.fastq

By default, the prior contamination rate is 0.05. This can be changed
(and one is encouraged to do so!) with:

    scythe -a adapter_file.fasta -p 0.1 -o trimmed_sequences.fastq sequences.fastq

If you'd like to use standard out, it is recommended you use the
--quiet option:

    scythe -a adapter_file.fasta --quiet sequences.fastq > trimmed_sequences.fastq

Also, more detailed output about matches can be obtained with:

    scythe -a adapter_file.fasta -o trimmed_sequences.fasta -m matches.txt sequences.fastq

By default, Illumina's quality scheme (pipeline > 1.3) is used. Sanger
or Solexa (pipeline < 1.3) qualities can be specified with -q:

    scythe -a adapter_file.fasta -q solexa -o trimmed_sequences.fasta sequences.fastq

Lastly, a minimum match length argument can be specified with -n <integer>:

    scythe -a adapter_file.fasta -n 4 -o trimmed_sequences.fasta sequences.fastq

## Notes

Scythe only checks for 3'-end contaminants, up to the adapter's length
into the 3'-end. For reads with contamination in *any* position, the
program TagDust (<http://genome.gsc.riken.jp/osc/english/dataresource/>)
is recommended. Scythe has the advantages of allowing fuzzier matching
and being base quality-aware, while TagDust has the advantages of very
fast matching (but allowing few mismatches, and not considering
quality) and FDR. TagDust also removes contaminated reads *entirely*, while
Scythe trims off contaminants. 

A possible pipeline would run FASTQ reads through Scythe, then
TagDust, then a quality-based trimmer, and finally through a read
quality statistics program such as qrqc
(<https://github.com/vsbuffalo/qrqc>) or FASTqc
(<http://www.bioinformatics.bbsrc.ac.uk/projects/fastqc/>).
