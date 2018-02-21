# vcf2maf

To convert a [VCF](http://samtools.github.io/hts-specs/) into a [MAF](https://wiki.nci.nih.gov/x/eJaPAQ), each variant must be mapped to only one of all possible gene transcripts/isoforms that it might affect. But even within a single isoform, a `Missense_Mutation` close enough to a `Splice_Site`, can be labeled as either in MAF format, but not as both. **This selection of a single effect per variant, is often subjective. And that's what this project attempts to standardize.** The `vcf2maf` and `maf2maf` scripts leave most of that responsibility to [Ensembl's VEP](http://useast.ensembl.org/info/docs/tools/vep/index.html), but allows you to override their "canonical" isoforms, or use a custom ExAC VCF for annotation. Though the most useful feature is the **extensive support in parsing a wide range of crappy MAF-like or VCF-like formats** we've seen out in the wild.

# Table of Contents
- [Quick start (using Docker)](#quick-start-using-docker)
  - [Create local VEP cache directory](#create-local-vep-cache-directory)
- [Quick start (manual installation)](#quick-start-manual-installation)
- [maf2maf](#maf2maf)
- [License](#license)

# Quick start (using Docker)

This documentation is specific for:
- Ensembl release 89
- Ensembl VEP 89 as part of `ensembl-tools`
- GRCh37 / hg19

Since release Ensembl release 90, VEP has moved from `ensembl-tools` to `ensembl-vep` and has significantly changed functionality. All [`ensembl-vep`](https://github.com/Ensembl/ensembl-vep) versions (even 88 and 89 from [DockerHub](https://hub.docker.com/r/ensemblorg/ensembl-vep/)) are not compatible with `vcf2maf`. This Docker approach uses the latest release from [`ensembl-tools`](https://github.com/Ensembl/ensembl-tools/tree/release/89/scripts/variant_effect_predictor).

## Create Docker image

### Option 1: Docker Hub
The easiest way to obtain the Docker image is to pull it from Docker Hub:
```bash
docker pull thehyve/vcf2maf
```

### Option 2: Build it from Git repository
If you would like to build your own Docker image, you could clone the Git repository and build it with Docker:
```bash
git clone --branch master https://github.com/thehyve/vcf2maf
cd vcf2maf
docker build -t vcf2maf .
```
When buiding from a local Git repository, substitute the `thehyve/vcf2maf` image name by `vcf2maf` in the subsequent commands, as well as in `vep_cache_preparation.sh`.

## Create local VEP cache directory
First create a directory for the VEP cache folder.
```bash
mkdir /<local_path>/<vep_cache_folder>
```

Add this path as environment variable to `~/.bash_profile` or `~/.bashrc`.
```bash
export VEP_CACHE=/<local_path>/<vep_cache_folder>/
```
Load this variable with `source ~/.bash_profile` or `source ~/.bashrc`.

Creating the cache directory includes downloading the Ensembl release, reference genome and ExAC VCF. **This will take several hours.**
```bash
/bin/bash vep_cache_preparation.sh
```

## Test VEP, vcf2maf and maf2maf
Tests can be found in the [Tests markdown file](docs/vcf2maf_tests.md).

## Run vcf2maf
To run `vcf2maf` with a local VCF file, use the [test example](docs/vcf2maf_tests.md#test-vcf2maf) and mount a directory for input and output using the `-v` command in `docker run`. For example: `-v /local_input_output/:/input_output/`. In the `vcf2maf.pl` command, direct to the input and output files, for example: `--input-vcf /input_output/input.vcf` and `--output-maf /input_output/output.maf`.

# Quick start (manual installation)

Find the [latest stable release](https://github.com/mskcc/vcf2maf/releases), download it, and view the detailed usage manuals for `vcf2maf` and `maf2maf`:

    export VCF2MAF_URL=`curl -sL https://api.github.com/repos/mskcc/vcf2maf/releases | grep -m1 tarball_url | cut -d\" -f4`
    curl -L -o mskcc-vcf2maf.tar.gz $VCF2MAF_URL; tar -zxf mskcc-vcf2maf.tar.gz; cd mskcc-vcf2maf-*
    perl vcf2maf.pl --man
    perl maf2maf.pl --man

If you don't have [VEP](http://useast.ensembl.org/info/docs/tools/vep/index.html) installed, then [follow this gist](https://gist.github.com/ckandoth/f265ea7c59a880e28b1e533a6e935697). Of the many annotators out there, VEP is preferred for its large team of active coders, and its CLIA-compliant [HGVS formats](http://www.hgvs.org/mutnomen/recs.html). After installing VEP, you can test the script like so:

    perl vcf2maf.pl --input-vcf tests/test.vcf --output-maf tests/test.vep.maf

To fill columns 16 and 17 of the output MAF with tumor/normal sample IDs, and to parse out genotypes and allele counts from matched genotype columns in the VCF, use options `--tumor-id` and `--normal-id`. Skip option `--normal-id` if you didn't have a matched normal:

    perl vcf2maf.pl --input-vcf tests/test.vcf --output-maf tests/test.vep.maf --tumor-id WD1309 --normal-id NB1308

VCFs from variant callers like [VarScan](http://varscan.sourceforge.net/somatic-calling.html#somatic-output) use hardcoded sample IDs TUMOR/NORMAL in the genotype columns of the VCF. To have this script correctly parse the correct genotype columns, while still printing the proper IDs in the output MAF:

    perl vcf2maf.pl --input-vcf tests/test_varscan.vcf --output-maf tests/test_varscan.vep.maf --tumor-id WD1309 --normal-id NB1308 --vcf-tumor-id TUMOR --vcf-normal-id NORMAL

If you have the VEP script in a different folder like `/opt/vep`, and its cache in `/srv/vep`, there are options available to use those instead:

    perl vcf2maf.pl --input-vcf tests/test.vcf --output-maf tests/test.vep.maf --vep-path /opt/vep --vep-data /srv/vep

# maf2maf

If you have a MAF or a MAF-like file that you want to reannotate, then use `maf2maf`, which simply runs `maf2vcf` followed by `vcf2maf`:

    perl maf2maf.pl --input-maf tests/test.maf --output-maf tests/test.vep.maf

After tests on variant lists from many sources, `maf2vcf` and `maf2maf` are quite good at dealing with formatting errors or "MAF-like" files. It even supports VCF-style alleles, as long as `Start_Position == POS`. But it's OK if the input format is imperfect. Any variants with a reference allele mismatch are kept aside in a separate file for debugging. The bare minimum columns that `maf2maf` expects as input are:

    Chromosome	Start_Position	Reference_Allele	Tumor_Seq_Allele2	Tumor_Sample_Barcode
    1	3599659	C	T	TCGA-A1-A0SF-01
    1	6676836	A	AGC	TCGA-A1-A0SF-01
    1	7886690	G	A	TCGA-A1-A0SI-01

See `data/minimalist_test_maf.tsv` for a sampler. Addition of `Tumor_Seq_Allele1` will be used to determine zygosity. Otherwise, it will try to determine zygosity from variant allele fractions, assuming that arguments `--tum-vad-col` and `--tum-depth-col` are set correctly to the names of columns containing those read counts. Specifying the `Matched_Norm_Sample_Barcode` with its respective columns containing read-counts, is also strongly recommended. Columns containing normal allele read counts can be specified using argument `--nrm-vad-col` and `--nrm-depth-col`.

# License

    Apache-2.0 | Apache License, Version 2.0 | https://www.apache.org/licenses/LICENSE-2.0
