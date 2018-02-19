vcf<img src="http://i.giphy.com/R6X7GehJWQYms.gif" width="30">maf
=======

To convert a [VCF](http://samtools.github.io/hts-specs/) into a [MAF](https://wiki.nci.nih.gov/x/eJaPAQ), each variant must be mapped to only one of all possible gene transcripts/isoforms that it might affect. But even within a single isoform, a `Missense_Mutation` close enough to a `Splice_Site`, can be labeled as either in MAF format, but not as both. **This selection of a single effect per variant, is often subjective. And that's what this project attempts to standardize.** The `vcf2maf` and `maf2maf` scripts leave most of that responsibility to [Ensembl's VEP](http://useast.ensembl.org/info/docs/tools/vep/index.html), but allows you to override their "canonical" isoforms, or use a custom ExAC VCF for annotation. Though the most useful feature is the **extensive support in parsing a wide range of crappy MAF-like or VCF-like formats** we've seen out in the wild.

# Table of Contents
- [Quick start (using Docker)](#quick-start-using-docker)
- [Quick start (manual installation)](#quick-start-manual-installation)
- [maf2maf](#maf2maf)
- [License](#license)

# Quick start (using Docker)

This documentation is specific for:
- Ensembl release 89
- Ensembl VEP 89 as part of `ensembl-tools`
- GRCh37 / hg19

Since release Ensembl release 90, VEP has moved from `ensembl-tools` to `ensembl-vep` and has significantly changed functionality. All [`ensembl-vep`](https://github.com/Ensembl/ensembl-vep) versions (even 88 and 89 from [DockerHub](https://hub.docker.com/r/ensemblorg/ensembl-vep/)) are not compatible with `vcf2maf`. This Docker approach uses the latest release from [`ensembl-tools`](https://github.com/Ensembl/ensembl-tools/tree/release/89/scripts/variant_effect_predictor).

### Create Docker image
```bash
git clone --branch docker_improvements https://github.com/thehyve/vcf2maf
cd vcf2maf
docker build -t vcf2maf .
```

### Create local cache directory
First create a directory for the VEP cache folder.
```bash
mkdir /<local_path>/vep_cache
```

Add this path as environment variable to `~/.bash_profile` or `~/.bashrc`.
```bash
export VEP_CACHE=/<local_path>/vep_cache/
```
Load this variable with `source ~/.bash_profile` or `source ~/.bashrc`.

### Download Ensembl release and reference genome
```bash
cd $VEP_CACHE
wget ftp://ftp.ensembl.org/pub/release-89/variation/VEP/homo_sapiens_vep_89_GRCh37.tar.gz
wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz

tar -xzf homo_sapiens_vep_89_GRCh37.tar.gz
gzip -d Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
```

### Download and prepare ExAC VCF
Download and modify the ExAC r0.3.1 VCF with germline variants called across thousands of normal samples excluding TCGA as described in as in https://gist.github.com/ckandoth/f265ea7c59a880e28b1e533a6e935697.

```bash
# Download the ExAC VCF
cd $VEP_CACHE
wget ftp://ftp.broadinstitute.org:/pub/ExAC_release/release0.3.1/subsets/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz

# Start Docker container
docker run --name vcf2maf -it -v $VEP_CACHE:/vep_cache/ -w /vep_cache/ --rm vcf2maf bash

# In the Docker container modify the ExAC VCF
echo "##FILTER=<ID=AC_Adj0_Filter,Description=\"Only low quality genotype calls containing alternate alleles are present\">" > header_line.tmp
curl -LO https://raw.githubusercontent.com/mskcc/vcf2maf/v1.6.14/data/known_somatic_sites.bed
bcftools annotate --header-lines header_line.tmp --remove FMT,^INF/AF,INF/AC,INF/AN,INF/AC_Adj,INF/AN_Adj,INF/AC_AFR,INF/AC_AMR,INF/AC_EAS,INF/AC_FIN,INF/AC_NFE,INF/AC_OTH,INF/AC_SAS,INF/AN_AFR,INF/AN_AMR,INF/AN_EAS,INF/AN_FIN,INF/AN_NFE,INF/AN_OTH,INF/AN_SAS /vep_cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz | bcftools filter --targets-file ^known_somatic_sites.bed --output-type z --output ExAC_nonTCGA.r0.3.1.sites.fixed.vcf.gz

# Rename VCF and remove temporary files
mv -f /vep_cache/ExAC_nonTCGA.r0.3.1.sites.fixed.vcf.gz /vep_cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
rm header_line.tmp
rm known_somatic_sites.bed

# Index with tabix
tabix -p vcf /vep_cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
```

### Prepare Ensembl release with tabix
Convert the offline cache for use with tabix, which significantly speeds up the lookup of known variants as described in https://gist.github.com/ckandoth/f265ea7c59a880e28b1e533a6e935697.
```bash
docker run --rm -it --name vcf2maf \
  -v $VEP_CACHE:/vep_cache/ \
  -w /opt/variant_effect_predictor_89/ensembl-tools-release-89/scripts/variant_effect_predictor  \
  vcf2maf \
  perl convert_cache.pl \
    --species homo_sapiens \
    --version 89_GRCh37 \
    --dir /vep_cache/
```

### Test VEP
```bash
# Adding --cache_version 89 is required, else VEP searches for 88
docker run --rm -it --name vcf2maf \
  -v $VEP_CACHE:/vep_cache/ \
  -w /opt/variant_effect_predictor_89/ensembl-tools-release-89/scripts/variant_effect_predictor  \
  vcf2maf \
  perl variant_effect_predictor.pl --species homo_sapiens --assembly GRCh37 --offline \
    --no_progress --no_stats --sift b --ccds --uniprot --hgvs --symbol --numbers \
    --domains --gene_phenotype --canonical --protein --biotype --tsl --pubmed \
    --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number \
    --no_escape --xref_refseq --failed 1 --vcf --minimal --flag_pick_allele \
    --polyphen b --gmaf --maf_1kg --maf_es --regulatory \
    --pick_order canonical,tsl,biotype,rank,ccds,length \
    --dir /vep_cache/ \
    --fasta /vep_cache/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    --input_file example_GRCh37.vcf \
    --output_file TEST_output_example_GRCh37.vep.vcf \
    --custom /vep_cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz,ExAC,vcf,exact,1,AC,AN \
    --cache_version 89
```

### Test vcf2maf
``` bash
# Adding --cache_version 89 is required, else VEP searches for 88
docker run --rm -it --name vcf2maf \
  -v $VEP_CACHE:/vep_cache/ \
  -w /opt/vcf2maf  \
  vcf2maf \
  perl vcf2maf.pl \
    --input-vcf tests/test.vcf \
    --output-maf tests/TEST_OUTPUT.vep.maf \
    --vep-path /opt/variant_effect_predictor_89/ensembl-tools-release-89/scripts/variant_effect_predictor \
    --vep-data /vep_cache/ \
    --ref-fasta /vep_cache/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    --filter-vcf /vep_cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
    --custom-enst data/isoform_overrides_uniprot \
    --cache-version 89
```

### Run vcf2maf
To run `vcf2maf` with a local VCF file, use the [test example](#test-vcf2maf) and mount a directory for input and output using the `-v` command in `docker run`. For example: `-v /local_input_output/:/input_output/`. In the `vcf2maf.pl` command, direct to the input and output files, for example: `--input-vcf /input_output/input.vcf` and `--output-maf /input_output/output.maf`.

### Start and browse the Docker container
To start and browse the Docker container, use the following run command:
```bash
docker run --name vcf2maf -it -v $VEP_CACHE/vep_cache/:/vep_cache/ -w /opt/ --rm vcf2maf bash
```

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
