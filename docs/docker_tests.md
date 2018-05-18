# vcf2maf tests
- [Test VEP](#test-vep)
- [Test vcf2maf](#test-vcf2maf)
- [Test maf2maf](#test-maf2maf)
- [Start and browse the Docker container](#start-and-browse-the-docker-container)

## Test VEP
Executing the following Docker command will test VEP.
```bash
docker run --rm -it --name vcf2maf \
  -v $VEP_CACHE:/vep_cache/ \
  -w /opt/variant_effect_predictor_89/ensembl-tools-release-89/scripts/variant_effect_predictor \
  thehyve/vcf2maf \
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
    --output_file example_GRCh37_vep.vcf \
    --custom /vep_cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz,ExAC,vcf,exact,1,AC,AN \
    --cache_version 89
```

## Test vcf2maf
Executing the following Docker command will test vcf2maf.
```bash
docker run --rm -it --name vcf2maf \
  -v $VEP_CACHE:/vep_cache/ \
  -w /opt/vcf2maf \
  thehyve/vcf2maf \
  perl vcf2maf.pl \
    --input-vcf tests/test.vcf \
    --output-maf tests/test_vcf2maf.maf \
    --vep-path /opt/variant_effect_predictor_89/ensembl-tools-release-89/scripts/variant_effect_predictor \
    --vep-data /vep_cache/ \
    --ref-fasta /vep_cache/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    --filter-vcf /vep_cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
    --custom-enst data/isoform_overrides_uniprot \
    --cache-version 89
```

## Test maf2maf
Executing the following Docker command will test maf2maf.
```bash
docker run --rm -it --name vcf2maf \
  -v $VEP_CACHE:/vep_cache/ \
  -w /opt/vcf2maf \
  thehyve/vcf2maf \
  perl maf2maf.pl \
    --input-maf tests/test.maf \
    --output-maf tests/test_maf2maf.maf \
    --vep-path /opt/variant_effect_predictor_89/ensembl-tools-release-89/scripts/variant_effect_predictor \
    --vep-data /vep_cache/ \
    --ref-fasta /vep_cache/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
    --filter-vcf /vep_cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz \
    --custom-enst data/isoform_overrides_uniprot \
    --cache-version 89
```

## Start and browse the Docker container
To start and browse the Docker container, use the following run command:
```bash
docker run --rm --name vcf2maf -it \
  -v $VEP_CACHE/vep_cache/:/vep_cache/ \
  -w /opt/ \
  thehyve/vcf2maf \
  bash
```
