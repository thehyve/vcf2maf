#!/bin/sh

# Check if VEP_CACHE is set
if [ -z ${VEP_CACHE+x} ]; then 
  echo "VEP_CACHE is unset."
  echo "Please set VEP_CACHE with 'export VEP_CACHE=/<local_path>/<vep_cache_folder>/' in '~/.bashrc' or '~/.bash_profile' and load this with 'source ~/.bash_profile' or 'source ~/.bashrc'."
else 
  echo "VEP_CACHE is set to '$VEP_CACHE'"
  
  # Check if VEP_CACHE exists
  if [ ! -d "$VEP_CACHE" ]; then
    echo "'$VEP_CACHE' does not exist yet, please create it.";
  else

    # Start of VEP cache preparation
    cd $VEP_CACHE

    # Download the Ensembl release and reference genome
    wget -N ftp://ftp.ensembl.org/pub/release-89/variation/VEP/homo_sapiens_vep_89_GRCh37.tar.gz
    wget -N ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz

    # Untar and unzip the Ensembl release and reference genome
    tar -xzf homo_sapiens_vep_89_GRCh37.tar.gz
    gzip -d Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz

    # Download the ExAC VCF
    wget -N ftp://ftp.broadinstitute.org:/pub/ExAC_release/release0.3.1/subsets/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz

    # In a Docker container:
    # - Modify the ExAC VCF
    # - Rename VCF and remove temporary files
    # - Index with tabix
    docker run --rm -it --name vcf2maf \
      -v $VEP_CACHE:/vep_cache/ \
      -w /vep_cache/ \
      thehyve/vcf2maf \
      /bin/bash -c '
        echo "Creating header_line.tmp";
        echo "##FILTER=<ID=AC_Adj0_Filter,Description=\"Only low quality genotype calls containing alternate alleles are present\">" > header_line.tmp;
        echo "Downloading known_somatic_sites.bed";
        curl -sLO https://raw.githubusercontent.com/mskcc/vcf2maf/v1.6.14/data/known_somatic_sites.bed;
        echo "Running bcftools";
        bcftools annotate --header-lines header_line.tmp --remove FMT,^INF/AF,INF/AC,INF/AN,INF/AC_Adj,INF/AN_Adj,INF/AC_AFR,INF/AC_AMR,INF/AC_EAS,INF/AC_FIN,INF/AC_NFE,INF/AC_OTH,INF/AC_SAS,INF/AN_AFR,INF/AN_AMR,INF/AN_EAS,INF/AN_FIN,INF/AN_NFE,INF/AN_OTH,INF/AN_SAS /vep_cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz | bcftools filter --targets-file ^known_somatic_sites.bed --output-type z --output ExAC_nonTCGA.r0.3.1.sites.fixed.vcf.gz;
        echo "Cleaning up";
        mv -f /vep_cache/ExAC_nonTCGA.r0.3.1.sites.fixed.vcf.gz /vep_cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz;
        rm header_line.tmp;
        rm known_somatic_sites.bed;
        echo "Runing tabix";
        tabix -p vcf /vep_cache/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz;
        '

    # In a Docker container use tabix to improve peformance of Ensembl release
    docker run --rm -it --name vcf2maf \
      -v $VEP_CACHE:/vep_cache/ \
      -w /opt/variant_effect_predictor_89/ensembl-tools-release-89/scripts/variant_effect_predictor  \
      thehyve/vcf2maf \
      perl convert_cache.pl \
        --species homo_sapiens \
        --version 89_GRCh37 \
        --dir /vep_cache/

    # In a Docker container, move the Plugins directory to $VEP_CACHE/Plugins
    docker run --rm -it --name vcf2maf \
      -v $VEP_CACHE:/vep_cache/ \
      -w /opt/variant_effect_predictor_89/tmp_plugins-cache \
      thehyve/vcf2maf \
      /bin/bash -c 'mv Plugins/ /vep_cache/'
  fi
fi
