FROM ubuntu:14.04
RUN apt-get update && apt-get install -y \
	autoconf \
	automake \
	make \
	g++ \
	gcc \
	build-essential \ 
	zlib1g-dev \
	libgsl0-dev \
	perl \
	curl \
	git \
	wget \
	unzip \
	tabix \
	libncurses5-dev

RUN apt-get install -y cpanminus
RUN apt-get install -y libmysqlclient-dev
RUN cpanm CPAN::Meta \
	Archive::Zip \
	DBI \
	DBD::mysql \ 
	JSON \
	DBD::SQLite \
	Set::IntervalTree \
	LWP \
	LWP::Simple \
	Archive::Extract \
	Archive::Tar \
	Archive::Zip \
	CGI \
	Time::HiRes \
	Encode \
	File::Copy::Recursive \
	Perl::OSType \
	Module::Metadata version \
	Bio::Root::Version \
	TAP::Harness \
	Module::Build

WORKDIR /opt
RUN curl -LOOO https://github.com/samtools/{samtools/releases/download/1.3.1/samtools-1.3.1,bcftools/releases/download/1.3.1/bcftools-1.3.1}.tar.bz2
RUN cat *tar.bz2 | tar -ijxf -
RUN cd samtools-1.3.1 && make && make install && ln -s /opt/samtools-1.3.1/samtools /usr/bin/samtools && cd ..
RUN cd bcftools-1.3.1 && make && make install && ln -s /opt/bcftools-1.3.1/bcftools /usr/bin/bcftools && cd ..
RUN rm *.tar.bz2

RUN wget https://github.com/Ensembl/ensembl-tools/archive/release/89.zip
RUN mkdir variant_effect_predictor_89
RUN unzip 89.zip -d variant_effect_predictor_89
RUN rm 89.zip
WORKDIR /opt/variant_effect_predictor_89/ensembl-tools-release-89/scripts/variant_effect_predictor/
RUN perl INSTALL.pl --AUTO ap --PLUGINS LoF --CACHEDIR /opt/variant_effect_predictor_89/tmp_plugins-cache

# During installation of VEP, we don't install the full cache, only the plugins.
# When preparing the VEP cache in a separate step (docker_vep_cache.sh) the 
# Plugins folder is moved to the local VEP cache directorty.

WORKDIR /opt/variant_effect_predictor_89/tmp_plugins-cache/Plugins
RUN wget https://raw.githubusercontent.com/konradjk/loftee/v0.3-beta/splice_module.pl

WORKDIR /opt
COPY . /opt/vcf2maf

LABEL maintainers=" \
 Michele Mattioni (Seven Bridges) <michele.mattioni@sbgenomics.com>, \
 Dionne Zaal (The Hyve) <dionne@thehyve.nl>, \
 Sander Tan (The Hyve) <sandertan@thehyve.nl> \
 "
