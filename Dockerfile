# Base Image
FROM biocontainers/biocontainers:latest

# Metadata
LABEL base.image="biocontainers:latest"
LABEL version="1"
LABEL software="ssnp"
LABEL software.version="20180622"
LABEL description="ssnp analysis"
LABEL tags="Proteomics"

# Maintainer
MAINTAINER lucky <panyunlai@126.com>

USER biodocker

RUN conda install -c r r-randomforest -y

WORKDIR /data/

CMD ["/bin/bash"]