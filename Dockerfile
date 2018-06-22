# Base Image
FROM centos:7

# Metadata
LABEL base.image="biocontainers:latest"
LABEL version="1"
LABEL software="ssnp"
LABEL software.version="20180622"
LABEL description="ssnp analysis"
LABEL tags="Proteomics"

# Maintainer
MAINTAINER lucky <panyunlai@126.com>

RUN yum install epel-release -y && \
    yum clean all

#安装R和wget以及linux基本工具
RUN yum install -y  \
    R \
    wget \
    bzip2 && \
    yum clean all

#安装配置conda环境
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://mirrors.ustc.edu.cn/anaconda/archive/Anaconda2-5.2.0-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda && \
    rm -f ~/anaconda.sh
ENV PATH /opt/conda/bin:$PATH
ADD .condarc /opt/conda/.condarc

#安装R包（randomForest）
RUN wget --quiet http://mirrors.tuna.tsinghua.edu.cn/CRAN/src/contrib/randomForest_4.6-14.tar.gz -O ~/randomForest_4.6-14.tar.gz && \
    R CMD INSTALL ~/randomForest_4.6-14.tar.gz  && \
    rm -f ~/randomForest_4.6-14.tar.gz

RUN mkdir /data

WORKDIR /data/

CMD ["/bin/bash"]