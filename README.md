# ssnp

# 使用说明

    nextflow run dggene/ssnp --input=XXXX/xx.vcf -N user@host.com

## 如果需要在pbs上运行，需要在本地目录下添加netxflow.config 配置文件

    process{
        executor='pbs'
        queue='nextflow'
    }

# 环境依赖

需要安装有anaconda环境