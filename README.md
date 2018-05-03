# gencore
A tool to generate consensus reads from paired-end data. This tool accepts input of a sorted BAM/SAM with its corresponding reference fasta, and outputs an unsorted BAM/SAM.

# Example
```shell
gencore -i input.sorted.bam -o output.bam -r hg19.fasta
```

# Get gencore
## download binary 
This binary is only for Linux systems: http://opengene.org/gencore/gencore
```shell
# this binary was compiled on CentOS, and tested on CentOS/Ubuntu
wget http://opengene.org/gencore/gencore
chmod a+x ./gencore
```
## or compile from source
```shell
# step 1: download and compile htslib from: https://github.com/samtools/htslib
# step 2: get gencore source (you can also use browser to download from master or releases)
git clone https://github.com/OpenGene/gencore.git

# step 3: build
cd gencore
make

# step 4: install
sudo make install
```

# Unique molecular indentifier (UMI) format
`gencore` supports calling consensus reads with or without UMI. Although UMI is not required, it is strongly recommended.   

The UMI should in the tail of query names. It can have a prefix like `UMI`, followed by an underscore. If the UMI has a prefix, it should be specified by `--umi_prefix` or `-u`. It can also have two parts, which are connected by an underscore.   

## Some valid UMI examples:
* Read query name = `NB551106:8:H5Y57BGX2:1:13304:3538:1404:UMI_GAGCATAC`, prefix = `UMI`, umi = `GAGCATAC`
* Read query name = `NB551106:8:H5Y57BGX2:1:13304:3538:1404:UMI_GAGC_ATAC`, prefix = `UMI`, umi = `GAGC_ATAC`
* Read query name = `NB551106:8:H5Y57BGX2:1:13304:3538:1404:GAGCATAC`, prefix = ``, umi = `GAGCATAC`
* Read query name = `NB551106:8:H5Y57BGX2:1:13304:3538:1404:GAGC_ATAC`, prefix = ``, umi = `GAGC_ATAC`

# Usage
```
options:
  -i, --in                   input sorted bam/sam file. STDIN will be read from if it's not specified (string [=-])
  -o, --out                  output bam/sam file. STDOUT will be written to if it's not specified (string [=-])
  -r, --ref                  reference fasta file name (should be an uncompressed .fa/.fasta file) (string)
  -u, --umi_prefix           the prefix for UMI, if it has. None by default. Check the README for the defails of UMI formats. (string [=])
  -s, --supporting_reads     only output consensus reads that merged by >= <supporting_reads> reads. Default value is 2. (int [=2])
      --quit_after_contig    stop when <quit_after_contig> contigs are processed. Only used for fast debugging. Default 0 means no limitation. (int [=0])
      --debug                output some debug information to STDERR.
  -?, --help                 print this message
```
