# gencore
A tool to generate consensus reads

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
# step 2: get source (you can also use browser to download from master or releases)
git clone https://github.com/OpenGene/gencore.git

# step 3: build
cd mutscan
make

# step 4: install
sudo make install
```

# Usage
```shell
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
