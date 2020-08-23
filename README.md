# AmpBinner
A barcode demultiplexer for Oxford Nanopore long-read amplicon sequencing data. 

## Features
- AmpBinner supports multiple barcoding strategies and custom-designed barcodes
- AmpBinner uses the sequence upstream/downstream of the barcode to help locate barcode position and eliminates random matching
- AmpBinner is able to demultiplex Oxford Nanopore sequencing data generated from **10X Genomics Chromium single cell libraries**. 

## <a name="Requirements"></a>Requirements
- Operating system: Linux or macOS
- [Python](https://www.python.org/) 2.7 or later
- [minimap2](https://github.com/lh3/minimap2) 2.8 or later

## <a name="Installation"></a>Installation

If you don't have `minimap2` in you system, you can install it following the instructions [here](https://github.com/lh3/minimap2#install).  
If you are using Linux, you can acquire precompiled binaries using the following commands:

```
wget https://github.com/lh3/minimap2/releases/download/v2.17/minimap2-2.17_x64-linux.tar.bz2
tar -jxvf minimap2-2.17_x64-linux.tar.bz2
./minimap2-2.17_x64-linux/minimap2
```

Next, you can clone the repository of AmpBinner using the following command.
```
git clone https://github.com/WGLab/AmpBinner.git
```
## <a name="Usage"></a>Usage

### Demultiplexing regular amplicons


### Demultiplexing 10X Genomics Chromium single cell libraries


