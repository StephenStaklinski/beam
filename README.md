# Bayesian Evolutionary Analysis of Metastasis (beam)
![beam Logo](logo.jpg)

beam is a BEAST2 package for Bayesian cancer migration graph inference from CRISPR cell lineage tracing data. beam provides a joint model of lineage reconstruction and migration history inference from the raw data, avoiding the need to condition on a single phylogeny while instead inferring a distribution of migration graphs.


## Installation

1. Install BEAST2 from [beast2.org](https://www.beast2.org/)
2. Install BEAGLE from [beagle-dev.github.io](https://beagle-dev.github.io/)
3. Install required BEAST2 packages based on xml declaration. The examples use the following:
   ```bash
   packagemanager -add NS
   packagemanager -add BDSKY
   packagemanager -add BEAST_CLASSIC
   packagemanager -add feast
   ```
4. Install tidetree:
   - Download the latest release from [tidetree GitHub](https://github.com/seidels/tidetree/releases)
   - Find your BEAST2 packages directory (check with `packagemanager -list`)
   - Create and set up tidetree in that directory:
     ```bash
     mkdir tidetree
     cd tidetree
     unzip <tidetree-release-file>.zip
     ```
5. Install beam:
   - Create a `beam` directory in your BEAST2 packages directory
   - Extract the beam release zip file contents into this directory:
   ```bash
     mkdir beam
     cd beam
     unzip <beam-release-file>.zip
     ```
   ```

Verify installation by running `packagemanager -list`. Both beam and tidetree should appear in the list.

## Usage

### Basic Usage

Run beam using the BEAST2 command line interface:

```bash
beast \
  -seed 12345 \
  -threads 5 \
  -overwrite \
  -working \
  examples/quinn_cp26/1.xml
```

