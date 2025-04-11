# Bayesian Evolutionary Analysis of Metastasis (BEAM)
![BEAM Logo](logo.jpg)

BEAM is a BEAST2 package for Bayesian cancer migration graph inference from CRISPR cell lineage tracing data. BEAM provides a joint model of lineage reconstruction and migration history inference from the raw data, avoiding the need to condition on a single phylogeny while instead inferring a distribution of migration graphs.


## Installation

1. Install beast2 from [beast2.org](https://www.beast2.org/)
1. Install beagle from [beagle-dev.github.io](https://beagle-dev.github.io/)
1. Use beast2 package manager to install any packages used in xml declaration. The examples depend on NS, BDSKY, BEAST_CLASSIC, and feast.
1. Install tidetree from source code by obtaining the zip file from [tidetree github](https://github.com/seidels/tidetree/releases). Find where beast2 packages are installed, which should be indicated from `packagemanager -list` and then move to that directory, run `mkdir tidetree && cd tidetree`, and unzip the release zip file here.
1. Install beam zip file in the same way into a beam directory in the packages directory.

Now, both beam and tidetree should be listed when running `packagemanager -list`.

## Inference

Several examples are provides in `examples/` but here is one of them:
```
beast \
-seed 12345 \
-threads 5 \
-overwrite \
-working \
examples/quinn_cp26/1.xml
```

