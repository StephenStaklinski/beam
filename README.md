# Bayesian Evolutionary Analysis of Metastasis (BEAM)
![BEAM Logo](logo.jpg)

BEAM is a BEAST2 package for Bayesian cancer migration graph inference from CRISPR cell lineage tracing data. BEAM provides a joint model of lineage reconstruction and migration history inference from the raw data, avoiding the need to condition on a single phylogeny while instead inferring a distribution of migration graphs.


## Inference

Install java+fx 17.0.9 and BEAGLE and set the path to it in `BEAGLE_LIB_PATH`. Then, the provided `beam.jar` file is all that is needed to run the method after installing dependent packages:

```
java -cp beam.jar beast.pkgmgmt.PackageManager -add NS
java -cp beam.jar beast.pkgmgmt.PackageManager -add BEAST_CLASSIC
java -cp beam.jar beast.pkgmgmt.PackageManager -add BEASTLabs
java -cp beam.jar beast.pkgmgmt.PackageManager -add feast

# specify path to the repo
REPO_DIR="/Users/staklins/projects/crispr-barcode-cancer-metastasis/beam_dev/beam"

# submit example data to run
java -Djava.library.path=$BEAGLE_LIB_PATH -jar beam.jar \
-seed 1724005224593 \
-D repo=$REPO_DIR \
-overwrite \
-working examples/pR_298.xml
```

## Development

An initial install will obtain dependencies for developing on the code base by running the following:
```
ant
```



