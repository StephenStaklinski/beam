# Metastabayes BEAST2 package for Bayesian cancer metastasis graph inference from CRISPR cell lineage tracing data

## Inference

Install BEAGLE and set the path to it in `BEAGLE_LIB_PATH`.

```
REPO_DIR="/Users/staklins/projects/crispr-barcode-cancer-metastasis/metastabayes_dev/metastabayes"

java -Djava.library.path=$BEAGLE_LIB_PATH -jar metastabayes.jar \
-seed 1724005224593 \
-D repo=$REPO_DIR \
-overwrite \
-working examples/pR_298.xml
```
