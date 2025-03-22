#!/bin/bash

java -Djava.library.path=$BEAGLE_LIB_PATH -Xmx10g -jar /Users/staklins/projects/crispr-barcode-cancer-metastasis/beam/beam.jar \
-seed 12345 \
-threads 5 \
-overwrite \
-working \
/Users/staklins/projects/crispr-barcode-cancer-metastasis/beam/examples/quinn_cp26/1.xml


