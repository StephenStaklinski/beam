#!/bin/bash

java -Djava.library.path=$BEAGLE_LIB_PATH -Xmx10g -jar /Users/staklins/projects/crispr-barcode-cancer-metastasis/beam_dev/beam/beam.jar \
-seed 12345 \
-threads 5 \
-overwrite \
-working \
/Users/staklins/projects/crispr-barcode-cancer-metastasis/beam_dev/beam/examples/mmus1875_cp01/beam_random_mcmc_template_1875_cp1_debug.xml


