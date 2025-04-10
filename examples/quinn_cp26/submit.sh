#!/bin/bash

java -Djava.library.path=$BEAGLE_LIB_PATH -Xmx10g -jar /grid/siepel/home/staklins/beam/beam.jar \
-seed 12345 \
-threads 5 \
-overwrite \
-working \
/grid/siepel/home/staklins/beam/examples/quinn_cp26/1.xml


