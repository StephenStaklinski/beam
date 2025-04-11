#!/bin/bash

java -Djava.library.path=$BEAGLE_LIB_PATH -Xmx10g -jar beam.jar \
-seed 12345 \
-threads 5 \
-overwrite \
-working \
examples/ns_random_mmus1495_cp112/45.xml
