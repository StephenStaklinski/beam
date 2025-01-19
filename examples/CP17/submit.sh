REPO_DIR="/Users/staklins/projects/crispr-barcode-cancer-metastasis/beam_dev/beam"

FILE_DIR="${REPO_DIR}/examples/CP17"

java  -Djava.library.path=$BEAGLE_LIB_PATH -Xmx10g -jar $REPO_DIR/beam.jar \
-threads 5 \
-overwrite \
-working \
-D fileDir="$FILE_DIR" \
-seed 1736521990304 $FILE_DIR/3.xml > $FILE_DIR/terminal_3.log


# random model
java  -Djava.library.path=$BEAGLE_LIB_PATH -Xmx10g -jar $REPO_DIR/beam.jar \
-threads 5 \
-overwrite \
-working \
-D fileDir="$FILE_DIR" \
-seed 1736521990304 $FILE_DIR/3_random.xml > $FILE_DIR/terminal_3_random.log