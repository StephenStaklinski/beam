REPO_DIR="/Users/staklins/projects/crispr-barcode-cancer-metastasis/beam_dev/beam"

EXAMPLE_DIR="${REPO_DIR}/examples/CP17"


java  -Djava.library.path=$BEAGLE_LIB_PATH -Xmx10g -jar $REPO_DIR/beam.jar \
-threads 5 \
-overwrite \
-working \
-seed 1736521990304 $EXAMPLE_DIR/3.xml > $EXAMPLE_DIR/terminal_3.log



# REPO_DIR="/Users/staklins/projects/crispr-barcode-cancer-metastasis/beam_dev/beam"

# EXAMPLE_DIR="${REPO_DIR}/examples/3457_Apc_T4"

# java  -Djava.library.path=$BEAGLE_LIB_PATH -Xmx10g -jar $REPO_DIR/beam.jar \
# -threads 5 \
# -overwrite \
# -working \
# -seed 2 $EXAMPLE_DIR/3.xml > $EXAMPLE_DIR/terminal_3.log