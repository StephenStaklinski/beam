# Metastabayes BEAST2 package for Bayesian cancer metastasis graph inference from CRISPR cell lineage tracing data

## Inference

Install java+fx 17.0.9 and BEAGLE and set the path to it in `BEAGLE_LIB_PATH`. Then, the provided `metastabayes.jar` file is all that is needed to run the method with:

```
# specify path to the repo
REPO_DIR="/Users/staklins/projects/crispr-barcode-cancer-metastasis/metastabayes_dev/metastabayes"

# submit example data to run
java -Djava.library.path=$BEAGLE_LIB_PATH -jar metastabayes.jar \
-seed 1724005224593 \
-D repo=$REPO_DIR \
-overwrite \
-working examples/pR_298.xml
```

## Development

An initial install will obtain dependencies for developing on the code base by running the following with the path to Java FX specified:
```
# make the working directory
mkdir metastabayes_dev
cd metastabayes_dev

# download dependencies
gh repo clone StephenStaklinski/metastabayes
gh repo clone CompEvol/beast2
gh repo clone BEAST2-Dev/beast-classic
gh repo clone CompEvol/BeastFX
gh repo clone BEAST2-Dev/BEASTLabs
gh repo clone seidels/tidetree
gh repo clone tgvaughan/feast

# specify path to Java FX and compile feast code, can also be an OpenJDK install
cd feast
JAVA_FX_HOME="/Users/staklins/bin/zulu17.46.19-ca-fx-jdk17.0.9-macosx_aarch64/lib" ant

# compile metastabayes code and related dependencies
cd ../metastabayes
ant compile
ant jar
```

Then, any further changes in the code base can be reflected by running the following in the `metastabayes_dev/metastabayes` directory:
```
ant compile
ant jar
```



