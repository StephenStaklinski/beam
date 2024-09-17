# Metastabayes BEAST2 package for Bayesian cancer metastasis graph inference from CRISPR cell lineage tracing data

## Inference

Install java+fx 17.0.9 and BEAGLE and set the path to it in `BEAGLE_LIB_PATH`.

```
REPO_DIR="/Users/staklins/projects/crispr-barcode-cancer-metastasis/metastabayes_dev/metastabayes"

java -Djava.library.path=$BEAGLE_LIB_PATH -jar metastabayes.jar \
-seed 1724005224593 \
-D repo=$REPO_DIR \
-overwrite \
-working examples/pR_298.xml
```

## Development

I have been using `build.xml` to run `ant compile && ant jar` after developing some changes in the code base. This `build.xml` will require the following repos for dependencies to be in the parent directory as such:
```
../beast2
../beast-classic
../BeastFx
../BEASTLabs
../feast
../tidetree
```

An initial install will obtain these by running the following with the path to Java FX specified:
```
mkdir metastabayes_dev
cd metastabayes_dev
gh repo clone StephenStaklinski/metastabayes
gh repo clone CompEvol/beast2
gh repo clone BEAST2-Dev/beast-classic
gh repo clone CompEvol/BeastFX
gh repo clone BEAST2-Dev/BEASTLabs
gh repo clone tgvaughan/feast
gh repo clone seidels/tidetree
cd feast
JAVA_FX_HOME="/Users/staklins/bin/zulu17.46.19-ca-fx-jdk17.0.9-macosx_aarch64/lib" ant
cd ../metastabayes
ant compile
ant jar
```



