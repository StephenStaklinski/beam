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

which can be obtained from:
- [beast2](https://github.com/CompEvol/beast2)
- [beast-classic](https://github.com/BEAST2-Dev/beast-classic)
- [BeastFX](https://github.com/CompEvol/BeastFX)
- [BEASTLabs](https://github.com/BEAST2-Dev/BEASTLabs)
- [feast](https://github.com/tgvaughan/feast)
- [tidetree](https://github.com/seidels/tidetree)



