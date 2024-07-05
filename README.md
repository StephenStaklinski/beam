# Metastabayes BEAST2 package for Bayesian cancer metastasis graph inference from CRISPR cell lineage tracing data

### Compilation

Use `ant compile` to compile the `.java` files to `.class` files if the `bin/` directory does not exist. Then, use `ant jar` to build a runnable jar file if the `metastabayes.jar` file does not exist

## Inference

Use `java -jar file.xml` to run the metastabayes package with a properly formatted xml file, such as the one provided in `examples/`.

