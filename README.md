# continuous-methylation

A bioinformatics pipeline to get methylation signals with [modkit](https://github.com/nanoporetech/modkit/) (in bed format) for phased alignments using long-read sequencing (third-gen) data.

This pipeline accepts only hifi OR ont, and it can 1) generate phased alignments using a phased assembly if parental illuminas are available, 2) generate phased alignments using heterozygous SNVs + SVs and will output a single bam that you can query using "HP:i:1 or HP:i:2 or ![HP]", 3) [only for autism cohort] perform differential methylation.

## Getting started
1. Clone the repo
2. Install Snakemake version >= 7.6.0
3. Fill in the config and manifests.
4. Start the analysis!
    * Begin with a dry-run
    ```
    ./runlocal 10 -np
    ```
    * If dry-run looks good, proceed with:
    ```
    ./runlocal 10
    ```

## To-do
- [ ] Put in CI tests
- [ ] Add conda enviroments
- [ ] Build bare-bone container to get minimal example running

## FAQ
1. What if I want to change the version of certain tools?
   1. You could change the defaults in the tools.yaml config, and if you want to change the resources needed for certain steps- then please edit the resources.yaml
   * 
   ```
   ./config
   ├── config.yaml
   ├── manifest.tab
   ├── resources.yaml
   └── tools.yaml
   ```
2. How do I use different tech?
   1. Currently, there is only support for `hifi` or `ont`. You could use the pipeline like the following:
   * 
   ```
   # ONT
   ./runsnake 30 --config tech=ont manifest=config/manifest.tab (for dry run: put -np)
    
   # HiFi
   ./runsnake 30 --config tech=hifi manifest=config/manifest-hifi.tab (for dry run: put -np)
   ```
3. Where can I get the singularity containers?
   1. Answer coming soon. I'm going to build a barebone container that is enough to run the pipeline, as well as adding a conda statement for all the rules so it can be fine tuned on the versioning.
4. What is an example config and manifest?
   1. Answer coming soon.

## Overview
![pipeline vector](https://github.com/projectoriented/continuous-methylation/blob/main/agg-dag.svg)
