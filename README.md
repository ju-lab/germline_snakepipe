# germline-snakepipe
Snakemake pipeline to process FASTQ -> Processed BAM -> basic germline SNV calling

Calls are made with `freebayes` with joint calling methods with following filters added. 

`--standard-filters -j --min-coverage 10 -F 0.2 -C 2 --read-snp-limit 3 --read-mismatch-limit 3 --ploidy 2`

Multisample VCF output of `freebayes` is then post-processed with `vcfallelicprimitives` to decompose complex events into SNVs and indels, and normalized with `vt`. 

# How to run 
```
$ snakemake -p done --profile profile -j 24 --immediate-submit
```

# How to configure `sampleConfig.yaml`

```yaml
project:
    [projectname]

samples:
    proband:
        normal_fastq: /path/to/fastq/prefix_R
        sex: XX
        status: 
          - diseased
    healthy1:
        normal_fastq: ../../../kjyi/Projects/temp_20180425_LeeEu
        sex: XY
        status:
          - healty
    ...
    ...

```

`project`: defines the project name
`normal_fastq`: path to fastq.gz minus '1_fastq.gz' string -> change to `normal_fastq1` `normal_fastq2` will make it more intuitive. 
`sex`: sex of the patient (not used for now)
`status`: status of the patient (not used for now)


