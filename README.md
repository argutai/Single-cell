# Single Cell analysis pipeline 
Specifically for BD Rhapsody data, using the BD pipeline as documented [here](https://scomix.bd.com/hc/en-us/articles/360023044512-BD-AbSeq-Bioinformatics-Guides).


### Converting Fastqs

Pull the BD pipeline docker image and download the BD human reference transcriptome

```
singularity pull fastqc.sif docker://bdgenomics_rhapsody:2.0.sif (run using srun if on KCL CREATE HPC)
wget http://bd-rhapsody-public.s3-website-us-east-1.amazonaws.com/Rhapsody-WTA/Pipeline-version2.x_WTA_references/
```

Add the paths to your forward and reverse read fastq files to Reads.file.location, and the refernce transcriptome to Reference_Archive.file.location and submit the job with `sbatch -p cpu job.sh` on KCL CREATE HPC.

This takes about 24 hours for f/r file sizes of ~6GB and returns BAM files, a html report, and a seurat/scanpy object.

### Quality control
Scp the Pipeline_Report.htmls, the Seurat.rds, and the h5ad files to your local. Take a look at the summary stats from the BD-pipeline on the html files. run `python QC.py path/to/h5ad` to apply quality control on scanpy object. 

You may find later that some scanpy commands do not work on the h5ad object, I found that converting the Seurat.rds to h5ad using the convert.R file creates a file which does not experience these problems.

