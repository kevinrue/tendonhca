## 28 Feb 2023

* Used Visual Studio Code to open `/project/tendonhca` folder for the first time.
* Created README.txt to record activity.
* Installed VSCode Live Preview extension to preview files. Works for:
  * HTML: YES
  * EPS: NO
* Carla shared link to <https://github.com/Botnar-MSK-Atlas/hamstring_atlas>.
* 10X Visium
  * Overview <https://www.10xgenomics.com/products/spatial-gene-expression>
  * Software <https://support.10xgenomics.com/spatial-gene-expression/software/overview/welcome>

## 02 March 2023

* Created sub-directory `001-SpatialExperiment-test`.
* Start new script verifying that the package `SpatialExperiment` can be attached.

## 03 May 2023

* Created sub-directory `002-spaceranger-test`.
* Wrote job script processing an arbitrary sample.
* Created own symbolic links to updated path to FASTQ files.
* Job successfully completed overnight (04 May 2023).

## 09 May 2023

* Create a Conda environment for `cell2location` following instructions at <https://pypi.org/project/cell2location/>.
* Learn to load Visium data following tutorial at <https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_short_demo.html#1.-Loading-Visium-data>.

## 12 May 2023

```
snakemake --slurm --profile default
```

(snakemake) [albrecht@cbrglogin3 003-snakemake]$ snakemake --slurm --profile default
SyntaxError:
Input and output files have to be specified as strings or lists of strings.
  File "/ceph/project/tendonhca/albrecht/003-snakemake/workflow/Snakefile", line 17, in <module>