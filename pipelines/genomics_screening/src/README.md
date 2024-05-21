This directory contains the source code for filtering the annotated variants produced by GenomicScreening workflow \
and outputs a pdf report of those variants which passed the filtering criteria.\
Criteria are described in [Specification for Variant Reporting for M42 TRE Use Case 1](https://docs.google.com/document/d/1C0dBQIyuvU15CO3rXcIkDTKAtzwSxaWpVH53PS2ZtsE/edit#heading=h.r6j6zw6qnxcb)

The source code in this directory was copied from `pipelines/cohort_screening/src` and the `variant_report.py` script has been updated for
general genomic screening (not just MODY).


## Docker
The script is run inside a Docker container which is then pulled and run by the workflow.  \
The Dockerfile for variants_report.py is in this src/ directory.\
Currently, the Docker image is built and pushed manually to both Azure Container Registry and GCR. 
Instructions for building and pushing the image are in the [Dockerfile](Dockerfile).