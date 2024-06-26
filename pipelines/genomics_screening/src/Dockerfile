FROM --platform=linux/amd64 python:3.9.2-slim-buster

RUN apt-get update \
    && apt-get install -y \
    wkhtmltopdf

RUN apt-get install -y \
    python3-pip

RUN python3 -m pip install --no-cache-dir --upgrade pip

COPY requirements.txt .

RUN pip3 install --no-cache-dir -r requirements.txt

COPY . /src

# Will build with GitHub Action "Main Validation and Release" ../.github/workflows/build-and-publish-main.yaml
# tag = us.gcr.io/broad-gotc-prod/genomics_variant_report:${{steps.get-artifact-slug.outputs.slug}}, us.gcr.io/broad-gotc-prod/genomics_variant_report:latest

# For now
# Google Artifact Registry
# to tag with the latest commit short SHA
# docker build --platform="linux/amd64" -t us.gcr.io/broad-gotc-prod/genomics_variant_report:<latest_commit_short_SHA> .
# docker build --platform="linux/amd64" -t us.gcr.io/broad-gotc-prod/genomics_variant_report:latest .
# Azure Container Registry
# docker build --platform="linux/amd64" -t terraworkflows.azurecr.io/genomics_variant_report:<latest_commit_short_SHA> .
# docker build --platform="linux/amd64" -t terraworkflows.azurecr.io/genomics_variant_report:latest .

# note that to run locally, assuming you are on a Mac, you will likely need to build without the platform flag

# docker run --rm -it us.gcr.io/broad-gotc-prod/genomics_variant_report:<new_version>

# to push to Google Artifact Registry
# make sure you are logged in to gcloud
# gcloud auth login
# gcloud auth configure-docker us.gcr.io
# and that application default credentials are set (https://docs.google.com/document/d/1b03-YphH6Uac5huBopLYTYjzgDAlwS6qf-orMqaph64/edit?usp=sharing)
# set the VERSION field in update_docker_image.sh in this directory and then build and push as needed
# for Azure Container Registery - you need to az login and then az acr login --name terraworkflows
