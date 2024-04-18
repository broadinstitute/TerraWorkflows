FROM ubuntu:latest

RUN apt-get update \
    && apt-get install -y \
    wkhtmltopdf

RUN apt-get install -y \
    python3-pip idle3

RUN python3 -m pip install --no-cache-dir --upgrade pip

COPY requirements.txt .

RUN python3 -m pip install -r requirements.txt

COPY . /src

# builds with GitHub Action "Main Validation and Release" ../.github/workflows/build-and-publish-main.yaml
# tag = us-east4-docker.pkg.dev/broad-gotc-prod/variantreport:${{steps.get-artifact-slug.outputs.slug}}, us-east4-docker.pkg.dev/broad-gotc-prod/variantreport:latest

# For local Dev
# Google Artifact Registry
# docker build -t us-east4-docker.pkg.dev/broad-gotc-prod/variantreport:<new_version> .
# then push, if you want to push
# Azure Container Registry
# docker build -t terraworkflows.azurecr.io/variantreport:<new_version> .
# then tag, then push, if you want to push
# docker run --rm -it us-east4-docker.pkg.dev/broad-gotc-prod/variantreport:<new_version>

# to build and push to Google Artifact Registry
# make sure you are logged in to gcloud and that application default credentials are set (https://docs.google.com/document/d/1b03-YphH6Uac5huBopLYTYjzgDAlwS6qf-orMqaph64/edit?usp=sharing)
# set the VERSION field in update_docker_image.sh in this directory and then build and push as needed, for local development
# for Azure Container Registery - you need to az login and then az acr login --name terraworkflows
