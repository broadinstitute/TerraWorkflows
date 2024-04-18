FROM ubuntu:latest

RUN apt-get update \
    && apt-get install -y \
    wkhtmltopdf

RUN apt-get install -y \
    python3-pip idle3

RUN python3 -m pip install --no-cache-dir --upgrade pip

COPY reporting/requirements.txt .

RUN python3 -m pip install -r requirements.txt

COPY src /src