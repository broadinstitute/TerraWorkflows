FROM python:3.8-slim

RUN pip install cyp2d6_parser
RUN pip install requests

COPY stargazer_parser.py .
COPY merge_calls.py .
COPY stargazer_converter.csv .