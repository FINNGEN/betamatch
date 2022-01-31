FROM ubuntu:18.04
RUN apt-get update && apt-get install -qqy wget unzip bzip2 curl python3.6 python3-venv python3-pip git nano \
&& apt-get install curl make tabix python3 python3-pip zlib1g-dev libjpeg-dev --yes && \
    apt-get clean && \
    pip3 install pytabix numpy pandas plotnine statsmodels
COPY betamatch.py /usr/local/bin
COPY corrplot.py /usr/local/bin
COPY beta_utils.py /usr/local/bin
RUN chmod +x /usr/local/bin/betamatch.py
RUN chmod +x /usr/local/bin/corrplot.py
