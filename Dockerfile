FROM python:3.6-jessie

ENV HGVS_VERSION 1.0.5

RUN mkdir -p /opt \
  && cd /opt \
  && wget https://github.com/lifeomic/hgvs/archive/v${HGVS_VERSION}.tar.gz \
  && tar -xvzf v${HGVS_VERSION}.tar.gz \
  && rm v${HGVS_VERSION}.tar.gz \
  && cd hgvs-${HGVS_VERSION} \
  && python setup.py install

# vt ---------------------------------------------------------------------------
# -> /usr/local/bin/vt
RUN cd /opt \
  && git clone --depth 1 https://github.com/atks/vt \
  && cd /opt/vt \
  && make \
  && cp vt /usr/local/bin \
  && rm -rf /opt/vt

RUN mkdir -p /opt/app
WORKDIR /opt/app
COPY . /opt/app
RUN pip install -r requirements.txt
RUN cp /opt/hgvs-${HGVS_VERSION}/pyhgvs/data/refGene.hg19.txt /opt/app/refGene.hg19.txt

ENTRYPOINT ["python", "/opt/app/src/convert.py"]
