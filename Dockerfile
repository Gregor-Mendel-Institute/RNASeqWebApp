FROM python:2
LABEL maintainer="Uemit Seren <uemit.seren@gmail.com>"

RUN useradd --uid=10372 -ms /bin/bash gwas-web

RUN mkdir -p /net/gmi.oeaw.ac.at/gwasapp
VOLUME '/net/gmi.oeaw.ac.at/gwasapp/'

WORKDIR /rnaseq-web-app
COPY backend/requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

RUN chown -R gwas-web:root /rnaseq-web-app && chmod -R 755 /rnaseq-web-app

# required because of https://github.com/docker/docker/issues/20240
RUN cd backend && ln -s static ../frontend/target/rnaseqwebapp-0.0.1-SNAPSHOT/ && cd public && ln -s resources ../../frontend/target/rnaseqwebapp-0.0.1-SNAPSHOT/WEB-INF/classes/com/gmi/rnaseqwebapp/client/resources/ && cd /rnaseq-web-app

ENV PYTHONPATH=/rnaseq-web-app/atgwas/src/:/rnaseq-web-app/
ENV GOTO_NUM_THREADS=1
ENV OMP_NUM_THREADS=1
USER gwas-web

EXPOSE 8090

ENTRYPOINT ["python" , "backend/__init__.py"]