# We should create an Open Babel implementation of pipeline_utils that handles the basic I/O for
# structure files so that the dependency on RDKit can be removed.
# See https://github.com/InformaticsMatters/pipelines-obabel/issues/1

FROM informaticsmatters/obabel:2.4.1
LABEL maintainer="Tim Dudgeon<tdudgeon@informaticsmatters.com>"

USER root

# Copy the obabel pipeline implementation into the image
COPY src/python /opt/python-obabel
RUN apt-get install -y --no-install-recommends python-setuptools gzip python-pip && pip install -e /opt/python-obabel


