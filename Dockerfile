# TODO - this image is a bit of a mongrel as it needs RDKit because of the need for pipelines_utils
# We should create an Open Babel implementation of pipeline_utils that handles the basic I/O for
# structure files so that the dependency on RDKit can be removed.
# See https://github.com/InformaticsMatters/pipelines-obabel/issues/1

FROM informaticsmatters/rdkit_pipelines:latest
LABEL maintainer="Tim Dudgeon<tdudgeon@informaticsmatters.com>"

USER root

RUN apt-get -y update && apt-get install -y --no-install-recommends openbabel python-openbabel

# Copy the obabel pipeline implementation into the image
COPY src/python /opt/python-obabel
RUN pip install -e /opt/python-obabel

# and the pipeline-utilities
RUN pip install im-pipelines-utils==1.0.*

