#!/usr/bin/env bash
# run locally with something like this:
# ./post-service-descriptors.sh http://localhost:8091/coreservices/rest/v1/services
# or
# docker run -it --rm -v $PWD:$PWD:Z -w $PWD --network deploy_squonk_back centos:7 ./post-service-descriptors.sh

set -e

POST=${1:-http://coreservices:8080/coreservices/rest/v1/services}
BASE_D='docker://github.com/InformaticsMatters/pipelines-obabel'
CT_DJ="application/x-squonk-service-descriptor-docker+json"
CT_MM="multipart/mixed"


for d in 'src/python/pipelines_obabel/docking'
do
    for file in $d/*.dsd.yml
    do
	    echo $file
	    curl -X POST \
         -T $file\
         -H "Content-Type: $CT_DJ"\
         -H "Base-URL: $BASE_D"\
         $POST
         echo ""
    done
done
