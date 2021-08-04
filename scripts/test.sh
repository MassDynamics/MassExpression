#!/bin/bash


RENV_PATHS_CACHE_HOST=$(pwd)/renv/cache
# the path to the cache within the container
RENV_PATHS_CACHE_CONTAINER=/usr/local/src/MassExpression/renv/cache

docker run -v "${RENV_PATHS_CACHE_HOST}:${RENV_PATHS_CACHE_CONTAINER}"  -e "RENV_PATHS_CACHE=${RENV_PATHS_CACHE_CONTAINER}" mass_expression R -e "renv::restore();sessionInfo();devtools::check()"

ls renv
ls -al $RENV_PATHS_CACHE_HOST
