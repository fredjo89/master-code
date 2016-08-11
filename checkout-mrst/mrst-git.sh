#!/bin/bash

mkdir mrst-bitbucket
cd mrst-bitbucket

repos='core autodiff multiscale internal solvers model-io visualization'
for r in $repos; do
 # git clone git@bitbucket.org:mrst/"mrst-${r}.git"
  git clone https://fredjo89@bitbucket.org/mrst/"mrst-${r}.git"
done

cp ../startup_user.m mrst-core/

