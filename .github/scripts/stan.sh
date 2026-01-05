#!/bin/bash
set -e

# Use $HOME to support both GitHub-hosted and self-hosted runners
JULIA_CMDSTAN_HOME="$HOME/cmdstan-2.34.1/"
OLDWD=`pwd`

cd ~
pwd
wget https://github.com/stan-dev/cmdstan/releases/download/v2.34.1/cmdstan-2.34.1.tar.gz
tar -xzpf cmdstan-2.34.1.tar.gz
ls -lia .
ls -lia ./cmdstan-2.34.1
ls -lia ./cmdstan-2.34.1/make
touch ./cmdstan-2.34.1/make/local
echo "STAN_THREADS=true" > ./cmdstan-2.34.1/make/local
make -C $JULIA_CMDSTAN_HOME build

# Export to GITHUB_ENV so subsequent steps can access it
if [ -n "$GITHUB_ENV" ]; then
    echo "JULIA_CMDSTAN_HOME=$JULIA_CMDSTAN_HOME" >> $GITHUB_ENV
fi

cd $OLDWD
echo "CmdStan installed at: $JULIA_CMDSTAN_HOME"
