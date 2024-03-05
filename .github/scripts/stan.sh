JULIA_CMDSTAN_HOME="/home/runner/cmdstan-2.29.2/"
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
cd $OLDWD
