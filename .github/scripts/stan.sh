JULIA_CMDSTAN_HOME="/home/runner/cmdstan-2.29.2/"
OLDWD=`pwd`
cd ~
pwd
wget https://github.com/stan-dev/cmdstan/releases/download/v2.29.2/cmdstan-2.29.2.tar.gz
tar -xzpf cmdstan-2.29.2.tar.gz
ls -lia .
ls -lia ./cmdstan-2.29.2
ls -lia ./cmdstan-2.29.2/make
touch ./cmdstan-2.29.2/make/local
echo "STAN_THREADS=true" > ./cmdstan-2.29.2/make/local
make -C $JULIA_CMDSTAN_HOME build
cd $OLDWD
