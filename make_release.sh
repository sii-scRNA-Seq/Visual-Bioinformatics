#! /usr/bin/bash

version=$1

source ~/.bashrc
eval `ssh-agent -s`
ssh-add ~/.ssh/github_deployment_key
git reset --hard
git fetch
git checkout tags/$version
git pull

cd front-end/
npm install
npm run build

pgrep -f back_end.py | awk '{print "kill -9 " $1}' | sh

cd ..
rm -r back-end/dist
cp -r front-end/dist back-end/

tempfile2=$(mktemp)
cat > $tempfile2 <<EOF
conda deactivate
conda env remove -n scampi
conda env create -f ../environment.yml
conda activate scampi
python back_end.py --production-mode
EOF
screen -X readbuf $tempfile2
screen -X paste .
rm -f $tempfile2
