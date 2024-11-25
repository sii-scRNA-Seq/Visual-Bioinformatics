#! /usr/bin/bash

version=$1

source ~/.bashrc
eval `ssh-agent -s`
ssh-add ~/.ssh/github_deployment_key
git reset --hard
git fetch

if [ $(git tag -l "$version") ]; then

git checkout tags/$version
git pull

cd front-end/
npm install
npm run build
cd ..

# Stop app
pgrep -f back_end.py | awk '{print "kill -9 " $1}' | sh
pgrep -f redis | awk '{print "kill -9 " $1}' | sh
pgrep -f celery | awk '{print "kill -9 " $1}' | sh

# Stop conda env
tempfile2=$(mktemp)
cat > $tempfile2 <<EOF
conda deactivate
EOF
screen -S scampi -X readbuf $tempfile2
screen -S scampi -X paste .
screen -S redis -X readbuf $tempfile2
screen -S redis -X paste .
screen -S celery -X readbuf $tempfile2
screen -S celery -X paste .
rm -f $tempfile2

rm -r back-end/dist
cp -r front-end/dist back-end/

# Rebuild conda environment
conda env remove -n scampi -y
conda env create -f back-end/environment.yml

# Restart Redis
# WD: back-end
tempfile2=$(mktemp)
cat > $tempfile2 <<EOF
conda activate scampi
chmod +x start_redis.sh
./start_redis.sh
EOF
screen -S redis -X readbuf $tempfile2
screen -S redis -X paste .
rm -f $tempfile2

# Restart Celery Worker
# WD: back-end/app
tempfile2=$(mktemp)
cat > $tempfile2 <<EOF
conda activate scampi
chmod +x ../start_celeryworker.sh
../start_celeryworker.sh
EOF
screen -S celery -X readbuf $tempfile2
screen -S celery -X paste .
rm -f $tempfile2

# Restart Server
# WD: back-end/app
tempfile2=$(mktemp)
cat > $tempfile2 <<EOF
conda activate scampi
python back_end.py --production-mode
EOF
screen -S scampi -X readbuf $tempfile2
screen -S scampi -X paste .
rm -f $tempfile2

echo "The app is being started in screen, you can monitor this via htop"

else

echo "The given version does not exist"

fi
