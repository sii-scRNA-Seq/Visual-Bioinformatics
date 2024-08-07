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

    pgrep -f back_end.py | awk '{print "kill -9 " $1}' | sh

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
    screen -S scampi -X readbuf $tempfile2
    screen -S scampi -X paste .
    rm -f $tempfile2

    echo "The app is being started in screen, you can monitor this via htop"

else

    echo "The given version does not exist"

fi
