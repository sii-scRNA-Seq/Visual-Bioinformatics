# SCAMPI

## Activating the Conda environment
In an Anaconda command line, navigate to Visual-Bioinformatics/back-end/ and run:
* `$ conda env remove -n scampi -y`
* `$ conda env create -f environment.yml`
* `$ conda activate scampi`

## Starting the back-end in development
From Visual-Bioinformatics/back-end/ run:
* `$ cd app/`
* `$ python back_end.py`

## Unit testing the back-end
From Visual-Bioinformatics/back-end/ run:
* `$ cd app/`
* `$ python -m pytest`

## Linting the back-end
From Visual-Bioinformatics/back-end/ run:
* `$ flake8 --ignore=E501 app`
