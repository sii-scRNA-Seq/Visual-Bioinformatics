# SCAMPI

## Activating the Conda environment
* Open an Anaconda command line
* Navigate to the Visual-Bioinformatics directory
* Create and activate the new Conda environment
  * `$ conda env remove -n scampi --all`
  * `$ conda env create -f back-end/environment.yml`
  * `$ conda activate scampi`

## Starting the back-end in development
* `$ cd back-end/app/`
* `$ python back_end.py`

## Unit testing the back-end
* `$ cd back-end/app/`
* `$ python -m pytest`

## Linting the back-end
* `$ cd back-end/`
* `$ flake8 --ignore=E501 app`
