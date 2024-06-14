# SCAMPI

## Starting the back-end in development mode
- Open an Anaconda terminal
- Make sure you are in the back-end directory
- Run `conda env create -f environment.yml`
- Run `conda activate scampi`
- Run `cd app`
- Run `python back_end.py` to launch the back-end on http://127.0.0.1:5000/

## Unit testing the back-end
- Run `cd app`
- Run `python -m pytest`

## Linting the back-end
- Run `flake8 --ignore=E501`
