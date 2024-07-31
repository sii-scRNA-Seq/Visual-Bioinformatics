# SCAMPI

Welcome to the SCAMPI project (Single Cell Analysis Methods Presented Interactively).

## Starting a development version of the application
### Starting the back-end
- Open an Anaconda terminal
- Navigate to the back-end directory
- Run `conda env create -f environment.yml`
- Run `conda activate scampi`
- Run `python back_end.py` to launch the back-end in development mode
### Starting the front-end
- Open a new terminal prompt
- Navigate to the front-end directory
- Run `npm install`
- Run `ng serve --open` to launch the front-end in development mode

## Compiling the application
- Run `npm run build`

## Starting the application in production

* Update package.json version number to release version
* Log into the production server (scampi.mvls.gla.ac.uk)
* Pull latest changes from `main` branch
* NPM install in frontend folder
* screen -r to reconnect to scampi session
* shut down server (ctrl-c)
* `$ conda deactivate`
* ctrl+a+d to come out of session
* `$ conda env remove scampi`
* `$ conda env create -f back-end/environment.yml`
* rebuild frontend
* copy into backend for serving
* Update `back_end.py` to run in production mode
* `$ screen -r`
* `$ python back-end/app/back_end.py`
* `ctrl+a+d`
* Verify on https://scampi.mvls.gla.ac.uk