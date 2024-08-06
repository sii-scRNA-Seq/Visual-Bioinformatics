# SCAMPI

Welcome to the SCAMPI project (Single Cell Analysis Methods Presented Interactively).


## Starting the application in development

### Starting the back-end
* Open an Anaconda command line
* Navigate to the Visual-Bioinformatics directory
* Create and activate the new Conda environment
  * `$ conda env remove -n scampi --all`
  * `$ conda env create -f back-end/environment.yml`
  * `$ conda activate scampi`
* Launch the back-end in development mode
  * `$ cd back-end/app/`
  * `$ python back_end.py`

### Starting the front-end
* Open a new command line
* Navigate to the Visual-Bioinformatics directory
* Launch the back-end in development mode
  * `$ cd front-end/`
  * `$ npm install`
  * `$ ng serve --open`
* Verify that the application is running on http://127.0.0.1:5000/


## Starting the application in production

* Update package.json version number to release version
* Create a tag with release version
  * `$ git tag x.x.x`
  * `$ git push origin x.x.x`
* Log into the production server (scampi.mvls.gla.ac.uk)
* Navigate to the Visual-Bioinformatics directory
* Pull latest changes from `main` branch
  * `$ source ~/.bashrc`
  * `` $ eval `ssh-agent -s` ``
  * `$ ssh-add ~/.ssh/github_deployment_key`
  * `$ git reset --hard`
  * `$ git pull`
* Build the front-end
  * `$ cd front-end`
  * `$ npm install`
  * `$ npm run build`
* Shut down the existing scampi session
  * `$ screen -r`
  * `Ctrl + c`
  * `Ctrl + a + d`
* Move the new front-end to the back-end directory
  * `$ cd ..`
  * `$ rm -r back-end/dist`
  * `cp -r front-end/dist back-end/`
* Create and activate the new Conda environment
  * `$ screen -r`
  * `$ conda deactivate`
  * `$ conda env remove -n scampi --all`
  * `$ conda env create -f back-end/environment.yml`
  * `$ conda activate scampi`
* Launch the server
  * `cd back-end/app/`
  * `$ python back_end.py --production-mode`
  * `Ctrl + a + d`
* Verify that the application is running on https://scampi.mvls.gla.ac.uk
