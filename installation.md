# Installation and upgrade protocol
Edit ~/.pypirc file to hold usernames for the repositories being used.

```
[distutils]
index-servers = 
	pypi
	testpypi

[pypi]
username: john_kendrick

[testpypi]
repository: https://test.pypi.org/legacy/
username: john_kendrick
```

Use keyring to store the passwords for the repositories.  Note, in the examples below there may be more repos than neccessary.

    keyring set https://test.pypi.org/legacy/ john_kendrick
    keyring set https://upload.pypi.org/simple/ john_kendrick
    keyring set https://upload.pypi.org/legacy/  john_kendrick
    keyring set https://upload.pypi.org  john_kendrick

## Git

Update the CHANGES file to reflect the version number.  If you are working on a development branch, merge the changes.

```
git status       # On development branch to make sure there are no pending commits needed
git checkout master
git pull
git merge --no-ff develop  # --no-ff keeps a track of all the commits in develop
git branch -d develop
git tag v6.3.1             # This should agree with what is in the CHANGES file

```

## Testing

### Setup.py
Edit setup.py and give the release a new version number, it is best to use a development number here as this can be changed.  Also add any new required modules to the install_requires option.

For example; version 6.3.1.dev0

Edit CHANGES file to record what changes have been made, don't use the developement number here, but use the final number you want used in the repository.

Run the following commands to remove the old installation create a new one;

    rm -r build dist
    python setup.py sdist bdist_wheel

### Install to test.pypi.org

    twine upload --repository testpypi dist/*

### Test installation in a new conda environment

    conda create -n test python
	conda activate test
    pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ PDielec
	pdgui

### Remove test installation

pip installs pdgui and preader command scripts into ~./local/bin and the PDielec modules are put into ~/.local/lib
To uninstall;

	pip uninstall PDielec

To remove the conda testing environment;
 
	conda deactivate
	conda env remove --name test

## Final Installation

Edit the setup.py file and remove the development designation from the project version

    rm -r build dist
    python setup.py sdist bdist_wheel
	twine upload --repository pypi dist/*
    # pip install PDielec


# Initial Installing to CONDA

Make sure .condarc is;

```
channel_priority: strict
channels:
  - conda-forge
  - defaults
```

Got https://https://github.com/conda-forge/staged-recipes
fork the repository

```
cd ~/Software
git clone git@github.com:conda-forge/staged-recipes.git
git checkout -b pdielec_6.3.1
```

make a new directory in recipes/ (called pdielec)

```
cd pdielec
switch to conda environment ?????
conda skeleton pypi PDielec
This will create a new meta.yaml file
Edit the yaml file
git commit -a
git push --set-upstream origin pdielec_3.6.1   # Not sure about this
git status
```

