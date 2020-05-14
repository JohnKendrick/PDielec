# Package Installation and upgrade protocol
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
git push                   # Check that the tags are valid on the remote github

```

Update github.com with a new release.  The process of releasing the code will give it a tag on the master branch, so issue a final pull on the local machine to download the new tag.

## Testing

### Setup.py
Edit setup.py and give the release a new version number, it is best to use a development number here as this can be changed.  Also add any new required modules to the install_requires option.

For example; version 6.3.1.0

Edit CHANGES file to record what changes have been made, don't use the developement number here, but use the final number you want used in the repository.

Run the following commands to remove the old installation create a new one;

    make pypy
    # rm -r build dist
    # python setup.py sdist bdist_wheel

This writes a PDielec/\_\_init\_\_.py file containing the version number.

### Install to test.pypi.org
You only need to do this if there have been significant changes to the packages required for installation.  If setup.py has not altered, or altered only in a small way then skip this step and go to *Final Installation*

    twine upload --repository testpypi dist/*

### Test installation in a new conda environment

	mkdir Test; cd Test
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

    make pypi
    # rm -r build dist
    # python setup.py sdist bdist_wheel
	twine upload --repository pypi dist/*
    # pip install PDielec

The make pypi command runs setup.py which creates a PDielec/\_\_init\_\_.py file with the version number in.

# Update the conda-forge installation
There is a bot (@regro-cf-autotick-bot) which seems to monitor pypi.org to see if there is a new release.  
This means that as long as meta.yaml has not changed there should be no reason to perform a manual installation
However, the manual work-flow for doing this is below, if it is necessary.  An alternative to the full manual process below is to clone the PR that the bot produces and edit the meta.yaml file there.

## Manual update of conda-forge installation  

Go the conda-forge pdielec feedstock entry at ;

    https://github.com/conda-forge/pdielec-feedstock

For the repository to your own github account, using the fork button in the top right. 
Clone a copy of your local fork of the feedstock.  I use a feedstock directory in my main Software/ directory.

    git clone git@github.com:JohnKendrick/pdielec-feedstock.git

If it has been a while since the fork was issued and you are returning to a local clone of the feedstock, issue the following;

    git checkout master
    git remote add upstream https://github.com/conda-forge/feedstock-pdielec
    git fetch upstream
    git rebase upstream/master

Create the changes on a new branch, use a branch name that reflects the version number.

    git checkout -b update_6_4_0

In pdielec-feedstock/recipe edit the meta.yaml file.  It has to have the right version number as recorderd in pypi.org.
The hash number for the package will also need updating.  You can get the hash number by going to the pypi.org.  
In the manage section for the new release there is an option against the PDielec-x.x.x.tar.gz file to view the hashes.  
Copy the sha256 hash and copy it into the meta.yaml file.
If there have been any changes to the dependencies that update the dependency list.

Finally check the build number, it should be 0 for a new release.  It only needs to change if changes were made to the packages but the source code has stayed the same.

Review all the changes and push the fork onto GitHub.

    git status
	git push origin update_6_4_0

Now create a pull request on https://github.com/JohnKendrick/pdielec-feedstock




# Initial Installing to CONDA - ONLY DO THIS IF THE FEEDSTOCK DOES NOT EXIST

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
# switch to conda environment ?????
conda skeleton pypi PDielec
# This will create a new meta.yaml file
# Edit the yaml file
git commit -a
git push --set-upstream origin pdielec_3.6.1   # Not sure about this
git status
```

