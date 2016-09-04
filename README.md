The *unRAVE* catalog
==================

[![Build Status](https://travis-ci.org/AnnieJumpCannon/RAVE.svg?branch=master)](https://travis-ci.org/AnnieJumpCannon/RAVE)

An independent data-driven analysis of *RAVE* spectra in anticipation of the *Gaia/TGAS* data release. You can access the **[most recent version of the PDF here](https://github.com/AnnieJumpCannon/RAVE/raw/master-pdf/article/unrave.pdf)**.


Instructions for making changes directly to the LaTeX
=====================================================

- [Sign up for a free GitHub account](https://github.com/join), then install `git` and [follow these instructions](https://help.github.com/articles/set-up-git/) to get your local `git` configuration set up with your GitHub account 

- Fork this repository on GitHub by clicking the 'Fork' button at the top of this page. Now you have your own parallel version of this repository on GitHub

- Create a local copy of your freshly-forked repository by using `git clone` in your terminal:

````
  git clone https://github.com/<your_username>/RAVE.git rave
````

- Make your changes to the `unrave.tex` document and ensure it compiles with `make`.

- Commit your changes using the following terminal commands (and if in doubt, commit often!):

````
git add article/unrave.tex
git commit -m "<short_summary_of_your_changes>"
````

- Now `push` your commited changes to your forked GitHub repository with this command:

````
git push
````

- Open a pull request on this repository by clicking the 'New pull request' button on the `AnnieJumpCannon/RAVE` repository page, which will ask me to accept the changes you have committed.


License
======= 
Copyright 2016 the authors. All rights reserved.
