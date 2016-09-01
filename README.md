*The Cannon* with RAVE
======================

Some info.

Compile the manuscript
======================

From your terminal:

````
  git clone https://github.com/AnnieJumpCannon/RAVE.git rave 
  cd rave/article
  make
````

The manuscript will now be compiled to `unrave.pdf`

If you have used `git` before, then **please make your changes to the manuscript
directly in LaTeX and create a pull request**. There are instructions below on
how to do this. If you have questions, Google or Andy will answer.

If you are not familiar with
`git`, then please make your manuscript suggestions in the form of a 
[GitHub issue](https://www.github.com/AnnieJumpCannon/RAVE/issues/new).
That way your comments can be included (on short notice!) in parallel, and
everyone can see the contributions from each co-author.


Instructions for making changes directly to the LaTeX
=====================================================

0. [Sign up for a free GitHub account](https://github.com/join), then install `git` and [follow these instructions](https://help.github.com/articles/set-up-git/) to get your local `git` configuration set up with your GitHub account 

1. Fork this repository on GitHub by clicking the 'Fork' button at the top of this page. Now you have your own parallel version of this repository on GitHub

2. Create a local copy of your freshly-forked repository by using `git clone` in your terminal:

````
  git clone https://github.com/<your_username>/RAVE.git rave
````

3. Make your changes to the `unrave.tex` document and ensure it compiles with `make`.

4. Commit your changes using the following terminal commands:

````
git add article/unrave.tex
git commit -m "<short_summary_of_your_changes>"
````

If in doubt, commit often!

5. Now `push` your commited changes to your forked GitHub repository with this command:

````
git push
````

6. Open a pull request on this repository by clicking the 'New pull request' button on the `AnnieJumpCannon/RAVE` repository page, which will ask me to accept the changes you have committed.


License
======= 
Copyright 2016. All rights reserved.
