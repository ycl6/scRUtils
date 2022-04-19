# Contributing to scRUtils development

Thank you for investing your time in contributing to scRUtils. This guide 
describes the contribution workflow for this package:

1. Filing a bug report or feature request in an issue.
1. Suggesting a change via a pull request.

scRUtils is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). 
By contributing to this project, you agree to abide by its terms.

## Issues

### Create a new issue

If you spot a problem with the package, search if an issue already exists. If 
a related issue doesn't exist, you can open a new issue using a relevant issue 
form.

When filing an issue, the most important thing is to include a minimal
reproducible example so that we can quickly verify the problem, and then figure
out how to fix it. There are three things you need to include to make your
example reproducible: 

1.  **required packages**: should be loaded at the top of the script, so it's 
    easy to see which ones the example needs.

1.  **data**: use `dput()` to generate the R code to recreate a data. For 
    example, to recreate the `mtcars` dataset in R, I'd perform the following 
    steps:

       1. Run `dput(mtcars)` in R
       2. Copy the output
       3. In my reproducible script, type `mtcars <- ` then paste.

     But even better is if you can create a `data.frame()` with just a handful
     of rows and columns that still illustrates the problem.

1.  **code**: Spend a little bit of time ensuring that your code is easy for 
    others to read, and use comments to indicate where your problem lies. Do 
    your best to make sure your code is concise and remove everything that is 
    not related to the problem.

You can check you have actually made a reproducible example by starting up a
fresh R session and pasting your script in.

### Solve an issue

Scan through the existing issues to find one that interests you. If you find an 
issue to work on, you are welcome to provide solutions or suggestions, either by 
replying to the issue itself and/or open a PR with a fix.

## Pull Request

To contribute a change to scRUtils, you follow these steps:

1. Fork the repository.
1. Create a working branch and start with your changes.
1. Commit the changes once you are happy with them.
1. Push the branch to your fork on GitHub and issue pull request (PR).
1. Discuss the pull request, and we may ask for changes to be made before a PR 
   can be merged.
1. Iterate until either the PR is accepted or decide that it's not a good fit.

