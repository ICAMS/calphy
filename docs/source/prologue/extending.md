
# Support, contributing and extending

calphy welcomes and appreciates contribution and extension to the
module. Rather than local modifications, we request that the
modifications be submitted through a pull request, so that the module
can be continuously improved.

## Reporting and fixing bugs

In case a bug is found in the module, it can be reported on the [issues
page of the repository](https://github.com/ICAMS/calphy/issues). Once a bug is reported, the status can once again monitored on
the issues page. Additionally, you are of course very welcome to fix any
existing bugs.

## New features

If you have an idea for new feature, you can submit a feature idea
through the [issues page of the
repository](https://github.com/ICAMS/calphy/issues). As much as
information as you can provide about the new feauture would be greatly
helpful. Additionally, you could also work on feature requests already
on the issues page. The following instructions will help you get started
with local feature development.

### Setting up local environment

1.  The first step is to fork pyscal. A detailed tutorial on forking can
    be found [here](https://help.github.com/en/articles/fork-a-repo).
    After forking, clone the repository to your local machine.
2.  We recommend creating a virtual environment to test new features or
    improvements to features. See this
    [link](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
    for help on managing environments.
3.  Once the environment is set up, you can create a new branch for your
    feature by `git checkout -b new_feauture`.
4.  Now implement the necessary feature.
5.  Once done, you can reinstall pyscal by `python setup.py install`.
    After that please make sure that the existing tests work by running
    `pytest tests/` from the main module folder.
6.  If the tests work, you are almost done! If the new feature is not
    covered in existing tests, you can to write a new test in the tests
    folder. calphy uses pytest for tests. [This
    link](http://doc.pytest.org/en/latest/getting-started.html) will
    help you get started.
7.  Add the necessary docstrings for the new functions implemented.
    calphy uses the [numpy docstring
    format](https://numpydoc.readthedocs.io/en/latest/format.html) for
    documentation.
8.  Bonus task: Set up few examples that document how the feature works
    in the `docs/source/` folder and link it to the examples section.
9.  Final step - Submit a pull request through github. Before you
    submit, please make sure that the new feature is documented and has
    tests. Once the request is submitted, automated tests would be done.
    Your pull request will fail the tests if - the unit tests fail, or
    if the test coverage falls below 80%. If all tests are successful,
    your feauture will be incorporated to calphy and your contributions
    will be credited.

If you have trouble with any of the steps, or you need help, please
[send an email](mailto:s.menon@mpie.de) and we will be happy to
help! All of the contributions are greatly appreciated and will be
credited in Developers/Acknowledgements page.

