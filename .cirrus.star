# The guide to programming cirrus-ci tasks using starlark is found at
# https://cirrus-ci.org/guide/programming-tasks/
#
# In this simple starlark script we simply check conditions for whether
# a CI run should go ahead. If the conditions are met, then we just
# return the yaml containing the tasks to be run.

load("cirrus", "env", "fs", "http")

def main(ctx):
    ######################################################################
    # Should wheels be built?
    # Only test on the scipy/scipy repository
    # Test if the run was triggered by:
    # - a cron job called "nightly". The cron job is not set in this file,
    #   but on the cirrus-ci repo page
    # - commit message containing [wheel build]
    # - a tag that begins with v, but doesn't end with dev0
    ######################################################################

    if env.get("CIRRUS_REPO_FULL_NAME") != "scipy/scipy":
        return []

    if env.get("CIRRUS_CRON", "") == "nightly":
        return fs.read("ci/cirrus_wheels.yml")

    tag = env.get("CIRRUS_TAG")

    # test if it's a tag?
    if tag:
        # if tag name contains "dev", then don't build any wheels
        # I tried using a regex for this, but I don't think go supports
        # the particular regex I wanted to use
        if "dev" in tag:
            return []
        elif tag.startswith("v"):
            # all other tags beginning with 'v' will build. cirrus_wheels.yml
            # will only upload wheels to staging if the following bash test is
            # True:
            # [[ $CIRRUS_TAG =~ ^v.*[^dev0]$ ]]
            return fs.read("ci/cirrus_wheels.yml")

    # Obtain commit message for the event. Unfortunately CIRRUS_CHANGE_MESSAGE
    # only contains the actual commit message on a non-PR trigger event.
    # For a PR event it contains the PR title and description.
    SHA = env.get("CIRRUS_CHANGE_IN_REPO")
    url = "https://api.github.com/repos/scipy/scipy/git/commits/" + SHA
    dct = http.get(url).json()
    if "[wheel build]" in dct["message"]:
        return fs.read("ci/cirrus_wheels.yml")

    return []
