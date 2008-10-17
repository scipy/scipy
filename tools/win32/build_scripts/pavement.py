options(
    setup=Bunch(
        name = "scipy-superpack",
    )
)

@task
def setup():
    print "Setting up package %s" % options.name
