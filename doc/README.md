# SciPy Documentation

## How to build the docs

To build the html docs for local development, SciPy itself needs to be built so your
environment needs to be set up for that.  For details on that, see the
[Contributor Guide](http://scipy.github.io/devdocs/dev/contributor/contributor_toc.html#development-environment)).

Also ensure to initialize and update submodules (this pulls in the SciPy Sphinx
theme and `numpydoc`):
```
git submodule update --init
```

Now to build both SciPy itself and the docs, use:
```
python3 runtests.py --doc html
```

Alternatively, if you prefer to build SciPy and the docs separately rather
than use `runtests.py`:
```
python setup.py develop  # in the root of the repo
cd doc && make html-scipyorg
```

In case the SciPy version found by the above command is different from that of the
latest commit in the repo, you will see a message like:
```
installed scipy 5fd20ec1aa != current repo git version '35fd20ec1a'
```

This indicates that you're likely picking up the wrong SciPy install, check
with `python -c "import scipy; print(scipy.__file__)"`.

If the build is successful, you can open it in your browser with `make show`
(which will open `build/html-scipyorg/index.html`).


## Building pdf docs

To build the pdf docs, which requires a LaTeX install and can be more fiddly
to get to work, replace the doc build commands in the section above with:
```
python3 runtests.py --doc latex
```
or:
```
make latex
```

That will use Sphinx to generate the LaTeX sources. To then produce a pdf,
navigate to `doc/build/latex/` and run:
```
make all-pdf
```

That will produce a file `scipy-ref.pdf` in `build/latex/`.


## Building documentation for a release

For building all the documentation artifacts for a release, run:
```
make dist
```

This will build SciPy in-place (to ensure the version is correct), build html
and pdf docs as well as create a zip archive of the html docs that can easily
be redistributed.


## Layout of the docs in this repository

- `source` is where most of the content lives.
  - `dev` contains the contributor and developer guides as well as the governance
    docs and the code of conduct.
  - `tutorial` contains all tutorial content.
- `release` contains the release notes. Note that those normally should not be
  updated as part of a PR; we keep releases notes for the upcoming releases
  on the wiki of the main SciPy repo.