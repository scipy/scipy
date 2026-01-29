# SciPy Documentation

## How work with the docs

Have a look at
https://scipy.github.io/devdocs/dev/contributor/rendering_documentation.html

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
  - `release` contains the release notes. Note that those normally should not
    be updated as part of a PR; we keep release notes for the upcoming releases
    on the wiki of the main SciPy repo.
