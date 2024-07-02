To get the `array-api` test suite, running one needs to do a few things.  First, you will need to get the [array-api test suite](https://github.com/data-apis/array-api-tests) and install this package into that environment:

```shell
git clone https://github.com/data-apis/array-api-tests.git
cd array-api-tests
# Simplest to just get the dev dependencies of scipy since we will build from source. 
# Also if you have a dev environment already created, you may will want to change the name in this environment.yml.
# Alternatively, you can use that dev environment directly and just skip to installing the array-api requirements.txt.
micromamba env create -f path/back/to/scipy/environment.yml
pip install -e path/back/to/scipy --no-build-isolation
pip install -r requirements.txt
```

Next, you will need to replace the runnable by setting the following variable:

```shell
export ARRAY_API_TESTS_MODULE=scipy.sparse._array_api
```

Once this is done, you can run the test suite via:

```shell
pytest array_api_tests --skips-file path/back/to/scipy/scipy/sparse/array_api/array-api-skips.txt
```

and make edits in real time (since we installed `scipy` in editable mode) which will then be tested for.

