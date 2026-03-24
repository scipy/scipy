# Running custom CodeQL queries in code scanning advanced setup

For running [CodeQL queries](.github/codeql-queries/scipy-c-cpp-queries/signed_shift.ql), we need CodeQL packs. CodeQL packs are used to package and version queries. There's already a CodeQL pack for C/C++ queries, so there's no need to create another one. If you'd ever want to write queries for another language, you'd need to create a CodeQL for that query.

[CodeQL query suites](.github/codeql-queries/scipy-code-scanning.qls) are used to specify which queries should be run, in case there are some queries you’d prefer not to run.

## Write and test a CodeQL query
The easiest way to write and run a query is by using VS Code with the VS Code CodeQL extension.

After installing the CodeQL extension:
- execute the command (with Ctrl/Cmd + Shift + P) "CodeQL: Create Query",
- then, choose a language (e.g. Python) in the dropdown that will appear. This will create a folder with an example query and a CodeQL pack.

To run the query, you'll need a CodeQL database of a codebase.
- Open the example query `example.ql`, right-click and choose "CodeQL: Run Query on Selected Database".
- You'll see a dropdown with options to select a database. Choose "Download from GitHub" and type `scipy/scipy`. This will download the prebuilt CodeQL database, and run the query on it.

To download a CodeQL database for another lanugage of the same codebase:
- execute the command (with Ctrl/Cmd + Shift + P) "CodeQL: Download Database"
- then, go to the CodeQL VS Code extension tab, "Databases" section, and choose the database you'd like to run the query against.

For running preexising queries, like the ones in `.github/codeql-queries/scipy-c-cpp-queries` folder, go directly to the query you'd like to run, and:
- execute the command (with Ctrl/Cmd + Shift + P) "CodeQL: Install Pack Dependencies"
- right-click and choose "CodeQL: Run Query on Selected Database". Then follow the steps to choose a database as above.

Note that after you have downloaded a CodeQL database once, you won't have to do that again, but if there are changes to the code upstream, you'll need to download the newest version of a CodeQL database.
