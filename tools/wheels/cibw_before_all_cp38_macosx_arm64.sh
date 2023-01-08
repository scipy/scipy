# https://cibuildwheel.readthedocs.io/en/stable/faq/#macos-building-cpython-38-wheels-on-arm64
curl -o /tmp/Python38.pkg https://www.python.org/ftp/python/3.8.10/python-3.8.10-macos11.pkg
sudo installer -pkg /tmp/Python38.pkg -target /
sh "/Applications/Python 3.8/Install Certificates.command"
