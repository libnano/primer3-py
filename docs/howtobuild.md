# Building for upload to PyPi
# Deps: setuptools, twine, wheel

# .pypirc file in ~ like: 

    [distutils]
    index-servers = pypi

    [pypi]
    repository = https://pypi.python.org/pypi
    username = <username>


1. python setup.py sdist
2. python setup.py bdist
3. twine upload dist/*

# repeat for other versions w/ bdist and `twine upload dist/<fn of wheel>
