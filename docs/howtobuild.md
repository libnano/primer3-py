### Building for upload to PyPi
#### Dependencies:
1. setuptools
2. twine
3. wheel

#### edit the `~/.pypirc` file:

```ini
    [distutils]
    index-servers = pypi

    [pypi]
    repository = https://pypi.python.org/pypi
    username = <username>
```
Then execute:

```bash
$ python setup.py sdist
$ python setup.py bdist
$ python setup.py bdist_wheel
$ twine upload dist/*
```

Repeat for other versions w/ `bdist` and `twine upload dist/<fn of wheel>`
