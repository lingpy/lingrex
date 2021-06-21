
# Releasing lingrex

- Do platform test via tox:
  ```shell
  tox -r
  ```
- test if the workflow scripts still work:
  ```shell
  pytest -m"workflow"
  ```

- Make sure statement coverage >= 99%
- Use black and flake8 to make the code unified:
  ```shell
  flake8 src
  black src/lingrex/*.py
  ```

- Update the version number, by removing the trailing `.dev0` in:
  - `setup.py`
  - `src/lingrex/__init__.py`

- Create the release commit:
  ```shell
  git commit -a -m "release <VERSION>"
  ```

- Create a release tag:
  ```shell
  git tag -a v<VERSION> -m"<VERSION> release"
  ```

- Release to PyPI:
  ```shell
  rm dist/*
  python setup.py sdist bdist_wheel
  twine upload dist/*
  ```

- Push to github:
  ```shell
  git push origin
  git push --tags
  ```

- Change version for the next release cycle, i.e. incrementing and adding .dev0
  - `setup.py`
  - `src/lingrex/__init__.py`

- Commit/push the version change:
  ```shell
  git commit -a -m "bump version for development"
  git push origin
  ```
