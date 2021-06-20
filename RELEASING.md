
Releasing linse
===============

- Do platform test via tox:
  ```shell script
  tox -r
  ```
- test if the workflow scripts still work:
  ```
  cd tests/workflows/list-2019
  make
  cd ../tests/workflows/bodt-2019
  make
  cd ../tests/workflows/wu-2020
  make
  ```

- Make sure statement coverage >= 99%
- Use black to make the code unified:
  ```
  black src/lingrex/*.py
  ```

- Update the version number, by removing the trailing `.dev0` in:
  - `setup.py`
  - `src/lingrex/__init__.py`

- Create the release commit:
  ```shell script
  git commit -a -m "release <VERSION>"
  ```

- Create a release tag:
  ```shell script
  git tag -a v<VERSION> -m"<VERSION> release"
  ```

- Release to PyPI (see https://github.com/di/markdown-description-example/issues/1#issuecomment-374474296):
  ```shell script
  rm dist/*
  python setup.py sdist
  twine upload dist/*
  rm dist/*
  python setup.py bdist_wheel
  twine upload dist/*
  ```

- Push to github:
  ```shell script
  git push origin
  git push --tags
  ```

- Change version for the next release cycle, i.e. incrementing and adding .dev0
  - `setup.py`
  - `src/lingrex/__init__.py`

- Commit/push the version change:
  ```shell script
  git commit -a -m "bump version for development"
  git push origin
  ```
