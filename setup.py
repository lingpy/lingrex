import distribute_setup
distribute_setup.use_setuptools()

try:
    from setuptools import setup, find_packages, Extension
except ImportError:
    from distribute_setup import use_setuptools
    use_setuptools()
    from setuptools import setup, find_packages, Extension

# setup package name etc as a default
pkgname = 'lingrex'
pkg_dir = {'':'.'}
pkg_location = '.'

setup(
        name=pkgname,
        version='0.1',
        packages=find_packages(pkg_location),
        package_dir=pkg_dir,
        install_requires=['lingpy', 'clldutils'],
        author='Johann-Mattis List'
        )
