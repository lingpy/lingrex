from setuptools import setup, find_packages


setup(
        name='lingrex',
        description="Linguistic reconstruction with LingPy",
        long_description=open('README.md').read(),
        long_description_content_type="text/markdown",
        author='Johann-Mattis List and Robert Forkel',
        url='https://github.com/lingpy/lingrex',
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
        ],
        version='1.4.0',
        packages=find_packages(where='src'),
        package_dir={'': 'src'},
        install_requires=['lingpy>=2.6.10'],
        extras_require={
            "test": ["pytest", "coverage", "pytest-mock", "pytest-cov"],
            "dev": ["wheel", "twine", "tox", "black", "flake8"],
        },
        keywords="historical linguistics, computational linguistics, computer-assisted language comparison"
        )
