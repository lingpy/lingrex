try:
    from setuptools import setup, find_packages
except ImportError:
    from distribute_setup import use_setuptools
    use_setuptools()
    from setuptools import setup, find_packages

setup(
        name='lingrex',
        description="Linguistic reconstruction with LingPy",
        author='Johann-Mattis List',
        url='https://github.com/lingpy/lingrex',
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
        ],
        version='0.1.2',
        packages=find_packages(where='src'),
        package_dir={'': 'src'},
        install_requires=['lingpy', 'sinopy'],
        keywords="historical linguistics, computational linguistics, computer-assisted language comparison"
        )
