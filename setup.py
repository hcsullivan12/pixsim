from setuptools import setup, find_packages
setup(
    name = 'pixsim',
    version = '0.0',
    package_dir = {'':'python'},
    packages = ['pixsim'],
    install_requires = [
        'Click',
        'meshio',
        'SQLAlchemy',
        'networkx',
    ],
    entry_points = dict(
        console_scripts = [
            'pixsim = pixsim.cli:main',
        ]
    )
)
