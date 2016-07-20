from setuptools import setup, find_packages

setup(
    name = "SMolPhot",
    version = "0.1",
    url = "https://bitbucket.org/ardiloot/smolphot-software",
    packages = find_packages(),
    install_requires = 
        ["numpy",
        "scipy",
        "pyyaml",
        "pillow",
        "scikit-image",
        "cython",
        "pyqt4",
        "pyqtgraph",
        "matplotlib",
        "appdirs"],
    entry_points = {
        "console_scripts": [],
        "gui_scripts": ["SMolPhot = SMolPhotGui.MainWindow:Run"]
        }
)
