import setuptools

# Update this for every PyPi release
version = "6.4.1"

with open("README.md", "r") as fh:
    long_description = fh.read()


init_lines = ["# Created by setup.py, whilst creating a new PyPi release\n","__version__ = \"{}\"\n".format(version)]
with open("PDielec/__init__.py", "w") as fh:
    fh.writelines(init_lines)

setuptools.setup(
    name="PDielec",
    version=version,
    author="John Kendrick and Andrew Burnett",
    author_email="john@kendrick.co.uk",
    description="PDielec package for the calculation of THz and infrared spectra",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JohnKendrick/PDielec",
    # packages=setuptools.find_packages(include=['PDielec','PDielec.*']),
    packages=['PDielec','PDielec/GUI'],
    license='MIT',
    keywords=['Infrared',' Infrared Spectroscopy','THz Spectroscopy','Terahertz Spectroscopy','Bruggeman','Maxwell-Garnett','Effective Medium'],
    install_requires=[
        'numpy',
        'scipy',
        'PyQt5',
        'PyOpengl',
        'PyMieScatt',
        'matplotlib',
        'xlsxwriter',
        'openpyxl',
        'imageio',
        'imageio-ffmpeg',
        'PyYAML',
        'termcolor',
        'bleach',
        'docutils',
        'Pygments',
        'beautifulsoup4',
        'psutil',
        'requests',
    ],
    include_package_data=True,   # Assume all information is in MANIFEST.in
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=2.7',
    entry_points= {
        'console_scripts': [
            'preader=preader:main',
        ],
        'gui_scripts': [
            'pdgui=pdgui:main',
        ],
    },
)
