import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PDielec",
    version="6.3.1",
    author="John Kendrick and Andrew Burnett",
    author_email="john@kendrick.co.uk",
    description="PDielec package for the calculation of THz and infrared spectra",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JohnKendrick/PDielec",
    # packages=setuptools.find_packages(include=['PDielec','PDielec.*']),
    packages=['PDielec'],
    license='MIT',
    keywords=['Infrared',' Infrared Spectroscopy','THz Spectroscopy','Terahertz Spectroscopy','Bruggeman','Maxwell-Garnett','Effective Medium'],
    install_requires=[
        'numpy',
        'scipy',
        'PyQt5',
        'PyMieScatt',
        'matplotlib',
        'xlsxwriter',
        'openpyxl',
        'imageio',
        'PyYAML',
        'termcolor',
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

