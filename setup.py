import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pdielc-pkg-john_kendrick", # Replace with your own username
    version="6.3.0",
    author="John Kendrick and Andrew Burnett",
    author_email="john@kendrick.co.uk",
    description="PDielec package for the calculation of THz and infrared spectra",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/JohnKendrick/PDielec",
    packages=setuptools.find_packages(include=['PDielec','PDielec.*']),
    install_requires=[
        'numpy',
        'scipy',
        'PyYAML',
        'matplotlib',
        'termcolor',
        'xlsxwriter',
        'openpyxl',
        'imageio',
        'PyMieScatt',
    ],
    include_package_data=True,   # Assume all information is in MANIFEST.in
    classifiers=[
        "Development Status :: 5 - Production/stable",
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
    package_dir={'':'Python'}
    entry_points= {
        'console_scripts': [
            'preader=preader:main',
        ],
        'gui_scripts': [
            'pdgui=pdgui:main',
        ],
    },
)

