import setuptools

setuptools.setup(
    name="mcsimnp",
    version="0.0.1",
    author="Usman Jamshed",
    author_email="author@example.com",
    description="A small example package",
    long_description="A sample Python package which performs MC simulation.",
    long_description_content_type="text/markdown",
    url="https://github.com/USERNAME/REPONAME",
    project_urls={
        "Bug Tracker": "https://github.com/pypa/sampleproject/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    package_dir={"": "."},
    packages=setuptools.find_packages(where="."),
    python_requires=">=3.9",
)