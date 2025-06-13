from setuptools import setup, find_packages

setup(
    name="dartsgpt",
    version="0.1.0",
    packages=find_packages(where=".", exclude=["tests*", "examples*"]),
    package_dir={"": "."},
    include_package_data=True,
)