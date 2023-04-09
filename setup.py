from setuptools import find_packages, setup

setup(
    name="PDBSequenceAlignment",
    version="1.0.0",
    author="Jack",
    packages=find_packages(),
    test_suite="test",
    setup_requires=[
        "pytest-runner",
        "opencv-python",
    ],
    tests_require=[
        "pytest",
        "pytest-timeout",
    ],
)
