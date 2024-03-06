import setuptools

setuptools.setup(
    name="gau_pes",
    version="0.1.0",
    author="mizu-bai",
    author_email="shiragawa4519@outlook.com",
    description="Interface for invoking PES from Gaussian",
    url="https://github.com/CQPES/Gaussian-PES",
    packages=setuptools.find_packages(),
    classifiers=[],
    python_requires=">=3.8",
    install_requires=[
        "numpy>1.20.0",
    ],
)
