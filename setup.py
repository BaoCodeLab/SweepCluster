from setuptools import setup, find_packages

setup(
    name="SweepCluster",
    version='1.2',
    description="Selective sweep detection method based on SNP clustering.",
    classifiers=['License :: GPL-3.0','Operating System :: POSIX :: Linux'],
    author="Junhui Qiu and Yun-Juan Bao",
    url="https://github.com/BaoCodeLab/SweepCluster",
    packages=find_packages(),
    python_requires=">3.7.0",
    install_requires=["numpy >= 1.17","scipy >= 1.4","pandas >= 1.1.4","PyVCF","scikit-learn >= 0.22"]
)
