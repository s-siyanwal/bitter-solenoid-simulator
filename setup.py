from setuptools import setup

setup(
    name="bittersim",
    version="1.0.0",
    description="Real-time 0.5 T water-cooled copper Bitter solenoid simulator (pre-July-2018 scientific stack)",
    package_dir={"": "package"},
    packages=["bittersim"],
    python_requires=">=3.6",
    install_requires=["numpy>=1.14", "scipy>=1.1", "numba>=0.38", "matplotlib>=2.2"],
    entry_points={"console_scripts": ["bittersim=bittersim.cli:main"]},
)
