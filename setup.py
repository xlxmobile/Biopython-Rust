from setuptools import setup, find_packages

setup(
    name="biopython-rust",  # 包名（确保 PyPI 上未被占用）
    version="0.1.0",         # 版本号
    description="High-performance Biopython implementation with Rust backend (WIP)",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Xuliang Xu",
    author_email="xlxmobile@outlook.com",
    url="https://github.com/xlxmobile/Biopython-Rust",
    license="MIT",
    packages=find_packages(),  # 自动找到 biopython_rust 包
    python_requires=">=3.7",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Development Status :: 3 - Alpha",  # 表示半成品
    ],
)