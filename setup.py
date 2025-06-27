from setuptools import setup, find_packages

setup(
    name='defect_tool',
    version='0.1',
    packages=find_packages(),  # This will automatically find 'defect' and 'scripts'
    install_requires=[
        'click',
        'MDAnalysis',
        'matplotlib',
        'numpy'
    ],
    entry_points={
        'console_scripts': [
            'defect-tool = scripts.cli:cli'
        ]
    },
)
