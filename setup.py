from setuptools import setup, find_packages
with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='siccus',
    version='1.0.0',
    description='Scripts relevant to the DRYSAT project',
    author='Samuel Massart',
    author_email='samuel.massart@geo.tuwien.ac.at',
    packages=find_packages(),
    install_requires=requirements
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
)
