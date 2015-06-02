from setuptools import setup, find_packages

# def readme():
#    with open('README.md') as f:
#        return f.read()

with open("requirements.txt", "r") as f:
        install_requires = [x.strip() for x in f.readlines() if not x.startswith("#")]


setup(name='ecov',
      version='0.99.1',
      description='exome coverage calculation',
      # long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      keywords='custom pipelines',
      package_data={'ecov': ['*.Rmd']},
      # url='http://github.com/lpantano/ich-wrapper',
      author='Lorena Pantano',
      author_email='lpantano@iscb.org',
      license='MIT',
      packages=find_packages(),
      scripts=['scripts/coverage.py'],
      install_requires=install_requires,
      include_package_data=True,
      zip_safe=False)
