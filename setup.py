from setuptools import setup

setup(name='scram2plot',
      version='0.1.0',
      description='scram2plot',
      author='Stephen Fletcher',
      author_email='s.fletcher@uq.edu.au',
      license='BSD3',
      packages=['scram2plot'],
      classifiers=[
    # How mature is this project? Common values are
    #   3 - Alpha
    #   4 - Beta
    #   5 - Production/Stable
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research',
    # Indicate who your project is intended for

    'Topic :: Scientific/Engineering :: Bio-Informatics',

    # Pick your license as you wish (should match "license" above)
     'License :: OSI Approved :: BSD3 License',

    # Specify the Python versions you support here. In particular, ensure
    # that you indicate whether you support Python 2, Python 3 or both.
    'Programming Language :: Python :: 3.11'],
      install_requires=['numpy','matplotlib','bokeh'],
      scripts=['scram_plot/scram2plot.py',
      'scram_plot/profile2plot.py',
      'scram_plot/compare2plot.py'
      ],
      zip_safe=False)