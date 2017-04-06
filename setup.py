from distutils.core import setup

setup(
    name='rbp-maps',
    version='0.0.2',
    packages=['peak', 'density', 'analysis'],
    url='',
    license='',
    author='brianyee',
    author_email='',
    description='RNA-binding protein maps for region/splicing',
    package_dir={
        'peak': 'maps/peak', 'density': 'maps/density',
        'analysis': 'maps/analysis'
    }
)
