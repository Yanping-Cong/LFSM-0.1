from setuptools import setup, find_packages
import os
def readme():
    absolute_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),'README.rst')
    with open(absolute_path) as f:
        return f.read()
setup(name='LFSM',
      version='0.1',
      description='The low-frequency sky model',
      long_description=readme(),
      classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Text Processing :: Linguistic',
      ],
      keywords='the LFSM is a temporary name',
      url='http://github.com/storborg/funniest',
      author='Yanping Cong',
      author_email='ypcong@bao.ac.cn',
      license='MIT',
      #packages=['funniest/funniest'],
      packages = find_packages(),
      #install_requires=[
      #    'markdown',
      #],
      dependency_links = [
          'https://github.com/zuoshifan/caput/tree/zuo/develop'
          ],
      include_package_data=True,
      zip_safe=False)
