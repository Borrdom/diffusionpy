from setuptools import setup, find_packages

setup(
  name='PyFusion',
  license='MIT',
  version='0.0.1',
  description='Software to model diffusion in/out of amorphous substances',
  author='Dominik Borrmann',
  author_email='dominik.borrmann@tu-dortmund.de',
  url='https://github.com/Borrdom/PyFusion',
  download_url='https://github.com/Borrdom/PyFusion.git',
  long_description="Software to model diffusion in/out of amorphous substances",
  packages=find_packages(),
  install_requires=['numpy', 'scipy','casadi'],
  platforms=["Windows", "Linux", "Mac OS", "Unix"],
  keywords=['PC-SAFT'],
  zip_safe=False
)
