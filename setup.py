import subprocess
from distutils.command.build import build as _build
import setuptools

with open("README.md", "r") as fh:
	long_description = fh.read()

# This class handles the pip install mechanism.
class build(_build):  # pylint: disable=invalid-name
  sub_commands = _build.sub_commands + [("CustomCommands", None)]

CUSTOM_COMMANDS = [
	["libdir=`ls -1 build|grep \"lib\"`; cd build/$libdir/silhouetteRank/ && Rscript --version"]]

class CustomCommands(setuptools.Command):
  """A setuptools Command class able to run arbitrary commands."""

  def initialize_options(self):
    pass

  def finalize_options(self):
    pass

  def RunCustomCommand(self, command_list):
    print("Running command: %s" % command_list)
    p = subprocess.Popen(
        command_list,
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    # Can use communicate(input='y\n'.encode()) if the command run requires
    # some confirmation.
    stdout_data, _ = p.communicate()
    print("Command output: %s" % stdout_data)
    if p.returncode != 0:
      raise RuntimeError(
          "Command %s failed: exit code: %s" % (command_list, p.returncode))

  def run(self):
    for command in CUSTOM_COMMANDS:
      self.RunCustomCommand(command)



setuptools.setup(
	name="silhouetteRank",
	version="1.0.5.13",
	author="Qian Zhu",
	author_email="zqian@jimmy.harvard.edu",
	description="silhouetteRank is a tool for finding spatially variable genes based on computing silhouette coefficient from binarized spatial gene expression data",
	long_description="",
	long_description_content_type="text/markdown",
	url="https://bitbucket.org/qzhu/silhouetteRank",
	packages=setuptools.find_packages(),
	entry_points = {
		"console_scripts": [
			"silhouette_rank_one = silhouetteRank.silhouette_rank_one:main",
			"silhouette_rank_main = silhouetteRank.evaluate_2b:main",
			"silhouette_rank_random = silhouetteRank.evaluate_exact_one_2b:main",
			"silhouette_rank_random_batch = silhouetteRank.evaluate_exact_2b:main",
			"slrank_fast = silhouetteRank.evaluate_fast_2b:main",
			"slrank_random_fast = silhouetteRank.evaluate_fast_one_2b:main",
			"slrank_prep_fast = silhouetteRank.prep_fast:main",
			"slrank_combine_fast = silhouetteRank.combine_fast:main",
		]
	},
	classifiers=(
		"Programming Language :: Python :: 3",
		"License :: OSI Approved :: MIT License",
		"Operating System :: OS Independent",
	),
	python_requires=">=3.5",
	package_data={"silhouetteRank":  ["do_gpd.R", "do_kmeans.R",
		"qval.R"]},
	install_requires=[
		"scipy", "numpy", "pandas", "seaborn", "scikit-learn", "matplotlib"],
	cmdclass={
		"build": build,
		"CustomCommands": CustomCommands,
		}	
)
	
