# from setuptools import setup, find_packages
# from setuptools.command.build_ext import build_ext
# import subprocess
# import glob
# import os
# import sys


# class CustomBuild(build_ext):
#     def run(self):

#         # g++ commands in the makefile
#         commands = [

#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/point.o', 'finalSrc/point.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/polyHelix.o', 'finalSrc/polyHelix.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/randomMolGen.o', 'finalSrc/randomMolGen.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/ktlMoleculeRandom.o', 'finalSrc/ktlMoleculeRandom.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/experimentalData.o', 'finalSrc/experimentalData.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/hydrationShellRandom.o', 'finalSrc/hydrationShellRandom.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/skmt.o', 'finalSrc/skmt.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/writheFP.o', 'finalSrc/writheFP.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/moleculeFitAndState.o', 'finalSrc/moleculeFitAndState.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/mainPredictionFinalQvar.o', 'finalSrc/mainPredictionFinalQvar.cpp'],

#             ['g++', '-O3', '-std=gnu++14', '-o', 'predictStructureQvary', 'finalSrc/point.o', 'finalSrc/polyHelix.o', 'finalSrc/randomMolGen.o', 'finalSrc/ktlMoleculeRandom.o', 'finalSrc/experimentalData.o', 'finalSrc/hydrationShellRandom.o', 'finalSrc/skmt.o', 'finalSrc/writheFP.o', 'finalSrc/moleculeFitAndState.o', 'finalSrc/mainPredictionFinalQvar.o'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/Flexible_generator.o', 'finalSrc/Flexible_generator.cpp'],
#             ['g++', '-O3', '-std=gnu++14', '-o', 'generate_structure', 'finalSrc/point.o', 'finalSrc/polyHelix.o', 'finalSrc/randomMolGen.o', 'finalSrc/ktlMoleculeRandom.o', 'finalSrc/Flexible_generator.o'] 
            
#             ]
        
#         try:
#             subprocess.check_call(['g++', '--version'])
#         except subprocess.CalledProcessError:
#             print("g++ is required to install this package. Please install g++ and try again.")
#             sys.exit(1)

#         for obj_file in glob.glob('finalSrc/*.o'):
#             os.remove(obj_file)

#         # Execute each g++ command
#         for cmd in commands:
#             try:
#                 subprocess.check_call(cmd)
#             except subprocess.CalledProcessError as e:
#                 print(f"Failed to execute {' '.join(cmd)}: {e}")
#                 sys.exit(1)

#         super().run()

# setup(
#     name='carbonarapy',
#     version='0.1',
#     packages=find_packages(),
#     cmdclass={ 'build_ext': CustomBuild },
#     scripts=['carbonarapy/carbonara_wrapper.py'],  # Path to your Python wrapper script
# )


from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
import subprocess
import glob
import os
import sys

# class CustomBuild(build_ext):
#     def run(self):
#         # Directory where the compiled binaries will be placed
#         build_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'build')
#         if not os.path.exists(build_dir):
#             os.makedirs(build_dir)

#         # g++ commands in the makefile
#         commands = [

#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/point.o', 'finalSrc/point.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/polyHelix.o', 'finalSrc/polyHelix.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/randomMolGen.o', 'finalSrc/randomMolGen.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/ktlMoleculeRandom.o', 'finalSrc/ktlMoleculeRandom.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/experimentalData.o', 'finalSrc/experimentalData.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/hydrationShellRandom.o', 'finalSrc/hydrationShellRandom.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/skmt.o', 'finalSrc/skmt.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/writheFP.o', 'finalSrc/writheFP.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/moleculeFitAndState.o', 'finalSrc/moleculeFitAndState.cpp'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/mainPredictionFinalQvar.o', 'finalSrc/mainPredictionFinalQvar.cpp'],

#             ['g++', '-O3', '-std=gnu++14', '-o', 'predictStructureQvary', 'finalSrc/point.o', 'finalSrc/polyHelix.o', 'finalSrc/randomMolGen.o', 'finalSrc/ktlMoleculeRandom.o', 'finalSrc/experimentalData.o', 'finalSrc/hydrationShellRandom.o', 'finalSrc/skmt.o', 'finalSrc/writheFP.o', 'finalSrc/moleculeFitAndState.o', 'finalSrc/mainPredictionFinalQvar.o'],
#             ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/Flexible_generator.o', 'finalSrc/Flexible_generator.cpp'],
#             ['g++', '-O3', '-std=gnu++14', '-o', 'generate_structure', 'finalSrc/point.o', 'finalSrc/polyHelix.o', 'finalSrc/randomMolGen.o', 'finalSrc/ktlMoleculeRandom.o', 'finalSrc/Flexible_generator.o'] 
            
#             ]

#         # Check if g++ is available
#         try:
#             subprocess.check_call(['g++', '--version'])
#         except subprocess.CalledProcessError:
#             print("g++ is required to install this package. Please install g++ and try again.")
#             sys.exit(1)

#         # Execute each g++ command
#         for cmd in commands:
#             print(f"Executing: {' '.join(cmd)}")
#             try:
#                 result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
#                 print(f"Output: {result.stdout}")
#                 if result.stderr:
#                     print(f"Errors: {result.stderr}", file=sys.stderr)
#             except subprocess.CalledProcessError as e:
#                 print(f"Failed to execute {' '.join(cmd)}: {e}", file=sys.stderr)
#                 sys.exit(1)

#         # After successful build, you might want to copy or move the binary to the package directory
#         # shutil.copy(os.path.join(build_dir, 'binary_name'), 'carbonarapy/')

#         super().run()

# setup(
#     name='carbonarapy',
#     version='0.1',
#     packages=find_packages(),
#     cmdclass={'build_ext': CustomBuild},
#     scripts=['carbonarapy/carbonara_wrapper.py'],  # Adjust if necessary
# )


from setuptools import setup, find_packages
from setuptools.command.install import install
import subprocess
import os

class CustomInstall(install):
    def run(self):
        # Compile the C++ code
        self.compile_cpp_code()
        # Run the standard install
        super().run()

    def compile_cpp_code(self):
        # Add your compilation commands here
        commands = [

            ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/point.o', 'finalSrc/point.cpp'],
            ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/polyHelix.o', 'finalSrc/polyHelix.cpp'],
            ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/randomMolGen.o', 'finalSrc/randomMolGen.cpp'],
            ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/ktlMoleculeRandom.o', 'finalSrc/ktlMoleculeRandom.cpp'],
            ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/experimentalData.o', 'finalSrc/experimentalData.cpp'],
            ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/hydrationShellRandom.o', 'finalSrc/hydrationShellRandom.cpp'],
            ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/skmt.o', 'finalSrc/skmt.cpp'],
            ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/writheFP.o', 'finalSrc/writheFP.cpp'],
            ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/moleculeFitAndState.o', 'finalSrc/moleculeFitAndState.cpp'],
            ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/mainPredictionFinalQvar.o', 'finalSrc/mainPredictionFinalQvar.cpp'],

            ['g++', '-O3', '-std=gnu++14', '-o', 'predictStructureQvary', 'finalSrc/point.o', 'finalSrc/polyHelix.o', 'finalSrc/randomMolGen.o', 'finalSrc/ktlMoleculeRandom.o', 'finalSrc/experimentalData.o', 'finalSrc/hydrationShellRandom.o', 'finalSrc/skmt.o', 'finalSrc/writheFP.o', 'finalSrc/moleculeFitAndState.o', 'finalSrc/mainPredictionFinalQvar.o'],
            ['g++', '-c', '-O3', '-std=gnu++14', '-o', 'finalSrc/Flexible_generator.o', 'finalSrc/Flexible_generator.cpp'],
            ['g++', '-O3', '-std=gnu++14', '-o', 'generate_structure', 'finalSrc/point.o', 'finalSrc/polyHelix.o', 'finalSrc/randomMolGen.o', 'finalSrc/ktlMoleculeRandom.o', 'finalSrc/Flexible_generator.o'] 
            
        ]
        for cmd in commands:
            print(f"Executing: {' '.join(cmd)}")
            subprocess.check_call(cmd)

        # Add any additional steps to copy/move the executable to the correct location

setup(
    name='carbonarapy',
    version='0.1',
    packages=find_packages(),
    cmdclass={
        'install': CustomInstall,
    },
    scripts=['carbonarapy/carbonara_wrapper.py'],
)
