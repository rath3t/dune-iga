import sys

# This file only exists to properly locally test the python package
# It takes care of adding the local package path to the path where python search for modules
# This is usually done in the configure step by Dune
# But e.g. in Clion and using a docker container ist runs every command in a new container
# Therefore in the configure step PYTHONPATH it set to contain the correct path
# But this is gone since Ctest is executed in a new container

# Every Python test needs to start with
# import setpath
#  setpath.set_path()
# Maybe there is a fix comiong from clion, see the bugreport https://youtrack.jetbrains.com/issue/CPP-32280

def set_path():
    sys.path.append('/tmp/dune-iga/build-cmake/python')