# SPDX-FileCopyrightText: 2023 The dune-iga developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

try:
    from dune.packagemetadata import metaData
except ImportError:
    from packagemetadata import metaData
from skbuild import setup

# When building a new package, update the version numbers below and run:
# access  docker container with mounted repo to /tmp/dune-iga
# build _iga
# cd /tmp/dune-iga
# /dune/dune-common/build-cmake/run-in-dune-env pip install twine scikit-build
# git config --global --add safe.directory /tmp/dune-iga
# /dune/dune-common/build-cmake/run-in-dune-env python setup.py sdist
# /dune/dune-common/build-cmake/run-in-dune-env python -m twine upload dist/* --verbose

duneigaVersion = "0.1.7"
duneVersion = "2.10.0"

metadata = metaData(duneVersion)[1]
metadata["version"] = duneigaVersion

metadata["name"] = "dune-iga-experimental"
setup(**metadata)
