# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2021-2024 Julian Rüth
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ********************************************************************

import os.path

from rever.activities.command import command

try:
  input("Are you sure you are on the master branch which is identical to origin/master? [ENTER]")
except KeyboardInterrupt:
  sys.exit(1)

$PROJECT = 'sage-flatsurf'

command('pixi', 'pixi install --manifest-path "$PWD/pyproject.toml" -e dev')

command('build', 'python -m build')
command('twine', 'twine upload dist/sage_flatsurf-' + $VERSION + '.tar.gz dist/sage_flatsurf-' + $VERSION + '-py3-none-any.whl')

$ACTIVITIES = [
    'version_bump',
    'pixi',
    'changelog',
    'tag',
    'push_tag',
    'build',
    'twine',
    'ghrelease',
]

$RELEASE_YEAR = $RELEASE_DATE.year

$VERSION_BUMP_PATTERNS = [
    ('flatsurf/version.py', r"version =", "version = \"$VERSION\""),
    ('README.md', r"tar zxf sage-flatsurf-.*.unix.tar.gz", "tar zxf sage-flatsurf-$VERSION.unix.tar.gz"),
    ('README.md', r"./sage-flatsurf-.*/jupyterlab  # or", "./sage-flatsurf-$VERSION/jupyterlab  # or"),
    ('README.md', r"./sage-flatsurf-.*/sage", "./sage-flatsurf-$VERSION/sage"),
    ('doc/index.rst', r' :target: https://mybinder.org/v2/gh/flatsurf/sage-flatsurf', r' :target: https://mybinder.org/v2/gh/flatsurf/sage-flatsurf/$VERSION?filepath=doc%2Fexamples'),
    ('doc/conf.py', r'copyright = ', "copyright = \"2016-$RELEASE_YEAR, the sage-flatsurf authors\""),
    ('pyproject.toml', r'version = ', 'version = "$VERSION"'),
    ('doc/install.rst', r"  curl -fsSL https://github.com/flatsurf/sage-flatsurf/releases/download/", r"  curl -fsSL https://github.com/flatsurf/sage-flatsurf/releases/download/$VERSION/sage-flatsurf-$VERSION.unix.tar.gz | tar zxf -"),
    ('doc/install.rst', r"  ./sage-flatsurf-.*/sage", r"  ./sage-flatsurf-$VERSION/sage"),
    ('doc/install.rst', r"  ./sage-flatsurf-.*/jupyterlab", r"  ./sage-flatsurf-$VERSION/jupyterlab"),
    ('installer/win/installer.iss', r"AppCopyright=", "AppCopyright=Copyright (C) 2016-$RELEASE_YEAR the sage-flatsurf authors"),
]

$CHANGELOG_FILENAME = 'ChangeLog'
$CHANGELOG_TEMPLATE = 'TEMPLATE.rst'
$CHANGELOG_NEWS = 'doc/news'
$PUSH_TAG_REMOTE = 'git@github.com:flatsurf/sage-flatsurf.git'

$PYPI_BUILD_COMMANDS = ['sdist', 'bdist_wheel']
$PYPI_NAME = "sage-flatsurf"

$GITHUB_ORG = 'flatsurf'
$GITHUB_REPO = 'sage-flatsurf'

$CHANGELOG_CATEGORIES = ('Added', 'Changed', 'Deprecated', 'Removed', 'Fixed', 'Performance')
