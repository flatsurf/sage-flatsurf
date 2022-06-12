# ********************************************************************
#  This file is part of sage-flatsurf.
#
#        Copyright (C) 2021-2022 Julian Rüth
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

try:
  input("Are you sure you are on the master branch which is identical to origin/master? [ENTER]")
except KeyboardInterrupt:
  sys.exit(1)

$PROJECT = 'sage-flatsurf'

$ACTIVITIES = [
    'version_bump',
    'changelog',
    'tag',
    'push_tag',
    'pypi',
    'ghrelease',
]

$VERSION_BUMP_PATTERNS = [
    ('recipe/meta.yaml', r"\{% set version =", r"{% set version = '$VERSION' %}"),
    ('recipe/meta.yaml', r"\{% set build_number =", r"{% set build_number = '0' %}"),
    ('flatsurf/version.py', r"version =", r"version = '$VERSION'"),
    ('setup.py', r"    version=", r"    version='$VERSION',"),
    ('flatsurf.yml', r"  - sage-flatsurf=", r"  - sage-flatsurf=$VERSION"),
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
