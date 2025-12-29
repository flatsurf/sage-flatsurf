# Configuration file for the Sphinx documentation builder.

import os
from pathlib import Path
from sage_docbuild.conf import skip_TESTS_block

BUILDDIR = Path(os.environ.get('ABS_BUILDDIR', '.')).absolute()

project = "sage-flatsurf"
copyright = "2016-2025, the sage-flatsurf authors"
author = 'the sage-flatsurf authors'

release = '0.8.0'

extensions = [
    # We need to use SageMath's autodoc to render nested classes in categories
    # correctly. Otherwise they just render as "alias for" in the
    # documentation.
    "sage_docbuild.ext.sage_autodoc",
    'sphinx.ext.intersphinx',
    "sphinx.ext.todo",
    'sphinx.ext.mathjax',
    "sphinx.ext.viewcode",
    'sphinx_book_theme',
    "myst_nb",
    "jupyter_sphinx",
]

# Extensions when rendering .ipynb/.md notebooks
myst_enable_extensions = [
    "dollarmath",
    "amsmath",
]

templates_path = ['_templates']

exclude_patterns = ["_build", "news"]

html_theme = 'sphinx_book_theme'

html_logo = 'static/logo.svg'
html_static_path = ["static"]

html_theme_options = {
    "repository_url": "https://github.com/flatsurf/sage-flatsurf",
    "icon_links": [{
        "name": "flatsurf",
        "url": "https://flatsurf.github.io",
        "icon": "https://flatsurf.github.io/assets/logo.svg",
        "type": "url",
    }, {
        "name": "GitHub",
        "url": "https://github.com/flatsurf/sage-flatsurf",
        "icon": "fa-brands fa-square-github",
        "type": "fontawesome",
    }, {
        "name": "Zulip",
        "url": "https://sagemath.zulipchat.com/#narrow/channel/271193-flatsurf",
        "icon": "fa-regular fa-comments",
        "type": "fontawesome",
    },
    ],
    "use_edit_page_button": True,
    "repository_branch": "master",
    "path_to_docs": "doc",
}

html_static_path = ['static']

html_css_files = ['extra.css']

intersphinx_mapping = {"sage": ("https://doc.sagemath.org/html/en/reference", None)}

html_css_files = [
    "jupyter_execute.css",
]

htmlhelp_basename = "sage-flatsurfdoc"

jupyter_execute_default_kernel = "sagemath"

nb_execution_mode = "cache"

linkcheck_ignore = [
    # Zulip channel requires signing up first.
    r'https://sagemath.zulipchat.com/#narrow/stream/271193-polygon/topic/hyperbolic.20geometry/near/284722650',
    # PyPI seems to have temporarily vanished from repology
    'https://repology.org/project/python:sage-flatsurf/packages',
    # ACM blocks GitHub runners
    'http://dx.doi.org/10.1145/2755996.2756664',
]

def setup(app):
    app.connect('autodoc-process-docstring', skip_TESTS_block)
