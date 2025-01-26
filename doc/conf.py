from flatsurf.version import version
from sage_docbuild.conf import html_theme, html_theme_options, pygments_style, pygments_dark_style, html_css_files, skip_TESTS_block, mathjax3_config, default_role

# -- General configuration ------------------------------------------------
extensions = [
    # We need to use SageMath's autodoc to render nested classes in categories
    # correctly. Otherwise they just render as "alias for" in the
    # documentation.
    "sage_docbuild.ext.sage_autodoc",
    "sphinx.ext.todo",
    "sphinx.ext.mathjax",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "myst_nb",
    "jupyter_sphinx",
]

# Extensions when rendering .ipynb/.md notebooks
myst_enable_extensions = [
    "dollarmath",
    "amsmath",
]

# The suffix of source filenames.
source_suffix = ".rst"

# The master toctree document.
master_doc = "index"

# General information about the project.
project = "sage-flatsurf"
copyright = "2016-2025, the sage-flatsurf authors"

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.
release = version

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ["_build", "news"]

# Allow linking to external projects, e.g., SageMath
intersphinx_mapping = {"sage": ("https://doc.sagemath.org/html/en/reference", None)}

# -- Options for HTML output ----------------------------------------------
if html_css_files != ["custom-furo.css", "custom-jupyter-sphinx.css", "custom-codemirror-monokai.css"]:
    raise NotImplementedError(
        "CSS customization has changed in SageMath. The configuration of sage-flatsurf documentation build needs to be updated."
    )

html_css_files = [
    "https://doc.sagemath.org/html/en/reference/_static/custom-furo.css",
    "https://doc.sagemath.org/html/en/reference/_static/custom-jupyter-sphinx.css",
    "https://doc.sagemath.org/html/en/reference/_static/custom-codemirror-monokai.css",
    "jupyter_execute.css"
]

html_theme_options["light_logo"] = html_theme_options["dark_logo"] = "logo.svg"
html_static_path = ["static"]

# Output file base name for HTML help builder.
htmlhelp_basename = "sage-flatsurfdoc"

# Options for jupyter-sphinx
jupyter_execute_default_kernel = "sagemath"

# Only rerender example notebooks when the cache is stale.
nb_execution_mode = "cache"

linkcheck_ignore = [
    # Zulip channel requires signing up first.
    r'https://sagemath.zulipchat.com/#narrow/stream/271193-polygon/topic/hyperbolic.20geometry/near/284722650',
    # PyPI seems to have temporarily vanished from repology
    'https://repology.org/project/python:sage-flatsurf/packages',
]

def setup(app):
    app.connect('autodoc-process-docstring', skip_TESTS_block)
