from flatsurf.version import version
import sage_docbuild.conf

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
copyright = "2016-2024, the sage-flatsurf authors"

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

# Imitate the look of the SageMath documentation.
html_theme = sage_docbuild.conf.html_theme
html_theme_options = sage_docbuild.conf.html_theme_options
pygments_style = sage_docbuild.conf.pygments_style
pygments_dark_style = sage_docbuild.conf.pygments_dark_style
html_css_files = sage_docbuild.conf.html_css_files

if html_css_files != ["custom-furo.css", "custom-jupyter-sphinx.css", "custom-codemirror-monokai.css"]:
    raise NotImplementedError(
        "CSS customization has changed in SageMath. The configuration of sage-flatsurf documentation build needs to be updated."
    )

html_css_files = [
    "https://doc.sagemath.org/html/en/reference/_static/custom-furo.css",
    "https://doc.sagemath.org/html/en/reference/_static/custom-jupyter-sphinx.css",
    "https://doc.sagemath.org/html/en/reference/_static/custom-codemirror-monoai.css"
]

html_theme_options["light_logo"] = html_theme_options["dark_logo"] = "logo.svg"
html_static_path = ["static"]

# Output file base name for HTML help builder.
htmlhelp_basename = "sage-flatsurfdoc"


# Only rerender example notebooks when the cache is stale.
nb_execution_mode = "cache"

