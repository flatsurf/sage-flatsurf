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
copyright = \"2016-2023, the sage-flatsurf authors\"

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

if html_css_files != ["custom-furo.css"]:
    raise NotImplementedError(
        "CSS customization has changed in SageMath. The configuration of sage-flatsurf documentation build needs to be updated."
    )

html_css_files = ["https://doc.sagemath.org/html/en/reference/_static/custom-furo.css"]

html_theme_options["light_logo"] = html_theme_options["dark_logo"] = "logo.svg"
html_static_path = ["static"]

# Output file base name for HTML help builder.
htmlhelp_basename = "sage-flatsurfdoc"

# -- Options for LaTeX output ---------------------------------------------

latex_elements = {}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (
        "index",
        "sage-flatsurf.tex",
        "sage-flatsurf Documentation",
        "the sage-flatsurf authors",
        "manual",
    ),
]

# -- Options for manual page output ---------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    (
        "index",
        "sage-flatsurf",
        "sage-flatsurf Documentation",
        ["the sage-flatsurf authors"],
        1,
    )
]

# -- Options for Texinfo output -------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    (
        "index",
        "sage-flatsurf",
        "sage-flatsurf Documentation",
        "the sage-flatsurf authors",
        "sage-flatsurf",
        "One line description of project.",
        "Miscellaneous",
    ),
]
