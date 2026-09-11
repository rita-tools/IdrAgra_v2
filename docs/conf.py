project = "IdrAgra"
author = "IdrAgra contributors"

extensions = [
    "myst_parser",
    "sphinx.ext.mathjax",
]

source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

root_doc = "pages/index"

html_theme = "furo"
html_static_path = ["_static"]
html_css_files = ["parameter-reference.css"]

myst_enable_extensions = [
    "attrs_block",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "substitution",
]

nitpicky = True
