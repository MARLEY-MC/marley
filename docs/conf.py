import sys, os, subprocess

from sphinx.highlighting import lexers
from pygments.lexers.web import PhpLexer

project = u'MARLEY'
copyright = u'2016-2026 Steven Gardiner'
master_doc = 'index'
templates_path = [ '_templates' ]
extensions = [ 'sphinxcontrib.bibtex', 'sphinxcontrib.newsfeed',
  'sphinx.ext.todo', 'sphinxcontrib.katex' ]

# KaTeX configuration: render math at build time (no JavaScript needed at
# runtime)
katex_prerender = True
source_suffix = '.rst'
version = '2.0.0'
exclude_patterns = ['_build']

highlight_language = 'none'

# TODO notes

#todo_include_todos = True
todo_include_todos = False

# -- HTML theme settings ------------------------------------------------

html_favicon = 'mar.png'
html_show_sourcelink = False
html_sidebars = {
    '**': ['logo-text.html',
           'globaltoc.html',
           'localtoc.html',
           'searchbox.html']
}

import guzzle_sphinx_theme

extensions.append("guzzle_sphinx_theme")
html_theme_path = guzzle_sphinx_theme.html_theme_path()
html_theme = 'guzzle_sphinx_theme'

# Guzzle theme options (see theme.conf for more information)
html_theme_options = {
    "base_url": "www.marleygen.org",
}

#html_add_permalinks = None
html_css_files = [ 'custom.css' ]
html_static_path = [ '_static' ]
html_extra_path = [ '../CITATION.bib', './CNAME' ]

# Bibliography files used by the sphinxcontrib.bibtex extension
bibtex_bibfiles = [ 'marley_pubs.bib', 'external_pubs.bib' ]


def setup(app):
    import re
    from docutils.nodes import FixedTextElement, Text, literal, math, Element
    from sphinx.transforms.post_transforms import SphinxPostTransform

    class MathDollarPostTransform(SphinxPostTransform):
        default_priority = 11

        _re_math = re.compile(r'\$([^$]+)\$')

        def run(self, **kwargs):
            for parent in list(self.document.traverse(Element)):
                if isinstance(parent, (FixedTextElement, literal, math)):
                    continue
                p = parent.parent
                skip = False
                while p:
                    if isinstance(p, (FixedTextElement, literal, math)):
                        skip = True
                        break
                    p = p.parent
                if skip:
                    continue

                children = list(parent.children)
                i = 0
                while i < len(children):
                    if not isinstance(children[i], Text):
                        i += 1
                        continue
                    j = i
                    while j < len(children) and isinstance(children[j], Text):
                        j += 1
                    merged = ''.join(str(c) for c in children[i:j])
                    parts = []
                    last_end = 0
                    has_math = False
                    for m in self._re_math.finditer(merged):
                        if m.start() > last_end:
                            parts.append(Text(merged[last_end:m.start()]))
                        parts.append(math(m.group(1), Text(m.group(1))))
                        has_math = True
                        last_end = m.end()
                    if has_math:
                        if last_end < len(merged):
                            parts.append(Text(merged[last_end:]))
                        for idx in range(j - 1, i - 1, -1):
                            parent.remove(children[idx])
                        for idx, p in enumerate(parts):
                            parent.insert(i + idx, p)
                    i = j

    app.add_post_transform(MathDollarPostTransform)
