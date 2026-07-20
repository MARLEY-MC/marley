import sys, os, subprocess, re

from sphinx.highlighting import lexers
from pygments.lexers.web import PhpLexer
from sphinx_math_dollar import split_dollars

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
bibtex_bibfiles = [ 'marley_pubs.bib', 'external_pubs.bib', 'external_pubs_refs.bib' ]


def setup(app):
    from docutils.nodes import FixedTextElement, Text, literal, math, Element
    from sphinx.transforms.post_transforms import SphinxPostTransform

    # ── Phase 1: protect $...$ in parsed bib entries before pybtex formats them ──

    app._math_placeholder_map = {}

    def _protect_bib_math_after_parse(app):
        dom = app.env.get_domain('cite')
        if dom is None:
            return
        bibdata = dom.data.get('bibdata')
        if bibdata is None:
            return

        counter = [0]
        math_map = {}

        def _hash_math(m):
            counter[0] += 1
            ph = '@@' + str(counter[0]) + '@@'  # digits-only, survives case changes
            math_map[ph] = m.group(0)
            return ph

        for entry in bibdata.data.entries.values():
            for field_name in list(entry.fields.keys()):
                value = entry.fields[field_name]
                if not isinstance(value, str):
                    continue
                protected = re.sub(r'\$\$[^$]+\$\$', _hash_math, value)
                protected = re.sub(r'\$[^$]+\$', _hash_math, protected)
                if protected != value:
                    entry.fields[field_name] = protected

        app._math_placeholder_map = math_map

    app.connect('builder-inited', _protect_bib_math_after_parse, priority=501)

    # ── Phase 2: restore placeholders and create math nodes ──

    class MathDollarPostTransform(SphinxPostTransform):
        default_priority = 11

        def run(self, **kwargs):
            math_map = getattr(
                self.document.settings.env.app, '_math_placeholder_map', {})

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

                    for ph, orig in math_map.items():
                        merged = merged.replace(ph, orig)

                    fragments = split_dollars(merged)
                    has_math = any(t != 'text' for t, _ in fragments)
                    if has_math:
                        parts = []
                        for typ, content in fragments:
                            if typ == 'text':
                                if content:
                                    parts.append(Text(content))
                            else:
                                parts.append(math(content, Text(content)))
                        if parts:
                            for idx in range(j - 1, i - 1, -1):
                                parent.remove(children[idx])
                            for idx, p in enumerate(parts):
                                parent.insert(i + idx, p)
                    i = j

    app.add_post_transform(MathDollarPostTransform)
