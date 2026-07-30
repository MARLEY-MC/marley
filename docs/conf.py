import sys, os, subprocess, re

from sphinx.highlighting import lexers
from pygments.lexers.web import PhpLexer
from sphinx_math_dollar import split_dollars

# ---------------------------------------------------------------------------
# Custom pybtex style: APS/RevTeX-inspired formatting
#
# Implements:
#   - Newest-first sorting (by year, then arXiv submission month)
#   - Author names as "F. M. Lastname" with abbreviated first names
#   - 3-author cap before italic "et al."
#   - Collaboration field rendered as "(Name Collaboration)" after author list
#   - DOI embedded as hyperlink on the journal/volume/pages/year block
#   - arXiv eprint rendered as separate "arXiv:XXXX.XXXXX [cat]" link
#   - Graceful fallback to bare "(year) + arXiv link" when no journal/DOI
#   - Numeric labels ([1], [B1], [C1] with labelprefix)
# ---------------------------------------------------------------------------
import pybtex.plugin
from pybtex.style.formatting.unsrt import Style as UnsrtStyle
from pybtex.style.sorting import BaseSortingStyle
from pybtex.style.template import node, join, sentence
from pybtex.richtext import Text, String, Symbol, HRef, Tag as RichTag


# ---------------------------------------------------------------------------
# Sorting: newest-first by year, then by arXiv submission month
# ---------------------------------------------------------------------------

def _bib_sort_key(entry):
    """Return a sort key that orders entries newest-first.

    Primary key  : publication year (descending → negate).
    Secondary key: arXiv month encoded in the eprint ID, e.g. "2502.xxxxx"
                   gives month 2 of year 25 (descending → negate).
                   Falls back to the ``month`` field if no eprint is present.
    """
    fields = entry.fields

    # Year (integer, descending)
    try:
        year = -int(fields.get('year', 0))
    except (ValueError, TypeError):
        year = 0

    # Sub-year ordering from arXiv eprint ID (YYMM.NNNNN format)
    eprint = fields.get('eprint', '') or ''
    month_num = 0
    if eprint:
        parts = eprint.strip().split('.')
        if len(parts) >= 2 and len(parts[0]) == 4 and parts[0].isdigit():
            # New-style arXiv ID: first 4 digits are YYMM
            try:
                month_num = -int(parts[0][2:4])  # descending
            except ValueError:
                pass

    # Fallback: named month field → convert to integer
    if month_num == 0:
        month_str = str(fields.get('month', '') or '').strip().lower()
        month_names = {
            'jan': 1, 'feb': 2, 'mar': 3, 'apr': 4,
            'may': 5, 'jun': 6, 'jul': 7, 'aug': 8,
            'sep': 9, 'oct': 10, 'nov': 11, 'dec': 12,
        }
        for abbr, num in month_names.items():
            if month_str.startswith(abbr):
                month_num = -num
                break
        if month_num == 0:
            try:
                month_num = -int(month_str)
            except (ValueError, TypeError):
                pass

    return (year, month_num)


class _NewestFirstSortStyle(BaseSortingStyle):
    """Sort bibliography entries from newest to oldest."""

    default_suffixes = None

    def sorting_key(self, entry):
        return _bib_sort_key(entry)


pybtex.plugin.register_plugin(
    'pybtex.style.sorting', 'newest_first', _NewestFirstSortStyle)


# ---------------------------------------------------------------------------
# Formatting: APS/RevTeX-inspired style
# ---------------------------------------------------------------------------

# Maximum number of authors to list before switching to "et al."
_MAX_AUTHORS = 3


def _get_field(entry, name):
    """Return entry field value as a stripped string, or '' if absent/empty."""
    val = entry.fields.get(name, '') or ''
    return val.strip()


class _MarleyStyle(UnsrtStyle):
    """APS/RevTeX-inspired bibliography style with newest-first sorting."""

    default_sorting_style = 'newest_first'
    default_label_style = 'number'

    def __init__(self, *args, **kwargs):
        # Force abbreviated first names: "F. Lastname" rather than "Firstname Lastname"
        kwargs.setdefault('abbreviate_names', True)
        super().__init__(*args, **kwargs)

    # ------------------------------------------------------------------
    # Author / name formatting
    # ------------------------------------------------------------------

    def format_names(self, role, as_sentence=True):
        """Format names capped at _MAX_AUTHORS, then italic 'et al.'

        Also handles the BibTeX convention of 'and others' as the last author,
        which pybtex parses as a Person whose last name is 'others'.  Such
        entries are treated as implying et al. regardless of count.
        """
        style_self = self

        @node
        def truncated_names(children, context):
            persons = context['entry'].persons.get(role, [])
            if not persons:
                return Text()

            # Detect "and others" — pybtex stores it as a Person with
            # last_name == 'others' and no first names.
            has_others = (
                persons
                and not persons[-1].first_names
                and persons[-1].last_names == ['others']
            )

            if has_others:
                # Strip the trailing 'others' Person; always use et al.
                named_persons = persons[:-1]
                use_etal = True
            elif len(persons) > _MAX_AUTHORS:
                named_persons = persons[:_MAX_AUTHORS]
                use_etal = True
            else:
                named_persons = persons
                use_etal = False

            if use_etal:
                formatted = [
                    style_self.format_name(p, style_self.abbreviate_names)
                    for p in named_persons
                ]
                name_text = join(sep=', ')[formatted].format_data(context)
                etal = RichTag('em', Text(String('et'), Symbol('nbsp'), String('al.')))
                result = Text(name_text, String(' '), etal)
                if as_sentence:
                    return sentence[result].format_data(context)
                return result
            else:
                return UnsrtStyle.format_names(
                    style_self, role, as_sentence=as_sentence
                ).format_data(context)

        return truncated_names()

    # ------------------------------------------------------------------
    # APS-specific helpers
    # ------------------------------------------------------------------

    def _format_collaboration(self, entry):
        """Return ' (Name Collaboration)' Text if the collaboration field is present.

        Ensures the word 'Collaboration' (capitalised) appears exactly once at
        the end.  If the raw field already ends with 'collaboration' (any case),
        it is capitalised in place; otherwise ' Collaboration' is appended.
        """
        collab = _get_field(entry, 'collaboration')
        if not collab:
            return Text()
        # Normalise any existing 'collaboration' variant → 'Collaboration'
        c = re.sub(r'(?i)\bcollaboration\b', 'Collaboration', collab)
        # Append ' Collaboration' if the word is not yet present
        if not re.search(r'\bCollaboration\b', c):
            c += ' Collaboration'
        return Text(String(' ('), Text.from_latex(c), String(')'))

    def _format_doi_href(self, inner_text, doi):
        """Wrap inner_text in a DOI hyperlink."""
        url = 'https://doi.org/' + doi
        return HRef(url, inner_text)

    def _format_arxiv_link(self, eprint, primary_class):
        """Return 'arXiv:XXXX.XXXXX [cat]' as a hyperlink, or empty Text."""
        if not eprint:
            return Text()
        url = 'https://arxiv.org/abs/' + eprint
        label_parts = [String('arXiv:'), String(eprint)]
        if primary_class:
            label_parts += [String(' ['), String(primary_class), String(']')]
        return HRef(url, Text(*label_parts))

    def _format_journal_block(self, entry):
        """Return the journal/volume/pages/year block, DOI-linked if available.

        APS format:        J. Name **vol**, pages (year)
        JHEP/JCAP style:   J. Name **vol** (number)   [no pages field]
        Preprint fallback: (year)

        Returns (journal_block_text, has_journal).  has_journal is True when
        a real journal name was rendered.
        """
        journal = _get_field(entry, 'journal')
        volume  = _get_field(entry, 'volume')
        pages   = _get_field(entry, 'pages')
        number  = _get_field(entry, 'number')
        year    = _get_field(entry, 'year')
        doi     = _get_field(entry, 'doi')

        # An empty string means no journal (journal = "" or journal = {} in BibTeX)
        has_journal = bool(journal and journal != '{}')

        if has_journal:
            # Decode LaTeX in journal name (e.g. '{JCAP}' → 'JCAP')
            journal_rich = Text.from_latex(journal)
            block_parts = [RichTag('em', journal_rich)]

            if volume:
                volume_rich = Text.from_latex(volume)
                block_parts.append(Symbol('nbsp'))
                block_parts.append(RichTag('strong', volume_rich))

            if pages:
                # Normalise '--' or '---' to an en-dash
                pages_normalised = re.sub(r'-{2,}', '\N{EN DASH}', pages)
                pages_rich = Text.from_latex(pages_normalised)
                block_parts.append(String(', '))
                block_parts.append(pages_rich)
            elif number:
                # JHEP/JCAP style: vol (number), no explicit pages field
                number_rich = Text.from_latex(number)
                block_parts.append(Symbol('nbsp'))
                block_parts.append(String('('))
                block_parts.append(number_rich)
                block_parts.append(String(')'))

            if year:
                year_rich = Text.from_latex(year)
                block_parts.append(String(' ('))
                block_parts.append(year_rich)
                block_parts.append(String(')'))

            journal_block = Text(*block_parts)

            if doi:
                journal_block = self._format_doi_href(journal_block, doi)

        else:
            # No journal — render just the year in parens as the block
            has_journal = False
            if year:
                year_rich = Text.from_latex(year)
                journal_block = Text(String('('), year_rich, String(')'))
                if doi:
                    journal_block = self._format_doi_href(journal_block, doi)
            else:
                journal_block = Text()

        return journal_block, has_journal

    # ------------------------------------------------------------------
    # Entry type templates
    # ------------------------------------------------------------------

    def get_article_template(self, e):
        """APS-style article:
        Authors[, (Collaboration)], Title, Journal vol, pages (year)[, arXiv:ID [cat]]
        """
        @node
        def aps_article(children, context):
            entry = context['entry']

            # 1. Authors (with et al. cap); gracefully empty if no author field
            if entry.persons.get('author'):
                author_text = self.format_names('author', as_sentence=False).format_data(context)
            else:
                author_text = Text()

            # 2. Collaboration suffix, e.g. " (DUNE collaboration)"
            collab_text = self._format_collaboration(entry)

            # 3. Title (sentence-capitalised, no quotes)
            title_raw = _get_field(entry, 'title')
            if title_raw:
                title_text = Text.from_latex(title_raw).capitalize()
            else:
                title_text = Text()

            # 4. Journal block (DOI-linked if available)
            journal_block, has_journal = self._format_journal_block(entry)

            # 5. arXiv eprint link
            eprint = _get_field(entry, 'eprint')
            pclass = _get_field(entry, 'primaryclass') or _get_field(entry, 'primaryClass')
            arxiv_link = self._format_arxiv_link(eprint, pclass)

            # --- Assemble ---
            # "Authors[, (Collab)], Title, JournalBlock[, arXiv:...]."
            result_parts = []

            # Author + collaboration block
            author_block = Text(author_text, collab_text)
            result_parts.append(author_block)

            if title_text:
                result_parts.append(String(', '))
                result_parts.append(title_text)

            if journal_block:
                result_parts.append(String(', '))
                result_parts.append(journal_block)

            if arxiv_link:
                result_parts.append(String(', '))
                result_parts.append(arxiv_link)

            # Add terminal period
            full = Text(*result_parts)
            return full.add_period()

        return aps_article()

    def get_phdthesis_template(self, e):
        """APS-style PhD thesis:
        Author, Title (PhD thesis, School, year)[, DOI].
        """
        @node
        def aps_thesis(children, context):
            entry = context['entry']

            author_text = self.format_names('author', as_sentence=False).format_data(context)

            title_raw = _get_field(entry, 'title')
            title_text = Text.from_latex(title_raw).capitalize() if title_raw else Text()

            school = _get_field(entry, 'school')
            year   = _get_field(entry, 'year')
            doi    = _get_field(entry, 'doi')

            # Build parenthetical: (PhD thesis, School, year)
            paren_parts = [String('PhD thesis')]
            if school:
                paren_parts += [String(', '), Text.from_latex(school)]
            if year:
                paren_parts += [String(', '), Text.from_latex(year)]
            paren = Text(String('('), Text(*paren_parts), String(')'))

            if doi:
                paren = self._format_doi_href(paren, doi)

            parts = [author_text]
            if title_text:
                parts += [String(', '), title_text]
            parts += [String(', '), paren]

            return Text(*parts).add_period()

        return aps_thesis()

    def get_mastersthesis_template(self, e):
        """APS-style Master's thesis — same structure as PhD."""
        @node
        def aps_mthesis(children, context):
            entry = context['entry']
            author_text = self.format_names('author', as_sentence=False).format_data(context)
            title_raw = _get_field(entry, 'title')
            title_text = Text.from_latex(title_raw).capitalize() if title_raw else Text()
            school = _get_field(entry, 'school')
            year   = _get_field(entry, 'year')
            doi    = _get_field(entry, 'doi')
            paren_parts = [String("Master's thesis")]
            if school:
                paren_parts += [String(', '), Text.from_latex(school)]
            if year:
                paren_parts += [String(', '), Text.from_latex(year)]
            paren = Text(String('('), Text(*paren_parts), String(')'))
            if doi:
                paren = self._format_doi_href(paren, doi)
            parts = [author_text]
            if title_text:
                parts += [String(', '), title_text]
            parts += [String(', '), paren]
            return Text(*parts).add_period()
        return aps_mthesis()

    def get_misc_template(self, e):
        """Misc entries: author, title (linked if DOI/URL present), (year)."""
        @node
        def aps_misc(children, context):
            entry = context['entry']

            author_text = Text()
            if entry.persons.get('author'):
                author_text = self.format_names('author', as_sentence=False).format_data(context)

            title_raw = _get_field(entry, 'title')
            doi       = _get_field(entry, 'doi')
            url       = _get_field(entry, 'url')
            howpub    = _get_field(entry, 'howpublished')
            year      = _get_field(entry, 'year')
            eprint    = _get_field(entry, 'eprint')
            pclass    = _get_field(entry, 'primaryclass') or _get_field(entry, 'primaryClass')

            parts = []
            if author_text:
                parts.append(author_text)

            if title_raw:
                title_text = Text.from_latex(title_raw).capitalize()
                link_url = doi and ('https://doi.org/' + doi) or url or None
                if link_url:
                    title_text = HRef(link_url, title_text)
                if parts:
                    parts.append(String(', '))
                parts.append(title_text)

            if howpub:
                if parts:
                    parts.append(String(', '))
                parts.append(Text.from_latex(howpub))

            if year:
                if parts:
                    parts.append(String(' '))
                parts.append(Text(String('('), Text.from_latex(year), String(')')))

            arxiv_link = self._format_arxiv_link(eprint, pclass)
            if arxiv_link:
                parts.append(String(', '))
                parts.append(arxiv_link)

            return Text(*parts).add_period() if parts else Text()

        return aps_misc()


pybtex.plugin.register_plugin(
    'pybtex.style.formatting', 'marley_style', _MarleyStyle)

# Use our custom style for all bibliography directives
bibtex_default_style = 'marley_style'

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
