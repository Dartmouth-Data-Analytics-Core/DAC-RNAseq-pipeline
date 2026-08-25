#!/usr/bin/env python3
"""
Rendering for the dashboard's QC verdict section.

Split out from build_dashboard.py so it depends on nothing but the standard library.
That matters because scripts/add_qc_report.py imports it: adding a verdict table to an
already-built dashboard should not require numpy and pandas, since the common case is a
delivered results folder holding little more than the HTML itself.
"""

import html
import os
import re


def _warn(msg):
    import sys
    print(f"[dashboard] WARNING: {msg}", file=sys.stderr)


def _first_existing(*paths):
    for p in paths:
        if p and os.path.exists(p):
            return p
    return None


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# QC ASSESSMENT SECTION
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The narrative QC write-up is authored outside this script - by hand, or by the
# rnaseq-qc skill - as Markdown, and rendered in here. Markdown rather than HTML
# because the author should be writing prose, not markup: a stray tag would
# otherwise land in a client deliverable, and escaping the text means an authored
# report can never break the page.

ASSESSMENT_START = "<!-- QC-ASSESSMENT:START -->"
ASSESSMENT_END = "<!-- QC-ASSESSMENT:END -->"

VERDICT_CLASS = {
    "PASS": "pass",
    "PASS WITH CAVEATS": "caveat",
    "PASS WITH CAVEAT": "caveat",
    "FAIL": "fail",
}

_SAFE_LINK = re.compile(r"^(https?://|mailto:|#)", re.I)
_BLOCK_START = re.compile(r"^(#{1,4}\s|```|[-*+]\s|\d+[.)]\s|>\s)")


def _md_inline(text):
    """Inline Markdown -> HTML. Escapes first, so authored text can never inject markup."""
    out = html.escape(text, quote=False)
    out = re.sub(r"`([^`]+)`", r"<code>\1</code>", out)
    out = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", out)
    out = re.sub(r"(?<!\*)\*([^*\n]+)\*(?!\*)", r"<em>\1</em>", out)

    def link(m):
        label, href = m.group(1), m.group(2).strip()
        if not _SAFE_LINK.match(href):
            return label
        return f'<a href="{html.escape(href, quote=True)}">{label}</a>'

    return re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)", link, out)


def _verdict_cell(cell):
    """A cell holding only a verdict renders as a badge, so the table scans at a glance."""
    key = re.sub(r"\s+", " ", cell.strip()).upper().rstrip(".")
    cls = VERDICT_CLASS.get(key)
    if not cls:
        return _md_inline(cell)
    return f'<span class="verdict {cls}">{html.escape(cell.strip())}</span>'


_NUMERIC = re.compile(r"^[<>~=+-]?[\d.,]+\s*(%|[KMG]|bp|reads|x)?$", re.I)


def _md_table(rows):
    """
    Render a pipe table, aligning each column by what it holds.

    The dashboard's own tables are all numeric so the stylesheet right-aligns cells;
    an authored table mixes prose and numbers, and right-aligned prose reads badly.
    A column is right-aligned only when every body cell in it looks numeric.
    """
    head, body = rows[0], rows[2:]
    ncol = len(head)
    numeric = []
    for c in range(ncol):
        vals = [r[c].strip() for r in body if c < len(r) and r[c].strip() not in ("", "-", "\u2014")]
        numeric.append(bool(vals) and all(_NUMERIC.match(v) for v in vals))

    def align(c):
        return "" if numeric[c] else " left"

    out = ['<div class="tablewrap"><table class="prose-table"><thead><tr class="cols">']
    out += [f'<th class="nosort{align(c)}">{_md_inline(h)}</th>' for c, h in enumerate(head)]
    out.append("</tr></thead><tbody>")
    for r in body:
        out.append("<tr>")
        for c, cell in enumerate(r):
            cls = ("sticky-l" if c == 0 else "") + (align(c) if c else "")
            attr = f' class="{cls.strip()}"' if cls.strip() else ""
            out.append(f"<td{attr}>{_verdict_cell(cell)}</td>")
        out.append("</tr>")
    out.append("</tbody></table></div>")
    return "".join(out)


def render_markdown(text):
    """
    A deliberately small Markdown subset: headings, paragraphs, lists, tables,
    blockquotes, fenced code, and inline emphasis/code/links. Enough for a QC
    write-up and a client email, with no dependency and no raw-HTML passthrough.
    """
    lines = text.replace("\r\n", "\n").split("\n")
    out, i = [], 0

    def is_row(s):
        return s.startswith("|") and "|" in s[1:]

    while i < len(lines):
        stripped = lines[i].strip()

        if stripped.startswith("```"):
            i += 1
            buf = []
            while i < len(lines) and not lines[i].strip().startswith("```"):
                buf.append(lines[i])
                i += 1
            i += 1
            out.append("<pre><code>" + html.escape("\n".join(buf)) + "</code></pre>")
            continue

        m = re.match(r"^(#{1,4})\s+(.*)$", stripped)
        if m:
            lvl = min(4, len(m.group(1)) + 1)   # '#' becomes h2; the page owns the h1
            out.append(f"<h{lvl}>{_md_inline(m.group(2))}</h{lvl}>")
            i += 1
            continue

        if is_row(stripped) and i + 1 < len(lines) and re.match(r"^\|[\s:|-]+\|?$", lines[i + 1].strip()):
            rows = []
            while i < len(lines) and is_row(lines[i].strip()):
                rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")])
                i += 1
            out.append(_md_table(rows))
            continue

        m = re.match(r"^([-*+]|\d+[.)])\s+(.*)$", stripped)
        if m:
            tag = "ul" if m.group(1) in ("-", "*", "+") else "ol"
            items = []
            while i < len(lines):
                mm = re.match(r"^([-*+]|\d+[.)])\s+(.*)$", lines[i].strip())
                if not mm:
                    break
                parts = [mm.group(2)]
                i += 1
                # Lazy continuation: a wrapped list item is normal Markdown, and without
                # this the remainder of the line escapes the list as its own paragraph.
                while (i < len(lines) and lines[i].strip()
                       and not _BLOCK_START.match(lines[i].strip())
                       and not is_row(lines[i].strip())):
                    parts.append(lines[i].strip())
                    i += 1
                items.append(f"<li>{_md_inline(' '.join(parts))}</li>")
            out.append(f"<{tag}>" + "".join(items) + f"</{tag}>")
            continue

        if stripped.startswith(">"):
            buf = []
            while i < len(lines) and lines[i].strip().startswith(">"):
                buf.append(lines[i].strip().lstrip(">").strip())
                i += 1
            out.append("<blockquote>" + _md_inline(" ".join(buf)) + "</blockquote>")
            continue

        if not stripped:
            i += 1
            continue

        buf = []
        while (i < len(lines) and lines[i].strip()
               and not _BLOCK_START.match(lines[i].strip())
               and not is_row(lines[i].strip())):
            buf.append(lines[i].strip())
            i += 1
        if buf:
            out.append("<p>" + _md_inline(" ".join(buf)) + "</p>")

    return "\n".join(out)


def find_qc_report(run_dir, explicit=None):
    """Explicit path wins; otherwise pick up a report sitting in the run directory."""
    if explicit == "":
        return None
    if explicit:
        if os.path.exists(explicit):
            return explicit
        _warn(f"QC report not found, section omitted: {explicit}")
        return None
    return _first_existing(
        os.path.join(run_dir, "qc_assessment.md"),
        os.path.join(run_dir, "qc_assessment.html"),
    )


def build_assessment(path):
    """Render the report into the section the template expects, or return ''."""
    if not path:
        return ""
    raw = open(path, encoding="utf-8", errors="ignore").read().strip()
    if not raw:
        return ""
    body = raw if path.lower().endswith((".html", ".htm")) else render_markdown(raw)
    return (
        ASSESSMENT_START
        + '\n  <section id="assessment">\n'
        '    <div class="sec-head"><div>'
        '<div class="eyebrow">Quality assessment</div>'
        "<h2>Sample QC verdicts</h2></div></div>\n"
        '    <div class="card card-pad prose">\n' + body + "\n    </div>\n"
        "  </section>\n  " + ASSESSMENT_END
    )
