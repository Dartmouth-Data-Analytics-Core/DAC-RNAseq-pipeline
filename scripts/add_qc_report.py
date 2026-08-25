#!/usr/bin/env python3
"""
Add (or replace) the QC report section in an already-built dashboard.

The pipeline's own path is `build_dashboard.py --qc-report`, which renders the report
at build time. This script exists for the case where the dashboard HTML is all you
have - a delivered results folder, or a report written after the run finished - so the
section can be added without re-running the pipeline.

The insert is bounded by HTML comment markers, so re-running replaces the previous
section rather than stacking another copy. Nothing else in the page is touched: the
embedded data payload and the chart runtime are never re-serialised.

Usage
-----
    python scripts/add_qc_report.py RNAseq_Dashboard.html qc_assessment.md
    python scripts/add_qc_report.py RNAseq_Dashboard.html qc_assessment.md -o Copy.html
    python scripts/add_qc_report.py RNAseq_Dashboard.html --remove
"""

import argparse
import os
import re
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dashboard_report import (  # noqa: E402
    ASSESSMENT_END,
    ASSESSMENT_START,
    build_assessment,
)


def strip_existing(html_text):
    """Remove a previously inserted section so repeat runs replace rather than stack."""
    pattern = re.compile(
        re.escape(ASSESSMENT_START) + ".*?" + re.escape(ASSESSMENT_END), re.S
    )
    return pattern.sub("", html_text), bool(pattern.search(html_text))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("dashboard", help="Built dashboard HTML to modify")
    ap.add_argument("report", nargs="?", help="Markdown (or HTML) QC write-up")
    ap.add_argument("-o", "--output", default=None,
                    help="Write here instead of editing in place")
    ap.add_argument("--remove", action="store_true",
                    help="Strip any existing QC report section and stop")
    ap.add_argument("--no-backup", action="store_true",
                    help="Skip the .bak copy taken before an in-place edit")
    args = ap.parse_args(argv)

    if not os.path.exists(args.dashboard):
        sys.exit(f"[add_qc_report] ERROR: no such dashboard: {args.dashboard}")
    if not args.remove and not args.report:
        sys.exit("[add_qc_report] ERROR: a report file is required (or pass --remove)")
    if args.report and not os.path.exists(args.report):
        sys.exit(f"[add_qc_report] ERROR: no such report: {args.report}")

    page = open(args.dashboard, encoding="utf-8", errors="ignore").read()
    page, replaced = strip_existing(page)

    if args.remove:
        section = ""
    else:
        section = build_assessment(args.report)
        if not section:
            sys.exit(f"[add_qc_report] ERROR: {args.report} is empty; nothing to insert")

    # Before </main> keeps the report inside the page body, above the footer. Falling
    # back to </body> covers a hand-edited page that has lost its <main>.
    anchor = "</main>" if "</main>" in page else "</body>"
    if anchor not in page:
        sys.exit("[add_qc_report] ERROR: could not find an insertion point in the page")
    page = page.replace(anchor, ("  " + section + "\n" if section else "") + anchor, 1)

    out = args.output or args.dashboard
    if out == args.dashboard and not args.no_backup:
        shutil.copy2(args.dashboard, args.dashboard + ".bak")

    with open(out, "w", encoding="utf-8") as fh:
        fh.write(page)

    verb = "removed" if args.remove else ("replaced" if replaced else "added")
    print(f"[add_qc_report] {verb} QC report section -> {out} "
          f"({os.path.getsize(out) / 1024:,.0f} KB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
