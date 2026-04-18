from typing import Dict, Any, Optional
import io

"""Compatibility wrapper for the Streamlit app.

This function keeps the existing signature used by `app.py` but adds the ability to
return PDF bytes when `path` is None. This lets callers request an in-memory PDF
for preview without duplicating the PDF generation logic.

Signature kept for backward-compatibility:
    export_report_pdf(report_text, blueprint, path, title, author)

If `path` is a string, the exported file path (str) is returned like before. If
`path` is None, the function returns bytes of the generated PDF.
"""


def export_report_pdf(report_text: str, blueprint: Optional[Dict[str, Any]], path: Optional[str] = None, title: str = "Helixar Produktionsbericht", author: Optional[str] = None):
    try:
        from report_generator import export_report_pdf as rg_export
    except Exception as e:
        raise ImportError("PDF-Export nicht verfügbar: report_generator.export_report_pdf konnte nicht geladen werden. Stellen Sie sicher, dass reportlab installiert ist (pip install reportlab) und versuchen Sie es erneut.") from e

    if blueprint is None:
        raise ValueError("Blueprint dict is required for PDF export.")

    # If caller requests in-memory PDF (path is None), use a BytesIO buffer.
    if path is None:
        buf = io.BytesIO()
        rg_export(blueprint, buf, title=title, author=author)
        buf.seek(0)
        return buf.read()

    # Otherwise, write to the provided filesystem path and return that path.
    return rg_export(blueprint, path, title=title, author=author)
