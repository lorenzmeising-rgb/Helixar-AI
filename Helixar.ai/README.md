Helixar — Produktions Blueprint Generator

Helixar ist ein kleines, regelbasiertes Entscheidungsunterstützungs‑Tool, das Gründern und Teams bei frühen Machbarkeitsbewertungen für mikrobiell hergestellte Stoffe hilft. Der Dienst liefert beschreibende, nicht‑operative Produktions‑Blueprints (keine Laborprotokolle).

Installation

1. Python 3.9+ verwenden
2. Abhängigkeiten installieren (Report-Export optional):

```bash
pip install -r requirements.txt
# optional: für PDF-Export
pip install reportlab
```

Starten der App

```bash
streamlit run app.py
```

PDF-Export

Der PDF-Export nutzt reportlab. Wenn reportlab nicht installiert ist, läuft die App weiterhin, zeigt aber eine Fehlermeldung beim Versuch, ein PDF zu erzeugen. Installieren Sie reportlab wie oben gezeigt, um PDF-Reports zu aktivieren.

Sicherheitshinweis

Die erzeugten Reports sind bewusst nicht-operativ und enthalten keine Laboranweisungen, Konzentrationen, Zeitangaben oder andere detaillierte Protocolle. Sie dienen ausschließlich strategischer Entscheidungsunterstützung und ersetzen keine Validierung, Sicherheits- oder regulatorische Beratung.
