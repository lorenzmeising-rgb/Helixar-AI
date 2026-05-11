"""
feedback_submission.py
======================
Sendet User-Feedback an Formsubmit.co (https://formsubmit.co).

Workflow:
- POST an https://formsubmit.co/ajax/<EMAIL> mit Form-Daten.
- Beim allererste Submit schickt Formsubmit eine Verifizierungs-Mail an
  die Empfaenger-Adresse; sobald der Link einmalig geklickt wurde, gehen
  alle weiteren Submits direkt durch.
- Keine Credentials, kein API-Key, kein Account noetig.

Bewusste Entscheidung gegen SMTP / Resend (siehe Memory-Files):
- Pitch-Phase, niedriges Feedback-Volumen
- Keine Secrets in Streamlit Cloud noetig (geringeres Leak-Risiko)
- Migration auf Resend / SMTP spaeter ohne UI-Aenderungen moeglich
"""

from typing import Dict, Any, Optional
import requests

# Empfaenger-Adresse fuer alle Feedback-Submits
FEEDBACK_RECIPIENT = "lorenzmeising@icloud.com"

# Formsubmit-AJAX-Endpoint (gibt JSON statt HTML zurueck)
FORMSUBMIT_ENDPOINT = f"https://formsubmit.co/ajax/{FEEDBACK_RECIPIENT}"

# Timeout fuer den POST-Request in Sekunden
REQUEST_TIMEOUT = 10

# Formsubmit prueft Origin/Referer und lehnt Requests ohne sie als
# "file://-aus-Browser" ab. Wir setzen die Live-Domain — funktioniert
# sowohl lokal als auch auf Streamlit Cloud.
ORIGIN_HEADER = "https://helixar-ai.streamlit.app"


def submit_feedback(
    feedback_text: str,
    category: str,
    name: Optional[str] = None,
    company: Optional[str] = None,
    user_email: Optional[str] = None,
    app_language: str = "de",
    extra_context: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """
    Submit user feedback to Formsubmit.co.

    Parameters
    ----------
    feedback_text : str
        The actual feedback content (required, non-empty).
    category : str
        Category label (e.g. "Bug", "Feature-Request", ...).
    name : Optional[str]
        Submitter name (optional).
    company : Optional[str]
        Submitter company / organisation (optional).
    user_email : Optional[str]
        Submitter email for follow-up (optional).
    app_language : str
        Current UI language ("de" or "en") — included as context.
    extra_context : Optional[Dict[str, Any]]
        Additional key/value pairs to include in the email body
        (e.g. current page, browser info).

    Returns
    -------
    dict with keys:
      - ok   : bool
      - message : str (human-readable status)
      - http_status : int or None
    """
    if not feedback_text or not str(feedback_text).strip():
        return {
            "ok": False,
            "message": "empty_feedback",
            "http_status": None,
        }

    payload = {
        # Formsubmit-Spezialfelder (alle mit _ prefixed)
        "_subject": f"[Helixar AI Feedback] {category}",
        "_template": "table",
        "_captcha": "false",
        # honeypot — Spam-Bots fuellen alle Felder, echte User lassen _honey leer
        "_honey": "",

        # Eigentliche Inhalte
        "Kategorie": category,
        "Sprache_App": app_language,
        "Name": (name or "").strip() or "(nicht angegeben)",
        "Firma": (company or "").strip() or "(nicht angegeben)",
        "Email_des_Absenders": (user_email or "").strip() or "(nicht angegeben)",
        "Feedback": str(feedback_text).strip(),
    }

    if extra_context:
        for k, v in extra_context.items():
            # Formsubmit zeigt Keys 1:1 im Mail-Body — Underscores statt Spaces
            safe_key = str(k).replace(" ", "_")
            payload[safe_key] = str(v)

    try:
        response = requests.post(
            FORMSUBMIT_ENDPOINT,
            data=payload,
            timeout=REQUEST_TIMEOUT,
            headers={
                "Accept": "application/json",
                # Formsubmit lehnt ohne Origin/Referer mit
                # "FormSubmit will not work in pages browsed as HTML files" ab.
                "Origin": ORIGIN_HEADER,
                "Referer": ORIGIN_HEADER + "/",
            },
        )
    except requests.exceptions.Timeout:
        return {
            "ok": False,
            "message": "timeout",
            "http_status": None,
            "upstream_message": None,
        }
    except requests.exceptions.ConnectionError:
        return {
            "ok": False,
            "message": "connection_error",
            "http_status": None,
            "upstream_message": None,
        }
    except Exception as e:
        return {
            "ok": False,
            "message": f"unexpected_error: {type(e).__name__}",
            "http_status": None,
            "upstream_message": None,
        }

    upstream_message = None
    if response.status_code == 200:
        # Formsubmit antwortet mit JSON {"success": "true"|"false", "message": "..."}
        try:
            data = response.json()
            success_flag = str(data.get("success", "")).lower() == "true"
            upstream_message = data.get("message")
        except Exception:
            # Falls kein gueltiges JSON: Status 200 als Erfolg werten
            success_flag = True

        # Special case: Form noch nicht aktiviert (erster Submit nach Setup).
        # Formsubmit hat eine Aktivierungs-Mail an FEEDBACK_RECIPIENT geschickt;
        # nach einem Klick auf den Link funktioniert das Formular dann sofort.
        if not success_flag and upstream_message and "Activation" in upstream_message:
            return {
                "ok": False,
                "message": "needs_activation",
                "http_status": response.status_code,
                "upstream_message": upstream_message,
            }

        return {
            "ok": success_flag,
            "message": "submitted" if success_flag else "rejected_by_service",
            "http_status": response.status_code,
            "upstream_message": upstream_message,
        }

    return {
        "ok": False,
        "message": "http_error",
        "http_status": response.status_code,
        "upstream_message": None,
    }
