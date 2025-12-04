from __future__ import annotations

import json
import os
from typing import Any, Dict, List

from jinja2 import BaseLoader, Environment, select_autoescape

from .chemistry import build_structure_payload

try:  # pragma: no cover - optional dependency
    import pdfkit
except Exception as exc:  # noqa: BLE001
    pdfkit = None  # type: ignore[assignment]
    _PDFKIT_IMPORT_ERROR: Exception | None = exc
    _PDFKIT_CONFIG = None
    _PDFKIT_CONFIG_ERROR: Exception | None = exc
else:
    _PDFKIT_IMPORT_ERROR = None
    wkhtmltopdf_path = os.getenv("WKHTMLTOPDF_PATH")
    try:
        _PDFKIT_CONFIG = (
            pdfkit.configuration(wkhtmltopdf=wkhtmltopdf_path)
            if wkhtmltopdf_path
            else pdfkit.configuration()
        )
        _PDFKIT_CONFIG_ERROR: Exception | None = None
    except (OSError, IOError) as exc:  # wkhtmltopdf binary missing
        _PDFKIT_CONFIG = None
        _PDFKIT_CONFIG_ERROR = exc

from .schemas import JobResponse

_ENV = Environment(
    loader=BaseLoader(),
    autoescape=select_autoescape(enabled_extensions=("html",)),
    trim_blocks=True,
    lstrip_blocks=True,
)

_REPORT_TEMPLATE = _ENV.from_string(
    """<!DOCTYPE html>
<html lang=\"en\">
  <head>
    <meta charset=\"utf-8\" />
    <title>PhaGen Agentic Report</title>
    <style>
      :root {
        font-family: 'Inter', 'Segoe UI', sans-serif;
        color: #0f1115;
      }
      body {
        margin: 0;
        padding: 32px 40px;
        background: #f7f7fb;
      }
      header {
        border-bottom: 2px solid #1c1f2b;
        margin-bottom: 24px;
        padding-bottom: 12px;
      }
      h1 {
        font-size: 28px;
        margin: 0;
      }
      .eyebrow {
        text-transform: uppercase;
        letter-spacing: 0.12em;
        font-size: 11px;
        color: #5f6476;
        margin: 0 0 4px;
      }
      .summary-box {
        background: #fff;
        border-radius: 16px;
        padding: 20px;
        border: 1px solid #e0e3ed;
        margin-bottom: 24px;
      }
      .structure-card {
        display: flex;
        flex-direction: column;
        gap: 12px;
      }
      .structure-wrapper {
        background: #fff;
        border-radius: 12px;
        padding: 12px;
        border: 1px solid #e0e3ed;
      }
      .structure-wrapper svg {
        width: 100%;
        height: auto;
      }
      .structure-meta {
        font-size: 11px;
        color: #5f6476;
      }
      .structure-error {
        color: #b42318;
        font-size: 12px;
      }
      .recommendation-chip {
        display: inline-block;
        padding: 4px 12px;
        border-radius: 999px;
        font-size: 12px;
        font-weight: 600;
        color: #0f1115;
        background: #d0f0ff;
        text-transform: uppercase;
      }
      .worker-section {
        margin-bottom: 28px;
        page-break-inside: avoid;
      }
      .worker-section h2 {
        margin: 0 0 8px;
        font-size: 20px;
      }
      .meta-list, .evidence-list {
        margin: 12px 0;
        padding-left: 18px;
      }
      .meta-list li, .evidence-list li {
        margin-bottom: 6px;
      }
      .meta-key {
        font-weight: 600;
      }
      .validation-list {
        margin: 16px 0 0;
        padding-left: 18px;
      }
      .validation-list li {
        margin-bottom: 10px;
      }
      .validation-status {
        font-weight: 600;
      }
      footer {
        margin-top: 40px;
        font-size: 11px;
        color: #6f7280;
      }
    </style>
  </head>
  <body>
    <header>
      <p class=\"eyebrow\">PhaGen Agentic · Innovation Story</p>
      <h1>{{ molecule }}</h1>
      <p>Job ID {{ job_id }} · Report V{{ report_version }} · Generated {{ generated_at }}</p>
    </header>

    <section class=\"summary-box\">
      <p class=\"eyebrow\">Recommendation</p>
      <div class=\"recommendation-chip\">{{ recommendation }}</div>
      <p style=\"margin-top: 12px;\">Market score: {{ market_score }}</p>
      <p style=\"margin-top: 12px;\">{{ innovation_story }}</p>
    </section>

    {% if structure %}
    <section class="summary-box structure-card">
      <p class="eyebrow">Molecular structure</p>
      {% if structure.svg %}
      <div class="structure-wrapper">
        {{ structure.svg | safe }}
      </div>
      {% endif %}
      {% if structure.path %}
      <p class="structure-meta">Asset saved to {{ structure.path }}</p>
      {% endif %}
      {% if structure.error %}
      <p class="structure-error">Structure unavailable: {{ structure.error }}</p>
      {% endif %}
    </section>
    {% endif %}

    {% if validation %}
    <section class="summary-box">
      <p class="eyebrow">Claim traceability · {{ validation.status_label }}</p>
      <p>{{ validation.claims_linked }} / {{ validation.claims_total }} claims linked to evidence.</p>
      <ol class="validation-list">
        {% for claim in validation.claims %}
        <li>
          <span class="validation-status">{{ claim.claim_text }}</span><br />
          Evidence: {{ claim.evidence_label }}
        </li>
        {% endfor %}
      </ol>
    </section>
    {% endif %}

    {% for worker in worker_sections %}
    <section class=\"worker-section\">
      <p class=\"eyebrow\">{{ worker.confidence_band }} confidence</p>
      <h2>{{ worker.name }}</h2>
      <p>{{ worker.summary }}</p>
      {% if worker.metadata %}
      <ul class=\"meta-list\">
        {% for meta in worker.metadata %}
        <li><span class=\"meta-key\">{{ meta.key }}:</span> {{ meta.value }}</li>
        {% endfor %}
      </ul>
      {% endif %}
      {% if worker.evidence %}
      <ol class=\"evidence-list\">
        {% for evidence in worker.evidence %}
        <li>
          <span class=\"meta-key\">{{ evidence.type }}</span> · {{ evidence.text }}
          {% if evidence.url %}
          <br />Source: {{ evidence.url }}
          {% endif %}
        </li>
        {% endfor %}
      </ol>
      {% endif %}
    </section>
    {% endfor %}

    <footer>
      Generated by PhaGen Agentic · Confidential · {{ generated_at }}
    </footer>
  </body>
</html>
""",
)

def _build_structure_block(job: JobResponse, payload: Dict[str, Any]) -> Dict[str, str] | None:
  existing = payload.get("structure")
  if existing:
    return existing

  smiles = payload.get("smiles") or payload.get("molecule_smiles")
  if not smiles:
    return None

  molecule_label = (payload.get("molecule") or "molecule").strip() or "molecule"
  structure = build_structure_payload(
    smiles=smiles,
    job_id=job.job_id,
    molecule_label=molecule_label,
  )
  payload["structure"] = structure
  return structure

_BAND_LABELS = {
    "low": "Low",
    "medium": "Medium",
    "high": "High",
}


def _truncate(value: str, limit: int = 480) -> str:
    if len(value) <= limit:
        return value
    return value[: limit - 3] + "..."


def _format_metadata(metadata: Dict[str, Any]) -> List[Dict[str, str]]:
    items: List[Dict[str, str]] = []
    for key, value in metadata.items():
        if isinstance(value, (dict, list)):
            value_str = json.dumps(value, ensure_ascii=False)
        else:
            value_str = str(value)
        items.append(
            {
                "key": key.replace("_", " ").title(),
                "value": _truncate(value_str, 240),
            }
        )
    return items


def _format_evidence(evidence: List[Dict[str, Any]]) -> List[Dict[str, str]]:
    entries: List[Dict[str, str]] = []
    for record in (evidence or [])[:5]:
        entries.append(
            {
                "type": record.get("type", "Evidence"),
                "text": _truncate(record.get("text", "")),
                "url": record.get("url", ""),
                "evidence_id": record.get("evidence_id", ""),
            }
        )
    return entries


def _format_claim_links(validation: Dict[str, Any]) -> Dict[str, Any] | None:
    if not validation:
        return None
    claims = []
    for claim in validation.get("claim_links", []):
        evidence_ids = claim.get("evidence_ids") or []
        label = ", ".join(evidence_ids) if evidence_ids else "Not linked"
        claims.append(
            {
                "claim_text": _truncate(claim.get("claim_text", ""), 320),
                "evidence_label": label,
            }
        )
    return {
        "status_label": (validation.get("status") or "Needs review").replace("_", " ").title(),
        "claims_total": validation.get("claims_total", 0),
        "claims_linked": validation.get("claims_linked", 0),
        "claims": claims,
    }


def _build_context(job: JobResponse) -> Dict[str, object]:
    if not job.payload:
        raise ValueError("Job payload is empty; cannot render report")
    payload: Dict[str, Any] = job.payload or {}
    workers: List[Dict[str, Any]] = []
    for name, data in (payload.get("workers") or {}).items():
        metadata = data.get("metadata") or {}
        evidence = data.get("evidence") or []
        workers.append(
            {
                "name": name.title(),
                "summary": data.get("summary", ""),
                "confidence_band": _BAND_LABELS.get(
                    (data.get("confidence_band") or "low").lower(),
                    "Confidence",
                ),
                "metadata": _format_metadata(metadata),
                "evidence": _format_evidence(evidence),
            }
        )
    workers.sort(key=lambda item: item["name"])

    molecule = payload.get("molecule") or "Molecule report"
    recommendation = payload.get("recommendation", "Investigate")
    if hasattr(recommendation, "value"):
        recommendation = getattr(recommendation, "value")
    market_score = payload.get("market_score", 0)
    try:
        market_score = int(market_score)
    except (ValueError, TypeError):
        market_score = 0

    validation_block = _format_claim_links(payload.get("validation") or {})
    report_version = payload.get("report_version") or job.report_version or 1
    structure_block = _build_structure_block(job, payload)

    return {
        "molecule": molecule,
        "job_id": job.job_id,
        "generated_at": job.updated_at.strftime("%b %d, %Y %H:%M UTC"),
        "report_version": report_version,
        "innovation_story": payload.get("innovation_story", ""),
        "recommendation": recommendation,
        "market_score": market_score,
        "worker_sections": workers,
        "validation": validation_block,
        "structure": structure_block,
    }


def render_report_html(job: JobResponse) -> str:
    context = _build_context(job)
    return _REPORT_TEMPLATE.render(**context)


def generate_report_pdf(job: JobResponse) -> bytes:
    html = render_report_html(job)

    if pdfkit and _PDFKIT_CONFIG:
        try:
            return pdfkit.from_string(html, False, configuration=_PDFKIT_CONFIG)
        except (OSError, IOError) as exc:  # noqa: BLE001
            raise RuntimeError(
                "wkhtmltopdf failed to generate a PDF. Ensure the wkhtmltopdf binary is installed "
                "and accessible (set WKHTMLTOPDF_PATH if it's not on PATH)."
            ) from exc

    raise RuntimeError(
        "wkhtmltopdf is not configured. Install wkhtmltopdf from https://wkhtmltopdf.org/downloads.html "
        "and ensure pdfkit can locate it (add to PATH or set WKHTMLTOPDF_PATH). "
        f"pdfkit import error: {_PDFKIT_IMPORT_ERROR}; configuration error: {_PDFKIT_CONFIG_ERROR}"
    )
