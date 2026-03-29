"""Compare prose word counts between INLAocc and spOccupancy vignettes."""
import re
from pathlib import Path


def extract_prose(rmd_text):
    """Strip YAML front matter, code chunks, and inline code from Rmd."""
    text = re.sub(r'^---.*?---', '', rmd_text, count=1, flags=re.DOTALL)
    text = re.sub(r'```\{.*?\}.*?```', '', text, flags=re.DOTALL)
    text = re.sub(r'```.*?```', '', text, flags=re.DOTALL)
    text = re.sub(r'`[^`]+`', '', text)
    text = re.sub(r'<[^>]+>', ' ', text)
    text = re.sub(r'!\[[^\]]*\]\([^)]*\)', '', text)
    text = re.sub(r'\[([^\]]*)\]\([^)]*\)', r'\1', text)
    text = re.sub(r'^#{1,6}\s+', '', text, flags=re.MULTILINE)
    text = re.sub(r'\s+', ' ', text).strip()
    return text


def count_dir(vdir):
    counts = {}
    for f in sorted(vdir.glob("*.Rmd")):
        prose = extract_prose(f.read_text(encoding="utf-8", errors="ignore"))
        counts[f.stem] = len(prose.split())
    return counts


inlaocc_dir = Path(__file__).resolve().parent.parent / "vignettes"
spocc_dir = Path("C:/Users/Gilles Colling/Documents/dev/spOccupancy_src/vignettes")

inla = count_dir(inlaocc_dir)
spocc = count_dir(spocc_dir)

# Topic mapping: (topic label, INLAocc vignette, spOccupancy vignette(s))
mapping = [
    ("Quickstart / model fitting",  "quickstart",            "modelFitting"),
    ("Data formatting",             "data-formatting",       "dataFormatting"),
    ("Random effects",              "random-effects",        "randomEffects"),
    ("Spatial models",              "spatial-models",        "spaceTimeModels"),
    ("Temporal models",             "temporal-models",       "spaceTimeModelsHTML"),
    ("SVC models",                  "svc-models",            "svcModels"),
    ("SVC interpretation",          None,                    "svcInterpretation"),
    ("Multi-species",               "multi-species",         "factorModels"),
    ("Integrated models",           "integrated-models",     "integratedMultispecies"),
    ("Identifiability",             "identifiability",       "identifiability"),
    ("Diagnostics",                 "diagnostics",           None),
    ("Algorithm details",           "algorithm-details",     None),
    ("Migration guide",             "spoccupancy-migration", None),
    ("MCMC samplers",               None,                    "mcmcSamplers"),
    ("MCMC factor models",          None,                    "mcmcFactorModels"),
    ("MCMC SVC models",             None,                    "mcmcSVCModels"),
]

print(f"{'Topic':<30s}  {'INLAocc':>8s}  {'spOcc':>8s}  {'Diff':>8s}")
print("-" * 60)

inla_total = 0
spocc_total = 0

for topic, inla_key, spocc_key in mapping:
    iw = inla.get(inla_key, 0) if inla_key else 0
    sw = spocc.get(spocc_key, 0) if spocc_key else 0
    inla_total += iw
    spocc_total += sw
    diff = iw - sw
    i_str = f"{iw:,d}" if iw else "—"
    s_str = f"{sw:,d}" if sw else "—"
    d_str = f"+{diff:,d}" if diff > 0 else (f"{diff:,d}" if diff < 0 else "—")
    print(f"{topic:<30s}  {i_str:>8s}  {s_str:>8s}  {d_str:>8s}")

print("-" * 60)
print(f"{'TOTAL':<30s}  {inla_total:>7,d}   {spocc_total:>7,d}   {inla_total - spocc_total:>+7,d}")
