"""Count prose words in Rmd vignettes (excludes code chunks and YAML front matter)."""
import re
from pathlib import Path


def extract_prose(rmd_text):
    """Strip YAML front matter, code chunks, and inline code from Rmd."""
    # Remove YAML front matter
    text = re.sub(r'^---.*?---', '', rmd_text, count=1, flags=re.DOTALL)
    # Remove fenced code chunks (```{r ...} ... ```)
    text = re.sub(r'```\{.*?\}.*?```', '', text, flags=re.DOTALL)
    # Remove any remaining fenced code blocks (``` ... ```)
    text = re.sub(r'```.*?```', '', text, flags=re.DOTALL)
    # Remove inline code (`...`)
    text = re.sub(r'`[^`]+`', '', text)
    # Remove HTML tags
    text = re.sub(r'<[^>]+>', ' ', text)
    # Remove markdown image/link syntax but keep link text
    text = re.sub(r'!\[[^\]]*\]\([^)]*\)', '', text)  # images
    text = re.sub(r'\[([^\]]*)\]\([^)]*\)', r'\1', text)  # links
    # Remove markdown heading markers
    text = re.sub(r'^#{1,6}\s+', '', text, flags=re.MULTILINE)
    # Collapse whitespace
    text = re.sub(r'\s+', ' ', text).strip()
    return text


root = Path(__file__).resolve().parent.parent
vignettes = root / "vignettes"

for rmd_file in sorted(vignettes.glob("*.Rmd")):
    text = extract_prose(rmd_file.read_text(encoding="utf-8", errors="ignore"))
    words = len(text.split())
    print(f"{rmd_file.stem:30s}  {words:6,d} words")

total = sum(
    len(extract_prose(f.read_text(encoding="utf-8", errors="ignore")).split())
    for f in vignettes.glob("*.Rmd")
)
print(f"\n{'TOTAL':30s}  {total:6,d} words")
