"""Ensure the curated corpus points to actual executed Rust test names."""
import json
from pathlib import Path
import re
import subprocess
root = Path(__file__).resolve().parents[1]
corpus = json.loads((root / "corpus.json").read_text())
listed = subprocess.check_output(
    ["cargo", "test", "--offline", "--", "--list"], cwd=root, text=True
)
executed_names = set(re.findall(r"^([a-zA-Z0-9_]+): test$", listed, re.M))
ids = set()
for row in corpus["entries"]:
    assert row["id"] not in ids, row
    ids.add(row["id"])
    assert row["test"] in executed_names, row
    text = (root / row["file"]).read_text()
    assert re.search(r"fn\s+" + re.escape(row["test"]) + r"\s*\(", text), row
for row in corpus["additional_checks"]:
    assert row["test"] in executed_names, row
print(json.dumps({
    "status": "PASS", "textbook_entries": len(ids),
    "listed_rust_tests": len(executed_names),
    "warning": "Source anchors and coverage are curated, not automatically certified.",
}))
