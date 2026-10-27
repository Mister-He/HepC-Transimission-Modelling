# =============================================================================
# registry.py — evidence and claim bookkeeping for the manuscript build, in the format of the
# scientific-writing skill: every factual or numeric line gets a claim ID (C...), bound to
# evidence IDs (E...) for the bibliography entries and repository files it rests on.
# Writes an audit copy of each document (with [claim:...] [evidence:...] markers), the clean
# Markdown for pandoc, claims.csv, source_manifest.json and the other registries.
# =============================================================================
import csv, hashlib, json, os, re

CITE_RE = re.compile(r"@([A-Za-z0-9_:-]+)")
ERRORS = []  # lines that need evidence; the build stops if any


def claim_hash(text):  # same normalisation as the skill's audit_claims.claim_text_hash
    return hashlib.sha256(" ".join(text.split()).encode("utf-8")).hexdigest()


def parse_bib(path):
    out = {}
    for m in re.finditer(r"@(\w+)\{([^,]+),(.*?)\n\}", open(path).read(), re.S):
        body = m.group(3) + "\n"  # the last field has no trailing newline inside the match
        def g(k):
            x = re.search(rf"\b{k}\s*=\s*[{{\"](.*?)[}}\"]\s*,?\s*\n", body, re.S | re.I)
            return re.sub(r"[{}]", "", x.group(1)).strip() if x else ""
        out[m.group(2).strip()] = {"type": m.group(1).lower(), "title": g("title"), "year": g("year"),
                                   "doi": g("doi"), "url": g("url"), "author": g("author")}
    return out


BIB_TYPE = {"article": "journal_article", "book": "book", "inproceedings": "conference_paper",
            "incollection": "chapter", "techreport": "report", "misc": "webpage", "online": "webpage"}


class Registry:
    def __init__(self, root, bib_path, files, restricted=()):
        self.root, self.bib, self.files = root, parse_bib(bib_path), files
        self.restricted = set(restricted)
        self.sources, self.claims = {}, []

    def eid(self, key):
        """Evidence ID for a bibliography key or a Facts.files key (created on first use)."""
        if key in self.sources:
            return self.sources[key]["evidence_id"]
        n = len(self.sources) + 1
        e = f"E{n:03d}"
        if key in self.bib:
            b = self.bib[key]
            src = {"source_type": BIB_TYPE.get(b["type"], "other"), "title": b["title"],
                   "authors": [a.strip() for a in b["author"].split(" and ") if a.strip()][:6],
                   "year": int(b["year"]) if b["year"].isdigit() else None,
                   "identifiers": {"doi": b["doi"], "isbn": "", "pmcid": "", "pmid": "", "url": b["url"]},
                   "locator": "[to be recorded by the verifying author]", "confidentiality": "public"}
        elif key in self.files:
            path = self.files[key]
            code = path.endswith((".R", ".cpp", ".py", ".stan"))
            src = {"source_type": "software" if code else "dataset", "title": f"Project repository file {path}",
                   "authors": [], "year": 2026,
                   "identifiers": {"doi": "", "isbn": "", "pmcid": "", "pmid": "", "url": ""},
                   "locator": path, "confidentiality": "restricted" if key in self.restricted else "public"}
        else:
            raise KeyError(f"unknown evidence key {key}")
        src.update({"evidence_id": e, "verification": {"source_opened": False, "status": "unverified",
                                                       "verified_by": "", "verified_on": ""}})
        self.sources[key] = src
        return e

    def write(self, outdir, docs, extra=None):
        os.makedirs(outdir, exist_ok=True)
        for d in docs:
            open(os.path.join(outdir, d.name + "_audit.md"), "w").write(d.audit_text())
            open(os.path.join(outdir, d.name + ".md"), "w").write(d.clean_text())
        for d in docs:
            with open(os.path.join(outdir, f"claims_{d.name}.csv"), "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["claim_id", "section", "claim_kind", "claim_text_sha256", "evidence_ids",
                            "verification_status", "uncertainty", "analysis_intent"])
                for c in self.claims:
                    if c["section"].startswith(d.name + ":"):
                        w.writerow([c["id"], c["section"], c["kind"], c["sha"], ";".join(c["ev"]), "unverified",
                                    c["unc"], c["intent"]])
        json.dump({"schema_version": "1.0", "sources": sorted(self.sources.values(), key=lambda s: s["evidence_id"])},
                  open(os.path.join(outdir, "source_manifest.json"), "w"), indent=2, ensure_ascii=False)
        for name, obj in (extra or {}).items():
            json.dump(obj, open(os.path.join(outdir, name), "w"), indent=2, ensure_ascii=False)


class Doc:
    """One manuscript file. Each physical line is a paragraph, heading, table row or caption."""

    def __init__(self, reg, name):
        self.reg, self.name, self.lines, self.section = reg, name, [], ""

    def h(self, level, text):
        self.section = text
        self.lines.append(("#" * level + " " + text, None, False))

    def raw(self, text):  # no factual content (front matter, table rules, placeholders)
        if re.search(r"[0-9]", text) and not text.startswith(("---", ":::", "![")):
            ERRORS.append(f"raw line with numbers needs evidence ({self.name}/{self.section}): {text[:70]}")
        self.lines.append((text, None, False))

    def p(self, text, ev=(), kind=None, intent="descriptive", unc="not_applicable", join=False):
        keys = list(ev) + [k for k in CITE_RE.findall(text) if k in self.reg.bib and k not in ev]
        if not keys and not re.search(r"[0-9]", text):
            self.lines.append((text, None, join)); return
        if not keys:
            ERRORS.append(f"numeric line without evidence in {self.name}/{self.section}: {text[:70]}")
            self.lines.append((text, None, join)); return
        eids = []
        for k in keys:
            e = self.reg.eid(k)
            if e not in eids: eids.append(e)
        if kind is None:
            kind = "result" if any(k in self.reg.files for k in keys) else "factual"
        cid = f"C{len(self.reg.claims) + 1:03d}"
        self.reg.claims.append({"id": cid, "section": f"{self.name}: {self.section}"[:120], "kind": kind,
                                "sha": claim_hash(text), "ev": eids, "unc": unc, "intent": intent})
        self.lines.append((text, (cid, eids), join))

    def table(self, header, rows, ev, unc="moderate"):
        self.p("| " + " | ".join(header) + " |", ev, kind="result", unc=unc)
        self.raw("|" + "|".join("---" for _ in header) + "|")
        for r in rows:
            self.p("| " + " | ".join(r) + " |", ev, kind="result", unc=unc)

    def audit_text(self):
        return _join([(t if c is None else f"{t} [claim:{c[0]}] [evidence:{','.join(c[1])}]", j) for t, c, j in self.lines])

    def clean_text(self):
        return _join([(t, j) for t, _, j in self.lines])


def _join(lines):  # blank line between blocks; table rows and joined sentences stay on adjacent lines
    out = []
    for text, join in lines:
        if out and (join or (text.startswith("|") and out[-1].split("\n")[-1].startswith("|"))):
            out[-1] += "\n" + text
        else:
            out.append(text)
    return "\n\n".join(out) + "\n"
