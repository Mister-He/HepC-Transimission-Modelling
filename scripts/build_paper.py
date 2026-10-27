#!/usr/bin/env python3
# =============================================================================
# build_paper.py — builds the PLOS Medicine submission files from the saved outputs:
#   docs/main.docx, docs/Supplementary Information.docx (S1 Text), docs/S1_Checklist.docx,
#   docs/paper/figures/Fig1-4.tif, and the evidence registries in docs/paper/
#   (claims.csv, source_manifest.json, *_audit.md, ...), then runs the scientific-writing
#   skill's local audits if the skill is installed.
# Prose lives in scripts/paper/main_text.py and si_text.py; every number is read from
# output/ and figures/ by scripts/paper/facts.py. Edit there, not in Word.
# Usage: python3 scripts/build_paper.py [ROOT] [--python PY_WITH_PILLOW]
# =============================================================================
import json, os, re, shutil, subprocess, sys, tempfile, zipfile

ROOT = os.path.abspath(next((a for a in sys.argv[1:] if not a.startswith("--")), os.getcwd()))
sys.path.insert(0, os.path.join(ROOT, "scripts", "paper"))
from facts import Facts
from registry import Registry, Doc
import main_text, si_text

DOCS, PAPER = os.path.join(ROOT, "docs"), os.path.join(ROOT, "docs", "paper")
SKILL = os.path.expanduser("~/.claude/skills/scientific-writing/scripts")
MAIN_FIGS = [("Fig1", "output/scenarios_2026/paper_fig1_schematic.png"), ("Fig2", "figures/calibration_agepop.png"),
             ("Fig3", "output/scenarios_2026/paper_fig2_trajectories.png"), ("Fig4", "output/scenarios_2026/paper_fig3_allocation.png")]

CHEERS = [  # item, topic, where it is addressed (item wording to be checked against the official checklist)
    ("1", "Title", "Title page"), ("2", "Abstract", "Abstract"), ("3", "Background and objectives", "Introduction, last paragraph"),
    ("4", "Health economic analysis plan", "Not available: no analysis plan was registered (Methods, Study design)"),
    ("5", "Study population", "Methods, Model structure; S1 Text Section 5"), ("6", "Setting and location", "Methods, Model structure (Singapore)"),
    ("7", "Comparators", "Methods, Policy scenarios"), ("8", "Perspective", "Methods, Costs (charges; payer perspective not available)"),
    ("9", "Time horizon", "Methods, Policy scenarios (2026–2030)"), ("10", "Discount rate", "Methods, Costs; S1 Text Section 7"),
    ("11", "Selection of outcomes", "Methods, Outcomes"), ("12", "Measurement of outcomes", "Methods, Outcomes; S1 Text Section 6"),
    ("13", "Valuation of outcomes", "Not applicable: outcomes were not valued (no utilities)"),
    ("14", "Measurement and valuation of resources and costs", "Methods, Costs; S1 Text Table H"),
    ("15", "Currency, price date, and conversion", "Methods, Costs; S1 Text Section 7"),
    ("16", "Rationale and description of model", "Methods, Model structure; S1 Text Section 1"),
    ("17", "Analytics and assumptions", "Methods; S1 Text Sections 1, 4 and 6"),
    ("18", "Characterising heterogeneity", "Results by group (D, J, F, X) and age; S1 Text Section 5"),
    ("19", "Characterising distributional effects", "Not addressed"),
    ("20", "Characterising uncertainty", "Methods, Uncertainty; S1 Text Section 8 (price sensitivity)"),
    ("21", "Approach to engagement with patients and others affected by the study", "Not addressed: no engagement (to be confirmed by the authors)"),
    ("22", "Study parameters", "S1 Text Tables A, B, D and H"), ("23", "Summary of main results", "Results; Table 1"),
    ("24", "Effect of uncertainty", "Results (uncertainty intervals and probabilities); S1 Text Table I"),
    ("25", "Effect of engagement with patients and others affected by the study", "Not applicable: no engagement"),
    ("26", "Study findings, limitations, generalisability, and current knowledge", "Discussion"),
    ("27", "Source of funding", "Funding statement (to be completed)"), ("28", "Conflicts of interest", "Competing interests (to be completed)"),
]


def pandoc(md, out, title=None):
    args = ["pandoc", md] + (["--metadata", f"title={title}"] if title else []) + [ "--citeproc", "--bibliography", os.path.join(DOCS, "references.bib"),
            "--csl", os.path.join(DOCS, "plos.csl"), "--resource-path", ROOT, "-o", out]
    subprocess.run(args, check=True, cwd=ROOT)
    plos_layout(out)


FOOTER = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?><w:ftr xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main">'
          '<w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:r><w:fldChar w:fldCharType="begin"/></w:r><w:r><w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
          '<w:r><w:fldChar w:fldCharType="separate"/></w:r><w:r><w:t>1</w:t></w:r><w:r><w:fldChar w:fldCharType="end"/></w:r></w:p></w:ftr>')


def plos_layout(path):
    """PLOS layout: double spacing, continuous line numbers and page numbers in the footer."""
    tmp = path + ".tmp"
    with zipfile.ZipFile(path) as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            data = zin.read(item.filename)
            if item.filename == "word/styles.xml":
                x = data.decode()
                x = re.sub(r'w:line="\d+"', 'w:line="480"', x)
                x = x.replace("<w:pPrDefault><w:pPr>", '<w:pPrDefault><w:pPr><w:spacing w:line="480" w:lineRule="auto"/>', 1) \
                    if 'w:line="480"' not in x else x
                data = x.encode()
            elif item.filename == "word/document.xml":
                x = data.decode()
                sect = '<w:footerReference w:type="default" r:id="rIdPlosFooter"/><w:lnNumType w:countBy="1" w:restart="continuous"/>'
                x = re.sub(r"<w:sectPr( [^>]*)?>", lambda m: m.group(0) + sect, x)
                if 'xmlns:r=' not in x[:2000]:
                    x = x.replace("<w:document ", '<w:document xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships" ', 1)
                data = x.encode()
            elif item.filename == "word/_rels/document.xml.rels":
                data = data.decode().replace("</Relationships>", '<Relationship Id="rIdPlosFooter" '
                    'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/footer" Target="footer_plos.xml"/></Relationships>').encode()
            elif item.filename == "[Content_Types].xml":
                data = data.decode().replace("</Types>", '<Override PartName="/word/footer_plos.xml" '
                    'ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.footer+xml"/></Types>').encode()
            zout.writestr(item, data)
        zout.writestr("word/footer_plos.xml", FOOTER)
    os.replace(tmp, path)


def tiff(src, dst, py):  # PLOS wants TIFF figures without captions in the image
    code = "import sys; from PIL import Image; Image.open(sys.argv[1]).convert('RGB').save(sys.argv[2], compression='tiff_lzw', dpi=(300, 300))"
    subprocess.run([py, "-c", code, src, dst], check=True)


def checklist(reg):
    d = Doc(reg, "checklist")
    d.p("Items follow the CHEERS 2022 statement [@cheers2022]; the item wording and the locations below are to be checked by the authors against the official checklist before submission.", kind="declaration")
    d.raw("| Item | Topic | Where addressed |")
    d.raw("|---|---|---|")
    for item, topic, where in CHEERS:
        d.p(f"| {item} | {topic} | {where} |", ("cheers2022",), kind="declaration")
    d.raw("::: {#refs}\n:::")
    return d


def registries(f):
    coverage = {"disclaimer": "Non-scoring coverage aid only; not a quality appraisal or compliance certificate.",
                "guideline_id": "cheers-2022", "schema_version": "1.0", "items": [
        {"topic_id": "population-setting-perspective", "status": "addressed", "locations": ["Methods: Model structure", "Methods: Policy scenarios", "Methods: Costs"]},
        {"topic_id": "effects-costs-and-currency", "status": "addressed", "locations": ["Methods: Outcomes", "Methods: Costs", "S1 Text Section 7"]},
        {"topic_id": "model-assumptions-and-uncertainty", "status": "addressed", "locations": ["Methods: Uncertainty", "S1 Text Sections 1, 4, 6"]},
        {"topic_id": "engagement-results-and-conflicts", "status": "missing", "locations": []}]}
    manifest = {"schema_version": "1.0", "document_id": "hcv-pwid-singapore-plos-medicine", "study_design": "health_economic_evaluation",
                "draft_status": "draft", "submission_ready": False,
                "registries": {"authorship_manifest": "authorship.json", "claim_evidence": "claims_main.csv",
                               "consistency_manifest": "consistency_manifest.json", "reporting_coverage": "reporting_coverage.json",
                               "source_manifest": "source_manifest.json"},
                "reporting_guidelines": ["cheers-2022"],
                "human_verification": {"completed": False, "verified_by": "", "verified_on": ""},
                "confidentiality_review": {"completed": False, "external_services_authorized": False, "policy_checked": False,
                                           "verified_by": "", "verified_on": ""},
                "required_statements": {k: "missing" for k in ("ai_disclosure", "author_contributions", "code_availability", "conflicts",
                                                               "consent", "data_availability", "ethics", "funding")}}
    s = f.summ
    E = lambda *keys: [f._reg.eid(k) for k in keys]
    facts = [("J only: fall in new infections, 2030", "J_only", "incidence_reduction"), ("J only: fall in HCV deaths, 2030", "J_only", "mortality_reduction"),
             ("5,000/year, 25% to F: fall in HCV deaths, 2030", "S5000_F25", "mortality_reduction"),
             ("8,000/year, 25% to F: fall in HCV deaths, 2030", "S8000_F25", "mortality_reduction"),
             ("Uncapped: fall in HCV deaths, 2030", "SInf", "mortality_reduction")]
    nf = [{"fact_id": f"N{i + 1:03d}", "concept": c, "section": "Abstract; Results", "value": round(100 * float(s[k][col]), 1), "unit": "percent",
           "numerator": None, "denominator": None, "sample_size": f.ndraw, "analysis_set": "posterior draws passing the 2015 equilibrium check",
           "evidence_ids": E("summary")} for i, (c, k, col) in enumerate(facts)]
    consistency = {"schema_version": "1.0", "numeric_facts": nf,
                   "methods": [{"method_id": "M001", "name": "Bayesian calibration (adaptive random-walk Metropolis)", "analysis_intent": "descriptive",
                                "protocol_status": "not_applicable", "outcome_ids": ["O001"]},
                               {"method_id": "M002", "name": "Policy scenario projection over posterior draws", "analysis_intent": "exploratory",
                                "protocol_status": "not_applicable", "outcome_ids": ["O002", "O003", "O004"]}],
                   "results": [{"result_id": "R001", "method_id": "M001", "outcome_id": "O001", "analysis_intent": "descriptive", "sample_size": int(float(f.strata["population_D"]["n_draws_used"])),
                                "evidence_ids": E("mcmc_diag", "overlap"), "reported_sections": ["Results: Calibration and the 2015 population"]},
                               {"result_id": "R002", "method_id": "M002", "outcome_id": "O002", "analysis_intent": "exploratory", "sample_size": f.ndraw,
                                "evidence_ids": E("summary"), "reported_sections": ["Abstract", "Results", "Table 1"]},
                               {"result_id": "R003", "method_id": "M002", "outcome_id": "O003", "analysis_intent": "exploratory", "sample_size": f.ndraw,
                                "evidence_ids": E("summary"), "reported_sections": ["Abstract", "Results", "Table 1"]},
                               {"result_id": "R004", "method_id": "M002", "outcome_id": "O004", "analysis_intent": "exploratory", "sample_size": f.ndraw,
                                "evidence_ids": E("summary", "annual"), "reported_sections": ["Results: Costs", "Table 1"]}]}
    authorship = json.load(open(os.path.join(os.path.dirname(SKILL), "assets", "authorship_template.json"))) \
        if os.path.exists(os.path.join(os.path.dirname(SKILL), "assets", "authorship_template.json")) else {}
    if authorship:
        authorship["authors"][0]["name"] = "Yichen He"
        authorship["ai_use"].update({"used": True, "tools": ["Claude (Anthropic)"], "disclosed_in": ["Use of artificial intelligence (draft)"]})
    return {"reporting_coverage.json": coverage, "manuscript_manifest.json": manifest,
            "consistency_manifest.json": consistency, "authorship.json": authorship}


UNRESOLVED = """# Unresolved issues before submission

Generated by scripts/build_paper.py. These items need a decision or a record from the authors;
the manuscript shows each as a bracketed placeholder or states it as a limitation.

1. Authors, affiliations, corresponding author, CRediT roles, funding, competing interests,
   acknowledgments and the AI-use statement are placeholders (authorship.json is a draft).
2. Ethics approval and consent arrangements for the Singapore Prison Service records are not
   recorded in the repository.
3. Data and code availability: repository link/DOI and the access route for the prison records.
4. Protocol registration: none recorded.
5. Every source in source_manifest.json and every claim in claims.csv is "unverified": an author
   must open each source, check the claim, record the locator, and mark it verified.
6. The stored seroprevalence targets differ from positive/screened in the current prison-record
   processing by up to {sero:.1f} percentage points; the provenance of the targets is not recorded.
7. Sources still to be cited: Singapore 2015 age-specific death rates; the base contact matrix;
   the WHO 2030 target wording (80% and 65% falls from 2015) against the WHO documents cited.
8. The CHEERS 2022 item wording in S1 Checklist was written from memory of the statement and must
   be checked against the official checklist; items 19 and 21 are not addressed.
9. The author summary headings follow the usual PLOS Medicine format; check them against the
   current author guidelines (the live guidelines page did not describe them).
10. The neural posterior cross-check failed its calibration check; its policy comparison was run
   before the recovered-DC model fix and is not reported.
"""


def main():
    py = next((sys.argv[i + 1] for i, a in enumerate(sys.argv) if a == "--python"), sys.executable)
    f = Facts(ROOT)
    reg = Registry(ROOT, os.path.join(DOCS, "references.bib"), f.files, restricted=("sps", "prison"))
    f._reg = reg
    docs = [main_text.build(reg, f), si_text.build(reg, f), checklist(reg)]
    import registry
    if registry.ERRORS: sys.exit("\n".join(registry.ERRORS))
    if os.path.isdir(PAPER): shutil.rmtree(PAPER)
    reg.write(PAPER, docs, registries(f))
    open(os.path.join(PAPER, "unresolved_issues.md"), "w").write(UNRESOLVED.format(sero=100 * f.sero_max_diff))
    # main text: figures are uploaded separately (PLOS), captions stay in the text
    pandoc(os.path.join(PAPER, "main.md"), os.path.join(DOCS, "main.docx"))
    pandoc(os.path.join(PAPER, "si.md"), os.path.join(DOCS, "Supplementary Information.docx"))
    pandoc(os.path.join(PAPER, "checklist.md"), os.path.join(DOCS, "S1_Checklist.docx"), title="S1 Checklist. CHEERS 2022 checklist")
    os.makedirs(os.path.join(PAPER, "figures"), exist_ok=True)
    for name, src in MAIN_FIGS:
        tiff(os.path.join(ROOT, src), os.path.join(PAPER, "figures", name + ".tif"), py)
    print(f"Wrote docs/main.docx, docs/Supplementary Information.docx, docs/S1_Checklist.docx; "
          f"{len(reg.claims)} claims, {len(reg.sources)} sources in docs/paper/")
    if os.path.isdir(SKILL):  # local, network-free audits; failures are expected until authors verify sources
        for cmd in (["validate_manifest.py", "source_manifest.json", "--kind", "source"],
                    ["audit_claims.py", "main_audit.md", "claims_main.csv", "source_manifest.json"],
                    ["audit_claims.py", "si_audit.md", "claims_si.csv", "source_manifest.json"],
                    ["audit_claims.py", "checklist_audit.md", "claims_checklist.csv", "source_manifest.json"],
                    ["check_consistency.py", "consistency_manifest.json"], ["check_references.py", "source_manifest.json"],
                    ["select_reporting_guidelines.py", "check", "reporting_coverage.json"],
                    ["validate_manifest.py", "manuscript_manifest.json", "--kind", "manuscript"],
                    ["lint_manuscript.py", "main_audit.md", "--manifest", "manuscript_manifest.json"]):
            r = subprocess.run([sys.executable, os.path.join(SKILL, cmd[0])] + cmd[1:], cwd=PAPER, capture_output=True, text=True)
            try:
                out = json.loads(r.stdout); codes = {}
                for i in out.get("issues", []): codes[i["code"]] = codes.get(i["code"], 0) + 1
                print(f"  {cmd[0]} {cmd[1]}: {out.get('status', '?')} {codes}")
            except json.JSONDecodeError:
                print(f"  {cmd[0]} {cmd[1]}: exit {r.returncode} {r.stdout[:200]} {r.stderr[:200]}")


if __name__ == "__main__":
    main()
