# =============================================================================
# fmt.py — number formatting shared by the main text and S1 Text.
# =============================================================================
F = float
WORDS = "zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen seventeen".split()


def pct(x, d=1):
    return f"{100 * F(x):.{d}f}%"


def ui(m, lo, hi, scale=1.0, d=1, pre="", suf="", first=False):
    f = lambda v: f"{pre}{F(v) * scale:,.{d}f}{suf}"
    tag = "95% uncertainty interval [UI]" if first else "95% UI"
    return f"{f(m)} ({tag} {f(lo)}–{f(hi)})"


def col(row, c, **k):
    return ui(row[c], row[c + "_lo"], row[c + "_hi"], **k)


def red(row, c, d=1, first=False):  # a reduction column as a percentage with its UI
    return col(row, c, scale=100, d=d, suf="%", first=first)


def money(row, c, d=1):
    v = lambda k: f"{F(row[k]) * 1e-6:,.{d}f}"
    return f"S${v(c)} million (95% UI {v(c + '_lo')}–{v(c + '_hi')})"


def n0(x):
    return f"{F(x):,.0f}"


def interval(m, lo, hi, d=0, scale=1.0):
    f = lambda v: f"{F(v) * scale:,.{d}f}"
    return f"{f(m)} ({f(lo)}–{f(hi)})"


def cap_word(n):
    return WORDS[n].capitalize() if n < len(WORDS) else f"{n}"


def scen_label(row):
    if row["scenario_id"] == "J_only": return "J only"
    if row["scenario_id"] == "SInf": return "J + uncapped screening"
    return f"J + {int(F(row['screening_total'])):,}/year, {round(100 * F(row['f_share']))}% to F"
