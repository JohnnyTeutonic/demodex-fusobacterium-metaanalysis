"""Fetch authoritative reference metadata from Crossref for manuscript DOIs,
so bibliography entries are accurate (no hallucinated authors/years/journals)."""
import json
import urllib.request

DOIS = [
    ("liang2021", "10.1007/s40123-021-00356-z"),
    ("yan2020", "10.3389/fmed.2020.592759"),
    ("kim2022", "10.3389/fcimb.2022.922753"),
    ("zou2024", "10.1186/s13071-024-06122-x"),
    ("kraken2", "10.1186/s13059-019-1891-0"),
    ("bracken", "10.7717/peerj-cs.104"),
    ("whitlock2005", "10.1111/j.1420-9101.2005.00917.x"),
    ("bh1995", "10.1111/j.2517-6161.1995.tb02031.x"),
    ("quast2013", "10.1093/nar/gks1219"),
    ("callahan2016", "10.1038/nmeth.3869"),
    ("martin2011", "10.14806/ej.17.1.200"),
    ("callahan2017", "10.1038/ismej.2017.119"),
    ("chang2017", "10.1016/j.jaad.2017.03.040"),
]

def fmt(key, doi):
    url = f"https://api.crossref.org/works/{doi}?mailto=jonathanreich100@gmail.com"
    try:
        with urllib.request.urlopen(url, timeout=40) as r:
            m = json.load(r)["message"]
    except Exception as e:
        return f"{key}\tERROR {e}"
    auth = m.get("author", []) or []
    names = ", ".join(
        f"{a.get('family','')} {''.join(p[0] for p in a.get('given','').split())}".strip()
        for a in auth[:8]
    )
    if len(auth) > 8:
        names += ", et al"
    title = (m.get("title") or [""])[0]
    ctr = (m.get("container-title") or [""])[0]
    yr = None
    for k in ("published-print", "published-online", "issued"):
        if m.get(k, {}).get("date-parts"):
            yr = m[k]["date-parts"][0][0]
            break
    vol = m.get("volume", "")
    iss = m.get("issue", "")
    pg = m.get("page", "")
    return f"{key}\t{names}. {title}. {ctr}. {yr};{vol}({iss}):{pg}. doi:{doi}"

for key, doi in DOIS:
    print(fmt(key, doi))
    print()
