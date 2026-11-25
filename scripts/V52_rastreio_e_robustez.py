#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V52 - rastreia os valores marcados e mede robustez por complexo.

PARTE A - os onze [CHECK]. Cada valor que o manuscrito traz marcado vem de uma
rodada antiga que eu nao reconferi. Aqui cada um e procurado nos arquivos do
projeto com padrao proprio E palavra de contexto obrigatoria, porque buscar
inteiro nu ja falhou: no V48k eu varri 600 mil arquivos e colhi 53 mil
ocorrencias de '112', quase todas em 'print("=" * 112)'. As raizes aqui sao as
do projeto e nao a home inteira.

Para cada acerto: arquivo, linha, o texto da linha e o SHA-256 do arquivo.

PARTE B - robustez por complexo. As proporcoes de 96,5 % e 65,3 % sao contagens
por residuo sobre treze complexos, e residuos do mesmo complexo nao sao
observacoes independentes. Aqui a metade retida e reavaliada deixando um
complexo de fora por vez, com o limiar FIXO, mais a sensibilidade de cada
complexo isoladamente. Mostra se a estimativa depende de um complexo so.

NAO recalibra o limiar. NAO roda simulacao. Somente leitura.

  python3 V52_rastreio_e_robustez.py
"""
import os, re, sys, csv, glob, math, hashlib
import numpy as np

BASE = "/home/juan"
RAIZES = [f"{BASE}/abertos", f"{BASE}/paper_artifacts", f"{BASE}/figuras_V28"]
RAIZES += sorted(glob.glob(f"{BASE}/entrega_virus_snca_ppdock*"))
SOLTOS = sorted(glob.glob(f"{BASE}/*.py") + glob.glob(f"{BASE}/*.sh")
                + glob.glob(f"{BASE}/*.tsv") + glob.glob(f"{BASE}/*.json")
                + glob.glob(f"{BASE}/*.md"))
METODO = f"{BASE}/V49_METODO.txt"
PISO = 106.0
EXTS = (".py", ".sh", ".tsv", ".csv", ".json", ".md", ".txt", ".log", ".ipynb")

# (rotulo, regex do valor, palavras de contexto obrigatorias na linha)
ALVOS = [
    ("retencao de SASA 91,2 %", r"(?<![\d.])(0\.912|91\.2)(?![\d])",
     r"sasa|reten|exposed|area|acess"),
    ("rho entre descritores -0,975", r"-\s?0\.975", r"spearman|rho|corr|casca|shell|bloq"),
    ("rho SASA x posicao radial +0,329", r"\+?0\.329", r"spearman|rho|corr|sasa|radial|casca"),
    ("sitio 22 em 37 de 100 poses", r"(?<![\d.])37(?![\d.])",
     r"pose|solucao|solution|top|contato|contact|100"),
    ("variacao de replicas 3,832", r"3\.832", r"replic|semente|seed|desvio|median|variac"),
    ("resolucao do 3VBS", r"(?<![\d.])(3\.0|3\.00)(?![\d])", r"3vbs|resolu"),
    ("atomos da montagem do 3VBS", r"391[,.\s]?740", r"atom|montagem|assembl|expand"),
    ("inventario de 29 a 30 linhas", r"(?<![\d.])(29|30)(?![\d.])",
     r"identifi|linhas|rows|inventar|validacao_filtro"),
    ("reamostragem de 200 divisoes", r"(?<![\d.])200(?![\d.])",
     r"split|divis|repart|range|rep\b|shuffle|boot"),
    ("versao do LightDock", r"lightdock", r"version|versao|0\.9|==|pip|setup"),
]


def novo(p):
    if not os.path.exists(p):
        return p
    r, e = os.path.splitext(p)
    n = 2
    while os.path.exists(f"{r}_{n:02d}{e}"):
        n += 1
    return f"{r}_{n:02d}{e}"


def ler(p):
    with open(p, newline="", encoding="utf8", errors="replace") as f:
        h = f.read(4096)
        f.seek(0)
        return list(csv.DictReader(f, delimiter="\t" if h.count("\t") >= h.count(",") else ","))


def num(v):
    try:
        return float(str(v).replace(",", ".").strip())
    except Exception:
        return float("nan")


def col(cs, ch):
    for c in cs:
        if c and c.strip().lower() in ch:
            return c
    for c in cs:
        if c and any(k in c.strip().lower() for k in ch):
            return c
    return None


def verdadeiro(v):
    return str(v).strip().lower() in ("1", "true", "sim", "yes", "t", "y")


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# =========================================================== PARTE A
print("=" * 78)
print("(A) rastreio dos valores marcados")
print("=" * 78)
arquivos = set()
for r in RAIZES:
    for p in glob.glob(os.path.join(r, "**", "*"), recursive=True):
        if p.endswith(EXTS) and os.path.isfile(p):
            arquivos.add(os.path.realpath(p))
arquivos |= {os.path.realpath(p) for p in SOLTOS}
arquivos.discard(os.path.realpath(__file__))     # senao casa com os proprios regex
arquivos = {p for p in arquivos if os.path.getsize(p) < 30_000_000}
print(f"{len(arquivos)} arquivos nas raizes do projeto\n")

hashes, achados = {}, []
for rot, rx_v, rx_c in ALVOS:
    rv, rc = re.compile(rx_v, re.I), re.compile(rx_c, re.I)
    hits = []
    for p in sorted(arquivos):
        try:
            with open(p, encoding="utf8", errors="replace") as f:
                for i, ln in enumerate(f, 1):
                    if len(ln) > 3000:
                        continue
                    if rv.search(ln) and rc.search(ln):
                        hits.append((p, i, ln.rstrip()[:118]))
                        if len(hits) >= 12:
                            break
        except Exception:
            continue
        if len(hits) >= 12:
            break
    print(f"--- {rot}: {len(hits)} acerto(s)")
    for p, i, ln in hits[:6]:
        if p not in hashes:
            hashes[p] = sha(p)
        aud = "  [auditoria]" if re.match(r"V4[89]|V5\d", os.path.basename(p)) else ""
        print(f"  {os.path.relpath(p, BASE)}:{i}{aud}")
        print(f"      {ln}")
        achados.append({"valor": rot, "arquivo": p, "linha": i, "texto": ln,
                        "sha256": hashes[p],
                        "origem": "script de auditoria" if aud else "projeto"})
    if len(hits) > 6:
        print(f"  ... e mais {len(hits)-6}")
    if not hits:
        print("  nenhum acerto com contexto -> nao localizado; retirar a alegacao")
        achados.append({"valor": rot, "arquivo": "", "linha": "", "texto": "",
                        "sha256": "", "origem": "NAO LOCALIZADO"})
    print()

# =========================================================== PARTE B
print("=" * 78)
print(f"(B) robustez por complexo na metade retida, piso fixo em {PISO:.0f}")
print("=" * 78)
cal = tes = []
if os.path.exists(METODO):
    txt = open(METODO, encoding="utf8", errors="replace").read()
    m = re.search(r"calibracao\s+([0-9A-Z \t]+)\n", txt)
    n_ = re.search(r"teste\s+([0-9A-Z \t]+)\n", txt)
    cal = m.group(1).split() if m else []
    tes = n_.group(1).split() if n_ else []
if not tes:
    sys.exit(f"nao consegui ler as metades de {METODO}; rode o V49 antes")

peg, vistos = {}, set()
for r in RAIZES + [f"{BASE}"]:
    for p in glob.glob(os.path.join(r, "**", "*_pegada.tsv"), recursive=True):
        rp = os.path.realpath(p)
        if rp in vistos:
            continue
        vistos.add(rp)
        pdb = os.path.basename(p).split("_")[0].upper()
        if pdb not in tes or pdb in peg:
            continue
        rows = ler(rp)
        if not rows:
            continue
        cs = list(rows[0].keys())
        cb = col(cs, {"blocking", "bloqueio", "blocking_atoms"})
        cf = col(cs, {"in_footprint", "pegada", "contato"})
        if not (cb and cf):
            continue
        c, n = [], []
        for row in rows:
            b = num(row.get(cb))
            if math.isnan(b):
                continue
            (c if verdadeiro(row.get(cf)) else n).append(b)
        peg[pdb] = (c, n, rp)

falta = [p for p in tes if p not in peg]
print(f"metade retida: {len(peg)} de {len(tes)} com pegada"
      + (f"   FALTA: {falta}" if falta else ""))


def metricas(pdbs):
    c = [x for p in pdbs if p in peg for x in peg[p][0]]
    n = [x for p in pdbs if p in peg for x in peg[p][1]]
    if not c or not n:
        return None
    c, n = np.array(c), np.array(n)
    return (float((c <= PISO).mean()), float((n > PISO).mean()), len(c), len(n))


base = metricas(tes)
print(f"\ncompleta         sensibilidade {base[0]:6.1%}   especificidade {base[1]:6.1%}"
      f"   ({base[2]} contato, {base[3]} nao)")
print(f"\ndeixando um complexo de fora por vez:")
print(f"{'fora':<8}{'sens':>8}{'espec':>8}{'n cont':>8}   e so esse complexo:")
linhas, sens_l, esp_l = [], [], []
for fora in tes:
    if fora not in peg:
        continue
    m = metricas([p for p in tes if p != fora])
    so = metricas([fora])
    if not m:
        continue
    sens_l.append(m[0])
    esp_l.append(m[1])
    s_txt = f"{so[0]:6.1%} / {so[1]:6.1%}" if so else "   sem par"
    print(f"{fora:<8}{m[0]:>8.1%}{m[1]:>8.1%}{m[2]:>8}   {s_txt}")
    linhas.append({"deixado_fora": fora, "sens_sem_ele": round(m[0]*100, 1),
                   "espec_sem_ele": round(m[1]*100, 1), "n_contato_restante": m[2],
                   "sens_so_dele": round(so[0]*100, 1) if so else "",
                   "espec_so_dele": round(so[1]*100, 1) if so else "",
                   "n_contato_dele": len(peg[fora][0]),
                   "arquivo": peg[fora][2]})
if sens_l:
    print(f"\namplitude deixando um de fora:")
    print(f"  sensibilidade  {min(sens_l):.1%} a {max(sens_l):.1%}"
          f"   (completa {base[0]:.1%})")
    print(f"  especificidade {min(esp_l):.1%} a {max(esp_l):.1%}"
          f"   (completa {base[1]:.1%})")
    print("\n  Uma amplitude estreita indica que a estimativa nao depende de um")
    print("  complexo so. Nao e intervalo de confianca e nao deve ser citada como")
    print("  um: e analise de influencia, nao inferencia.")

s1 = novo(f"{BASE}/V52_RASTREIO_CHECK.tsv")
with open(s1, "w", newline="", encoding="utf8") as f:
    w = csv.DictWriter(f, fieldnames=["valor", "origem", "arquivo", "linha",
                                      "texto", "sha256"], delimiter="\t",
                       extrasaction="ignore")
    w.writeheader()
    for r in achados:
        w.writerow(r)
if linhas:
    s2 = novo(f"{BASE}/V52_ROBUSTEZ_COMPLEXO.tsv")
    with open(s2, "w", newline="", encoding="utf8") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()), delimiter="\t")
        w.writeheader()
        for r in linhas:
            w.writerow(r)
    print(f"\n-> {s2}")
print(f"-> {s1}")
print("\nOs acertos marcados [auditoria] vem dos meus proprios scripts e nao")
print("servem de fonte primaria: um valor que so aparece la nao esta verificado.")
