#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V50 - a montagem oficial do RCSB tem o mesmo raio que a nossa?

PERGUNTA. O 5NEU declara 60 rotacoes puras em torno da origem (t = 0 em todos os
operadores) e tem as coordenadas centradas a 601 A da origem, o que poe a
expansao numa casca de 546 a 666 A. O RCSB exibe para ele uma Biological
Assembly 1 icosaedrica hetero-360-merica, e a nossa expansao produz 360 cadeias:
a contagem concorda, o raio e que esta em questao.

Restam duas leituras e este script decide entre elas:

  (a) o deposito e assim        -> a montagem oficial tambem tem ~600 A de raio,
                                   nossa leitura esta correta e a exclusao se
                                   sustenta pelo que o arquivo contem
  (b) o RCSB aplica algo a mais -> a montagem oficial tem ~165 A, e entao o nosso
                                   leitor perde alguma coisa e o 5NEU NAO deve
                                   ser excluido

CONTROLES. O 3J9F e o controle que importa, porque compartilha a convencao do
5NEU: t = 0 em todos os sessenta operadores. O 6I2K usa a convencao com centro
fora da origem, |t| ate 935 A. Os dois reconstroem capsideos corretos pela nossa
rota, entao se ambos reproduzirem a montagem oficial e so o 5NEU divergir, a
leitura (a) esta estabelecida.

Baixa os arquivos de montagem para um diretorio NOVO e nao toca em nada do
projeto. Se o download falhar, diz e nao adivinha.

  pip install gemmi
  python3 V50_assembleia_oficial.py
"""
import os, sys, math, gzip, shutil, urllib.request, urllib.error

BASE = "/home/juan"
ALVOS = [("5NEU", "em exame"), ("3J9F", "controle, t = 0 como o 5NEU"),
         ("6I2K", "controle, centro fora da origem")]
LOCAIS = [f"{BASE}/cif_proveniencia", f"{BASE}/figuras_V28", f"{BASE}/cif_criterio",
          f"{BASE}/abertos"]

try:
    import gemmi
except ImportError:
    sys.exit("falta gemmi:  pip install gemmi")


def novo_dir(p):
    if not os.path.exists(p):
        return p
    n = 2
    while os.path.exists(f"{p}_{n:02d}"):
        n += 1
    return f"{p}_{n:02d}"


DST = novo_dir(f"{BASE}/V50_assembleias_rcsb")
os.makedirs(DST)
print(f"destino dos downloads: {DST}\n")

URLS = ["https://files.rcsb.org/download/{p}-assembly1.cif.gz",
        "https://files.rcsb.org/download/{p}-assembly1.cif",
        "https://models.rcsb.org/{p}-assembly1.cif"]


def baixar(pdb):
    for u in URLS:
        url = u.format(p=pdb.lower() if "files.rcsb" in u else pdb)
        alvo = os.path.join(DST, os.path.basename(url))
        try:
            with urllib.request.urlopen(url, timeout=90) as r, open(alvo, "wb") as f:
                shutil.copyfileobj(r, f)
        except Exception as e:
            print(f"    {url}  -> {type(e).__name__}: {e}")
            continue
        if alvo.endswith(".gz"):
            fora = alvo[:-3]
            with gzip.open(alvo, "rb") as a, open(fora, "wb") as b:
                shutil.copyfileobj(a, b)
            alvo = fora
        return alvo, url
    return None, None


def achar_local(pdb):
    for d in LOCAIS:
        for ext in (".cif", ".cif.gz", ".pdb"):
            p = os.path.join(d, f"{pdb}{ext}")
            if os.path.exists(p):
                return p
    return None


def ca(st):
    return [(a.pos.x, a.pos.y, a.pos.z) for c in st[0] for r in c for a in r
            if a.name == "CA"]


def medidas(st):
    P = ca(st)
    if not P:
        return None
    n = len(P)
    cx = sum(p[0] for p in P) / n
    cy = sum(p[1] for p in P) / n
    cz = sum(p[2] for p in P) / n
    d = sorted(math.dist(p, (cx, cy, cz)) for p in P)
    return {"cadeias": len(st[0]), "ca": n, "min": d[0],
            "med": d[n // 2], "max": d[-1],
            "centroide": (cx, cy, cz),
            "dist_origem": math.dist((cx, cy, cz), (0.0, 0.0, 0.0))}


def linha(rot, m):
    if not m:
        print(f"  {rot:<28} sem CA")
        return
    print(f"  {rot:<28}{m['cadeias']:>5} cad {m['ca']:>7} CA   "
          f"raio {m['min']:6.1f} {m['med']:6.1f} {m['max']:6.1f} A   "
          f"|centroide| {m['dist_origem']:6.1f} A")


resumo = []
for pdb, papel in ALVOS:
    print("=" * 74)
    print(f"{pdb}   [{papel}]")
    print("=" * 74)

    # 1. a montagem oficial
    print("  baixando a montagem oficial...")
    arq, url = baixar(pdb)
    m_of = None
    if arq:
        print(f"  de {url}")
        try:
            st = gemmi.read_structure(arq)
            st.setup_entities()
            m_of = medidas(st)
        except Exception as e:
            print(f"  ilegivel: {e}")
    else:
        print("  download indisponivel; sem a montagem oficial nao ha comparacao")

    # 2. a nossa expansao, do arquivo local
    loc = achar_local(pdb)
    m_nos = None
    if loc:
        try:
            s2 = gemmi.read_structure(loc)
            s2.setup_entities()
            nome = max(s2.assemblies,
                       key=lambda a: sum(len(g.operators) for g in a.generators)).name
            s2.transform_to_assembly(nome, gemmi.HowToNameCopiedChain.AddNumber)
            m_nos = medidas(s2)
        except Exception as e:
            print(f"  expansao local falhou: {e}")
    else:
        print("  arquivo local nao encontrado")

    linha("montagem oficial RCSB", m_of)
    linha("nossa expansao", m_nos)
    if m_of and m_nos:
        dmax = abs(m_of["max"] - m_nos["max"])
        igual = dmax < 5.0 and abs(m_of["ca"] - m_nos["ca"]) <= max(1, 0.02 * m_of["ca"])
        print(f"  -> {'COINCIDEM' if igual else 'DIVERGEM'}  "
              f"(diferenca no raio maximo: {dmax:.1f} A)")
        resumo.append((pdb, igual, m_of["max"], m_nos["max"]))
    print()

print("=" * 74)
if len(resumo) == len(ALVOS):
    so_5neu_grande = all(
        (m_of > 400) == (p == "5NEU") for p, _, m_of, _ in resumo)
    todos_coincidem = all(ok for _, ok, _, _ in resumo)
    print("leitura dos resultados:")
    for p, ok, a, b in resumo:
        print(f"  {p:<6} oficial {a:6.1f} A   nossa {b:6.1f} A   "
              f"{'coincidem' if ok else 'DIVERGEM'}")
    print()
    if todos_coincidem and so_5neu_grande:
        print("  (a) A montagem oficial do 5NEU tem o mesmo raio que a nossa, e os")
        print("      dois controles tambem coincidem. Nossa leitura esta correta e a")
        print("      exclusao do 5NEU se sustenta pelo que o deposito contem.")
    elif todos_coincidem:
        print("  As tres coincidem, mas o padrao de raios nao e o esperado.")
        print("  Olhar os numeros caso a caso antes de concluir.")
    else:
        div = [p for p, ok, _, _ in resumo if not ok]
        print(f"  (b) Divergencia em {div}. Se o 5NEU divergir e os controles nao,")
        print("      o RCSB aplica algo que o nosso leitor nao aplica, e o 5NEU NAO")
        print("      deve ser excluido ate isso ser entendido.")
else:
    print("comparacao incompleta; sem os tres pares nao da para decidir.")
print()
print("Nada do projeto foi alterado. Os downloads estao em")
print(f"  {DST}")
