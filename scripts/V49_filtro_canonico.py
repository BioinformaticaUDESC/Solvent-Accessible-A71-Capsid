#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V49 - o filtro de registro, como o paper vai descrever.

DECISOES TOMADAS (24/11/2026), que este script implementa e nada mais:

  regra        bloqueio <= piso          UM eixo. O eixo radial sai do filtro
                                         e fica como coluna descritiva.
  piso         percentil 95 do bloqueio nos residuos de pegada da metade de
               calibracao, estimador np.quantile, alvo de sensibilidade
               declarado antes da medida
  divisao      ordem do validacao_filtro.tsv, semente 20261116, metade a metade
               (convencao do piso_calibra_testa.py; deixa o 6I2K retido)
  excluidos    23BG e 23PS - o deposito nao fornece os operadores necessarios
               para reconstruir a particula

Desempenho medido no V48l, divisao A, metade retida: sensibilidade 96,5 %,
especificidade 65,3 %. Este script recalcula tudo em vez de repetir numero.

O QUE ELE ESCREVE. A coluna `passes` das pegadas atuais foi gerada pelo
filtro_externo.py com `shell_pct >= 70 and blocking < 50`, que nao e o filtro
publicado. Este script recomputa a decisao a partir das colunas que ja existem
-- nao refaz medida estrutural nenhuma -- e grava COPIAS NOVAS num diretorio
novo, com a coluna `passes_V49` ao lado da original. Nenhum arquivo teu e
alterado ou sobrescrito.

Tambem imprime o veredito para o residuo 22 do alvo, lido do perfil radial.

  python3 V49_filtro_canonico.py
"""
import os, re, sys, csv, glob, random, shutil, datetime
import numpy as np

BASE = "/home/juan"
RAIZES = [f"{BASE}/abertos", f"{BASE}/paper_artifacts", f"{BASE}/figuras_V28", BASE]
SEMENTE = 20261116
SENS = 0.95
FORA = {"23BG", "23PS"}
SAIDA = f"{BASE}/V49_pegadas"
ALVO_RES = "22"


def novo_dir(p):
    if not os.path.exists(p):
        return p
    n = 2
    while os.path.exists(f"{p}_{n:02d}"):
        n += 1
    return f"{p}_{n:02d}"


def novo(p):
    if not os.path.exists(p):
        return p
    raiz, ext = os.path.splitext(p)
    n = 2
    while os.path.exists(f"{raiz}_{n:02d}{ext}"):
        n += 1
    return f"{raiz}_{n:02d}{ext}"


def le(p):
    with open(p, encoding="utf8", errors="replace") as f:
        ln = f.read().splitlines()
    if not ln:
        return [], []
    cab = ln[0].split("\t")
    return cab, [dict(zip(cab, l.split("\t"))) for l in ln[1:] if l.strip()]


def col(cs, ch):
    for c in cs:
        if c and c.strip().lower() in ch:
            return c
    for c in cs:
        if c and any(k in c.strip().lower() for k in ch):
            return c
    return None


def f(v):
    try:
        return float(str(v).replace(",", ".").strip())
    except Exception:
        return None


def dir_trabalho():
    for r in RAIZES:
        g = glob.glob(os.path.join(r, "**", "validacao_filtro.tsv"), recursive=True)
        if g:
            return os.path.dirname(sorted(g, key=len)[0])
    return None


D = dir_trabalho()
if not D:
    sys.exit("validacao_filtro.tsv nao localizado")
print(f"diretorio de trabalho: {D}")

# ------------------------------------------------------------- a divisao
cab_v, rows_v = le(os.path.join(D, "validacao_filtro.tsv"))
comp = []
for d in rows_v:
    if d.get("pdb") in FORA or not d.get("auc_casca"):
        continue
    p = os.path.join(D, f"{d.get('pdb')}_{d.get('cadeia')}_pegada.tsv")
    if os.path.exists(p):
        comp.append((d["pdb"], p))
random.Random(SEMENTE).shuffle(comp)
meio = len(comp) // 2
cal, tes = comp[:meio], comp[meio:]
print(f"\ndivisao: semente {SEMENTE}, ordem do validacao_filtro.tsv")
print(f"  calibracao ({len(cal)}): {' '.join(p for p, _ in cal)}")
print(f"  teste      ({len(tes)}): {' '.join(p for p, _ in tes)}")
print(f"  excluidos: {' '.join(sorted(FORA))}  (deposito sem operadores)")


def carrega(path):
    cs, rows = le(path)
    cb = col(cs, {"blocking", "bloqueio", "blocking_atoms"})
    cs_ = col(cs, {"shell_pct", "casca_pct", "percentil_radial"})
    cf = col(cs, {"in_footprint", "pegada", "contato"})
    out = []
    for r in rows:
        b = f(r.get(cb)) if cb else None
        s = f(r.get(cs_)) if cs_ else None
        if b is None:
            continue
        peg = str(r.get(cf, "")).strip().lower() in ("yes", "sim", "true", "1") if cf else None
        out.append((s, b, peg))
    return cs, rows, out, cb, cf


def junta(grupo):
    blo, peg = [], []
    for _, p in grupo:
        _, _, dados, _, _ = carrega(p)
        for s, b, pg in dados:
            if pg is None:
                continue
            blo.append(b)
            peg.append(pg)
    return np.array(blo), np.array(peg, bool)


bc, pc = junta(cal)
bt, pt = junta(tes)
PISO = float(np.quantile(bc[pc], SENS))
print(f"\nresiduos: calibracao {pc.sum()} de pegada em {len(pc)}   "
      f"teste {pt.sum()} de {len(pt)}")
print(f"\nPISO PUBLICADO:  bloqueio <= {PISO:.0f}")
print(f"  percentil {SENS:.0%} do bloqueio nos residuos de pegada da calibracao")
print(f"  estimador np.quantile; operador <= ; um eixo")


def mede(b, p, piso):
    passa = b <= piso
    return (float(passa[p].mean()) if p.any() else float("nan"),
            float((~passa[~p]).mean()) if (~p).any() else float("nan"))


sc, ec = mede(bc, pc, PISO)
st, et = mede(bt, pt, PISO)
print(f"\ndesempenho")
print(f"  calibracao     sensibilidade {sc:6.1%}   especificidade {ec:6.1%}")
print(f"  TESTE RETIDO   sensibilidade {st:6.1%}   especificidade {et:6.1%}")

# ------------------------------------------- recomputa passes em copias novas
dst = novo_dir(SAIDA)
os.makedirs(dst)
todas = set()
for r in RAIZES:
    for p in glob.glob(os.path.join(r, "**", "*_pegada.tsv"), recursive=True):
        todas.add(os.path.realpath(p))
em_cal = {p for p, _ in cal}
em_tes = {p for p, _ in tes}
manifesto, n_arq = [], 0
for p in sorted(todas):
    base = os.path.basename(p)
    pdb = base.split("_")[0].upper()
    cs, rows, dados, cb, cf = carrega(p)
    if not cb or not rows:
        continue
    cpass = col(cs, {"passes", "passa"})
    novo_cab = cs + (["passes_V49"] if "passes_V49" not in cs else [])
    saida_p = os.path.join(dst, base.replace("_pegada.tsv", "_pegada_V49.tsv"))
    if os.path.exists(saida_p):            # mesmo nome em raizes diferentes
        saida_p = novo(saida_p)
    n_sim = n_nao = 0
    conc = disc = 0
    with open(saida_p, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=novo_cab, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            b = f(r.get(cb))
            v = "" if b is None else ("yes" if b <= PISO else "no")
            r = dict(r)
            r["passes_V49"] = v
            if v == "yes":
                n_sim += 1
            elif v == "no":
                n_nao += 1
            if cpass and v:
                antigo = str(r.get(cpass, "")).strip().lower() in ("yes", "sim", "true", "1")
                if antigo == (v == "yes"):
                    conc += 1
                else:
                    disc += 1
            w.writerow(r)
    n_arq += 1
    manifesto.append({"pdb": pdb, "arquivo_origem": p,
                      "arquivo_V49": saida_p,
                      "particao": "calibracao" if p in em_cal else
                                  ("teste" if p in em_tes else
                                   ("excluido" if pdb in FORA else "fora do conjunto")),
                      "n_residuos": len(rows), "passa_V49": n_sim,
                      "reprova_V49": n_nao,
                      "concorda_com_passes_antigo": conc,
                      "discorda": disc})

print(f"\n{n_arq} pegadas recomputadas em copias novas")
print(f"  -> {dst}/")
print("  originais intactos; coluna passes_V49 ao lado da antiga")
d_tot = sum(m["discorda"] for m in manifesto)
c_tot = sum(m["concorda_com_passes_antigo"] for m in manifesto)
if c_tot + d_tot:
    print(f"  contra a coluna passes antiga (70; 50): {d_tot} divergencias em "
          f"{c_tot + d_tot} residuos ({100*d_tot/(c_tot+d_tot):.1f} %)")

man = novo(f"{BASE}/V49_FILTRO_CANONICO.tsv")
with open(man, "w", newline="", encoding="utf8") as fh:
    campos = ["pdb", "particao", "n_residuos", "passa_V49", "reprova_V49",
              "concorda_com_passes_antigo", "discorda", "arquivo_V49", "arquivo_origem"]
    w = csv.DictWriter(fh, fieldnames=campos, delimiter="\t", extrasaction="ignore")
    w.writeheader()
    for m in manifesto:
        w.writerow(m)

reg = novo(f"{BASE}/V49_METODO.txt")
with open(reg, "w", encoding="utf8") as fh:
    fh.write(f"""V49 - filtro de registro
gerado em {datetime.datetime.now().isoformat(timespec='seconds')}

regra            bloqueio <= {PISO:.0f}   (um eixo)
piso             percentil {SENS:.0%} do bloqueio nos residuos de pegada da
                 metade de calibracao, estimador np.quantile
divisao          semente {SEMENTE}, ordem do validacao_filtro.tsv, metade a metade
calibracao       {' '.join(p for p, _ in cal)}
teste            {' '.join(p for p, _ in tes)}
excluidos        {' '.join(sorted(FORA))} - deposito sem os operadores de simetria
eixo radial      shell_pct permanece como coluna descritiva, fora do filtro

desempenho       calibracao   sens {sc:.1%}  espec {ec:.1%}
                 teste retido sens {st:.1%}  espec {et:.1%}
""")

# -------------------------------------------------------------- o alvo
print("\nresiduo 22 do alvo, lido do perfil radial")
achou = False
perfis = set()
for r in RAIZES:
    perfis |= {os.path.realpath(p)
               for p in glob.glob(os.path.join(r, "**", "*perfil_radial*"), recursive=True)
               if p.endswith((".tsv", ".csv"))}
for p in sorted(perfis):
    if True:
        cs, rows = le(p)
        cb = col(cs, {"blocking_atoms", "blocking", "bloqueio"})
        cs_ = col(cs, {"shell_pct", "casca_pct"})
        if not cb:
            continue
        for row in rows:
            if str(row.get("residue", "")).strip() != ALVO_RES:
                continue
            b, s = f(row.get(cb)), f(row.get(cs_)) if cs_ else None
            if b is None:
                continue
            achou = True
            print(f"  {os.path.basename(p):<34} {row.get('resname','')}  "
                  f"bloqueio {b:5.0f}  casca {s if s is None else f'{s:.2f}'}  "
                  f"-> {'passa' if b <= PISO else 'REPROVA'}  (piso {PISO:.0f})")
if not achou:
    print("  nenhum perfil radial com residuo 22 localizado")

print(f"\n-> {man}")
print(f"-> {reg}")
