#!/usr/bin/env python3
"""Valida o filtro de acessibilidade com a pegada MEDIDA de um complexo.

Ordem correta das operacoes:
  1. expande o complexo INTEIRO (capsideo + receptor) por BIOMT/MTRIX
  2. mede a pegada da cadeia alvo contra TODAS as copias do receptor
     (o canyon fica entre protomeros vizinhos: na unidade assimetrica so
      aparece um fragmento da pegada)
  3. remove todas as copias do receptor
  4. calcula casca radial e bloqueio sobre o capsideo montado

Uso: python3 filtro_externo.py 3J8F --cadeia 1
     python3 filtro_externo.py 6I2K --cadeia A
"""
import argparse, sys, urllib.request
from pathlib import Path
import numpy as np, gemmi

AA = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER "
         "THR TRP TYR VAL MSE".split())
RAIO, CONE, EXT, MAXB = 30.0, 0.7071, 70.0, 50
REC = ("RECEPTOR", "CD155", "PVR", "ICAM", "SCARB2", "LIMP", "LYSOSOM",
       "PSGL", "NECTIN", "CD55", "SELECTIN", "HEPARAN", "INTEGRIN",
       "SIALOADHESIN", "DC-SIGN", "ANNEXIN", "KREMEN", "FCGRT", "CAR ")
CAP = ("VP1", "VP2", "VP3", "VP4", "VP0", "CAPSID", "POLYPROTEIN", "COAT")


def compnd(p):
    m, cur = {}, None
    for ln in p.read_text(errors="ignore").splitlines():
        if not ln.startswith("COMPND"):
            if ln.startswith(("ATOM", "HETATM", "SEQRES")):
                break
            continue
        b = ln[10:].strip().rstrip(";")
        if b.upper().startswith("MOLECULE:"):
            cur = b.split(":", 1)[1].strip()
        elif b.upper().startswith("CHAIN:") and cur:
            for c in b.split(":", 1)[1].split(","):
                m[c.strip()] = cur
    return m


def auc(a, b):
    a, b = np.asarray(a, float)[:, None], np.asarray(b, float)[None, :]
    return float(((a > b).sum() + .5 * (a == b).sum()) / a.size / b.size)


def xyz(ch):
    v = [[a.pos.x, a.pos.y, a.pos.z] for r in ch for a in r]
    return (np.array(v, dtype=np.float32) if v
            else np.zeros((0, 3), dtype=np.float32))


def expande(st):
    """Expande o complexo inteiro, devolvendo [(cadeia, nome_original)].

    Constroi a montagem a mao a partir dos geradores para nunca perder o
    rastro de qual cadeia original gerou cada copia: sem isso nao da para
    remover o receptor depois de expandir.
    """
    orig = {c.name: c for c in st[0]}
    if st.assemblies:
        for asm in st.assemblies:
            saida, n_ops = [], 0
            for gen in asm.generators:
                alvos = [c for c in gen.chains if c in orig] or list(orig)
                for op in gen.operators:
                    n_ops += 1
                    for nome in alvos:
                        c = orig[nome].clone()
                        for r in c:
                            for a in r:
                                a.pos = gemmi.Position(op.transform.apply(a.pos))
                        saida.append((c, nome))
            if len(saida) > len(orig):
                return saida, f"BIOMT '{asm.name}' ({n_ops} ops)"
    ops = [o for o in st.ncs if not o.given]
    if ops:
        saida = [(orig[n].clone(), n) for n in orig]
        for op in ops:
            for nome in orig:
                c = orig[nome].clone()
                for r in c:
                    for a in r:
                        a.pos = gemmi.Position(op.tr.apply(a.pos))
                saida.append((c, nome))
        return saida, f"MTRIX ({len(ops)} ops)"
    return [(orig[n], n) for n in orig], None


ap = argparse.ArgumentParser()
ap.add_argument("entrada")
ap.add_argument("--cadeia", help="cadeia do capsideo a analisar (ex: 1, A)")
ap.add_argument("--receptor", default=None,
                help="cadeias do receptor, separadas por virgula. Quando "
                     "dada, substitui o reconhecimento por palavra-chave no "
                     "COMPND, que falha quando o nome vem por extenso "
                     "(1D3E escreve 'intercellular adhesion molecule-1'). "
                     "A lista sai do organismo declarado no mmCIF.")
ap.add_argument("--proteina", default=None,
                help="alternativa a --cadeia, se o COMPND nomear VP1")
ap.add_argument("--corte", type=float, default=4.5)
A = ap.parse_args()

p = Path(A.entrada)
if not p.exists():
    p = Path(f"{A.entrada.lower()}.pdb")
    if not p.exists():
        u = f"https://files.rcsb.org/download/{A.entrada.upper()}.pdb"
        print("baixando", u)
        urllib.request.urlretrieve(u, p)
pref = p.stem.upper()
print(f"\n=== {pref} ===")

st = gemmi.read_structure(str(p))
st.setup_entities(); st.remove_hydrogens(); st.remove_ligands_and_waters()
cab = compnd(p)
if A.receptor:
    pedidas = [x.strip() for x in A.receptor.split(",") if x.strip()]
    existentes = {c.name for c in st[0]}
    faltam = [x for x in pedidas if x not in existentes]
    if faltam:
        sys.exit(f"\ncadeias de receptor ausentes: {faltam}; "
                 f"no arquivo existem: {sorted(existentes)}")
    rec = [c.name for c in st[0] if c.name in pedidas
           and sum(1 for r in c if r.name in AA) > 0]
else:
    rec = [c.name for c in st[0]
           if any(k in cab.get(c.name, "").upper() for k in REC)
           and sum(1 for r in c if r.name in AA) > 0]
cap = {}
for c in st[0]:
    up = cab.get(c.name, "").upper()
    if c.name in rec or not any(k in up for k in CAP):
        continue
    if sum(1 for r in c if r.name in AA) == 0:
        continue
    cap.setdefault(next((k for k in CAP if k in up), "CAPSIDE"), []).append(c.name)

print("\ncadeias:")
for c in st[0]:
    n = sum(1 for r in c if r.name in AA)
    papel = ("RECEPTOR" if c.name in rec else
             next((k for k, v in cap.items() if c.name in v), "-" if n == 0 else "?"))
    print(f"  {c.name:>3}: {n:4d} res [{papel:>11}] {cab.get(c.name,'(sem COMPND)')}")
if not rec:
    sys.exit("\nnenhum receptor reconhecido no COMPND. Passe --receptor com "
             "as cadeias (veja acha_receptor.py, que as deduz do organismo "
             "declarado no mmCIF) em vez de acrescentar palavra-chave em REC.")

# cadeia alvo: explicita, ou pelo COMPND se ele nomear VP1 de fato
if A.cadeia:
    alvo = A.cadeia
    if alvo not in {c.name for c in st[0]}:
        sys.exit(f"\ncadeia '{alvo}' nao existe")
elif A.proteina and A.proteina in cap:
    alvo = cap[A.proteina][0]
else:
    gen = [k for k in cap if k in ("POLYPROTEIN", "CAPSIDE", "COAT")]
    print("\nnao da para saber qual cadeia e VP1 so pelo COMPND.")
    if gen:
        print(f"  o COMPND diz apenas '{gen[0]}' para as cadeias do capsideo.")
    print("  tamanhos (EV-A71 maduro: VP1~297 VP2~254 VP3~238 VP4~69):")
    for c in st[0]:
        n = sum(1 for r in c if r.name in AA)
        if n and c.name not in rec:
            print(f"    {c.name}: {n} res")
    sys.exit("  rode de novo com --cadeia <letra>. Nao vou adivinhar por "
             "tamanho: foi assim que VP0 virou VP1 antes.")
print(f"\ncadeia alvo: {alvo}  ({cab.get(alvo,'?')})")

# ---- 1. expande o complexo INTEIRO -------------------------------------
copias, origem = expande(st)
n_rec = sum(1 for _, o in copias if o in rec)
n_cap = len(copias) - n_rec
print(f"\n-- montagem --")
if origem is None:
    print("  SEM BIOMT/MTRIX: a casca radial e o bloqueio NAO sao validos")
    print("  (um protomero isolado nao tem geometria de capsideo)")
else:
    print(f"  {origem}: {len(copias)} cadeias "
          f"({n_cap} capsideo + {n_rec} receptor)")

# ---- 2. pegada contra TODAS as copias do receptor ----------------------
ra = [xyz(c) for c, o in copias if o in rec]
ra = np.concatenate(ra) if ra else np.zeros((0, 3), np.float32)
print(f"\n-- pegada (contato < {A.corte} A) --")
print(f"  receptor: {len(ra)} atomos em {n_rec} copia(s)")

# copia identidade da cadeia alvo
pr = next(r for r in st[0][alvo] if r.name in AA)
p0 = np.array([pr[0].pos.x, pr[0].pos.y, pr[0].pos.z])
ident, dm = None, 1e9
for c, o in copias:
    if o != alvo:
        continue
    for r in c:
        if r.name in AA:
            d = np.linalg.norm([r[0].pos.x, r[0].pos.y, r[0].pos.z] - p0)
            if d < dm:
                ident, dm = c, d
            break
if ident is None:
    sys.exit(f"cadeia {alvo} nao sobreviveu a montagem")
print(f"  copia identidade de {alvo}: desvio {dm:.3f} A")

peg = []
if len(ra):
    for r in ident:
        if r.name not in AA:
            continue
        at = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in r], dtype=np.float32)
        if np.sqrt(((at[:, None] - ra[None]) ** 2).sum(-1)).min() < A.corte:
            peg.append(r.seqid.num)
peg.sort()
if not peg:
    sys.exit(f"\no receptor nao toca a cadeia {alvo} nem apos a montagem. "
             f"Tente outra --cadeia.")
print(f"  pegada em {alvo} ({len(peg)} res): {peg}")

# ---- 3. capsideo sem receptor ------------------------------------------
cap_ch = [c for c, o in copias if o not in rec
          and any(r.name in AA for r in c)]
if not cap_ch:
    sys.exit("nenhuma cadeia de capsideo com residuos; confira --cadeia")
co = np.concatenate([xyz(c) for c in cap_ch])
ce = co.mean(axis=0)
car = np.sort(np.array([np.linalg.norm([a.pos.x, a.pos.y, a.pos.z] - ce)
                        for c in cap_ch for r in c if r.name in AA
                        for a in r if a.name == "CA"], dtype=np.float32))
print(f"\n-- capsideo ({len(cap_ch)} cadeias, receptor removido) --")
print(f"  centro sobre {len(co)} atomos")
print(f"  casca {car[0]:.1f}-{car[-1]:.1f} A (espessura {car[-1]-car[0]:.1f} A)")
if car[-1] < 80:
    print("  ATENCAO: raio externo < 80 A. Isto nao e um capsideo de "
          "picornavirus montado; nao leia os numeros abaixo.")

# ---- 4. filtro ---------------------------------------------------------
L = []
for r in ident:
    if r.name not in AA:
        continue
    ca = next((a for a in r if a.name == "CA"), r[0])
    pos = np.array([ca.pos.x, ca.pos.y, ca.pos.z])
    rr = float(np.linalg.norm(pos - ce))
    casca = 100. * float(np.searchsorted(car, rr)) / len(car)
    u = (pos - ce) / rr
    d = co - pos; dist = np.linalg.norm(d, axis=1)
    m = (dist > 1e-6) & (dist < RAIO)
    bloq = int(((d[m] @ u) / dist[m] > CONE).sum())
    bf = float(np.mean([a.b_iso for a in r]))
    L.append((r.seqid.num, r.name, rr, casca, bloq, bf))
print(f"  perfil em {len(L)} residuos")

S = set(peg)
den = [x for x in L if x[0] in S]
fora = [x for x in L if x[0] not in S]
print(f"\n-- EXPOSICAO --")
print(f"  pegada no perfil: {len(den)}/{len(peg)}")
print(f"  casca mediana  pegada {np.median([x[3] for x in den]):6.1f}%  "
      f"resto {np.median([x[3] for x in fora]):6.1f}%")
print(f"  bloqueio median pegada {np.median([x[4] for x in den]):6.0f}   "
      f"resto {np.median([x[4] for x in fora]):6.0f}")
print(f"\n  AUC casca:            {auc([x[3] for x in den], [x[3] for x in fora]):.3f}")
print(f"  AUC bloqueio invert.: {auc([-x[4] for x in den], [-x[4] for x in fora]):.3f}")
ordem = sorted(L, key=lambda x: (-x[3], x[4]))
po = {x[0]: i + 1 for i, x in enumerate(ordem)}
pp = sorted(po[x[0]] for x in den)
n = len(L)
print(f"\n  posto entre {n}: mediana {int(np.median(pp))}, faixa {pp[0]}-{pp[-1]}")
print(f"  no decil superior: {sum(1 for r in pp if r <= n/10)}/{len(pp)} "
      f"(acaso {len(pp)/10:.1f})")
print(f"  passam exposicao (casca>={EXT:.0f}%, bloq<{MAXB}): "
      f"{sum(1 for x in den if x[3] >= EXT and x[4] < MAXB)}/{len(den)}")

bs = np.array([x[5] for x in L])
print(f"\n-- ORDEM (fator B) --")
if bs.std() < 1e-6 or bs.mean() <= 0:
    print("  fator B constante/ausente: eixo de ordem NAO avaliavel aqui")
else:
    bp = float(np.mean([x[5] for x in den]))
    print(f"  B medio {bs.mean():.1f} (dp {bs.std():.1f}, CV {bs.std()/bs.mean():.2f})")
    print(f"  B da pegada {bp:.1f} ({bp/bs.mean():.2f}x a media)")
    print(f"  AUC B invertido: {auc([-x[5] for x in den], [-x[5] for x in fora]):.3f}")

out = Path(f"{pref}_{alvo}_pegada.tsv")
with open(out, "w") as fh:
    fh.write("residue\tresname\tradius_A\tshell_pct\tblocking\tbfactor\t"
             "rank\tin_footprint\tpasses\n")
    for x in sorted(L):
        fh.write(f"{x[0]}\t{x[1]}\t{x[2]:.2f}\t{x[3]:.2f}\t{x[4]}\t{x[5]:.1f}\t"
                 f"{po[x[0]]}\t{'yes' if x[0] in S else 'no'}\t"
                 f"{'yes' if x[3] >= EXT and x[4] < MAXB else 'no'}\n")
print(f"\ngravado {out}")
print("""
So leia os AUC se a casca acima tiver raio externo plausivel (>100 A para
picornavirus). Sem montagem, os numeros nao significam nada.
AUC casca >0.75 -> replica o EV-A71 em sistema independente
AUC casca ~0.5  -> nao generaliza
pegada que NAO passa -> limiar de 70% rejeita sitio real, esta alto demais
O filtro EXCLUI sitio impossivel; nao preve sitio de ligacao.""")
