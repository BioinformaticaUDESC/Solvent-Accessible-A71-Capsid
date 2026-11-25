#!/usr/bin/env python3
"""Piso de um eixo contra piso de dois, mesma divisao de piso_calibra_testa.py."""
import csv, glob, os, random
SEMENTE=20261116; SENS=0.95; FORA=set()
SO=['1V9U', '3DPR', '3J8F', '6ADL', '6ADM', '6ADR', '6CX1', '6EIT', '6I2K', '6ILJ', '6ILK', '6PPO', '6PSF', '6SNW', '7BZT', '7BZU', '7C9W', '7VXZ', '7VY5', '7VY6', '7VYK', '7VYL', '7VYM', '7W14', '7XXA', '8B8R']

comp={}
for p in sorted(glob.glob('*_pegada.tsv')):
    pdb=os.path.basename(p).split('_')[0].upper()
    if pdb not in SO: continue
    rows=[r for r in csv.DictReader(open(p),delimiter='\t') if 'in_footprint' in r]
    if rows: comp[pdb]=rows
ids=sorted(comp); random.Random(SEMENTE).shuffle(ids)
cal,tes=ids[:len(ids)//2], ids[len(ids)//2:]
print(f"calibracao ({len(cal)}): {' '.join(cal)}")
print(f"teste      ({len(tes)}): {' '.join(tes)}\n")

def junta(sub):
    P=[];N=[]
    for i in sub:
        for r in comp[i]:
            try: t=(float(r['shell_pct']), float(r['blocking']))
            except ValueError: continue
            (P if r['in_footprint'].strip().lower()=='yes' else N).append(t)
    return P,N
Pc,Nc=junta(cal); Pt,Nt=junta(tes)

def quantil(v,q):
    v=sorted(v); return v[max(0,min(len(v)-1,int(q*(len(v)-1))))]
piso_s=quantil([s for s,_ in Pc], 1-SENS)      # casca: reter 95% acima
piso_b=quantil([b for _,b in Pc], SENS)        # bloqueio: reter 95% abaixo
print(f"pisos calibrados p/ {SENS:.0%} de sensibilidade:")
print(f"  casca    >= {piso_s:.1f}%")
print(f"  bloqueio <= {piso_b:.0f}\n")

def mede(P,N,regra):
    sens=sum(1 for x in P if regra(x))/len(P)
    espe=1-sum(1 for x in N if regra(x))/len(N)
    return sens,espe
regras={
 'so casca   ': lambda x: x[0]>=piso_s,
 'so bloqueio': lambda x: x[1]<=piso_b,
 'os dois (E)': lambda x: x[0]>=piso_s and x[1]<=piso_b,
 'um ou outro': lambda x: x[0]>=piso_s or  x[1]<=piso_b,
}
for rot,(P,N) in (('calibracao',(Pc,Nc)),('TESTE     ',(Pt,Nt))):
    print(rot)
    for nome,f in regras.items():
        s,e=mede(P,N,f)
        print(f"   {nome}   sensibilidade {s:6.1%}   especificidade {e:6.1%}")
    print()
print("sitio 22 de 6I2K: casca 1,68%  bloqueio 846")
print(f"   so casca    -> {'passa' if 1.68>=piso_s else 'REPROVA'}")
print(f"   so bloqueio -> {'passa' if 846<=piso_b else 'REPROVA'}")
