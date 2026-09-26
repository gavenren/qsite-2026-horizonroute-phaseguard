"""Create figures and a local route viewer from saved, measured results."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

ROOT=Path(__file__).resolve().parent
DATA=json.loads((ROOT/'results/benchmark_results.json').read_text())
OUT=ROOT/'figures';OUT.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False})
names=list(DATA['public'])
methods=['Official baseline','Static placement search','HorizonRoute']
colors=['#a2adb8','#6d87a2','#008879']
fig,ax=plt.subplots(figsize=(11,5.5))
for j,(method,color) in enumerate(zip(methods,colors)):
 values=[next(r['score'] for r in DATA['public'][name]['results'] if r['method']==method) for name in names]
 bars=ax.bar([i+(j-1)*.25 for i in range(len(names))],values,width=.24,label=method,color=color)
 ax.bar_label(bars,fmt='%g',padding=3,fontsize=9)
ax.set_xticks(range(len(names)),[n.replace('_','\n') for n in names]);ax.set_ylabel('Score: SWAPs + 0.5 × depth (lower is better)')
ax.set_title('HorizonRoute | Verified public benchmark results',loc='left',fontsize=18,fontweight='bold',pad=18)
ax.legend(frameon=False);ax.set_ylim(0,145);ax.grid(axis='y',alpha=.12);ax.set_axisbelow(True)
fig.tight_layout();fig.savefig(OUT/'score_comparison.png',dpi=180);plt.close(fig)

case=DATA['public']['ghz_star'];route=next(r for r in case['results'] if r['method']=='HorizonRoute')
steps=[0]+[i+1 for i,op in enumerate(route['route']) if op[0]=='SWAP']
steps=steps[:4]
positions={int(k):v for k,v in DATA['hardware_positions'].items()}
fig,axes=plt.subplots(1,len(steps),figsize=(4*len(steps),5),squeeze=False)
for ax,step in zip(axes[0],steps):
 inverse={int(p):int(q) for q,p in route['placement'].items()}
 for op in route['route'][:step]:
  if op[0]=='SWAP':
   a,b=op[1:];inverse[a],inverse[b]=inverse.get(b),inverse.get(a)
 for a,b in DATA['hardware_edges']:
  ax.plot([positions[a][0],positions[b][0]],[positions[a][1],positions[b][1]],color='#c7d0d8',lw=2,zorder=0)
 for p,(x,y) in positions.items():
  logical=inverse.get(p)
  ax.scatter(x,y,s=700,color='#008879' if logical==0 else '#d2e9e6' if logical is not None else '#eef1f4',edgecolors='white',zorder=2)
  ax.text(x,y, f'L{logical}' if logical is not None else '·',ha='center',va='center',color='white' if logical==0 else '#153c43',fontsize=11,fontweight='bold')
  ax.text(x+.22,y-.24,str(p),ha='center',va='center',color='#83929e',fontsize=7)
 if step:
  _,a,b=route['route'][step-1];x,y=positions[a];u,v=positions[b]
  ax.add_patch(FancyArrowPatch((x,y),(u,v),arrowstyle='<->',mutation_scale=20,color='#d88135',lw=3,zorder=3))
 ax.set_title('Initial placement' if step==0 else f'After operation {step}',fontsize=12)
 ax.set_xlim(-.5,3.5);ax.set_ylim(-.5,4.5);ax.set_aspect('equal');ax.axis('off')
fig.suptitle(f"GHZ route: {route['swap_count']} SWAPs, depth {route['depth']}, score {route['score']:g}",fontsize=17,fontweight='bold')
fig.tight_layout();fig.savefig(OUT/'ghz_route.png',dpi=180);plt.close(fig)

TEMPLATE=r'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>HorizonRoute | Route viewer</title>
<style>body{font:16px system-ui;background:#edf3f4;color:#16383f;margin:0}main{max-width:1050px;margin:40px auto;padding:30px;background:white;border-radius:24px}h1{font-size:42px;margin:0}small{color:#597079}select,button{font:inherit;padding:10px;border:1px solid #bcccd0;border-radius:10px;background:white;margin:8px 8px 8px 0}button{cursor:pointer}section{display:grid;grid-template-columns:1fr 1fr;gap:28px}.card{padding:24px;background:#f1f7f6;border-radius:16px}svg{width:100%;height:470px}input{width:100%}pre{height:350px;overflow:auto;font-size:14px;line-height:1.65;white-space:pre-wrap}.active{background:#cbebe4;font-weight:bold}strong{font-size:23px}@media(max-width:700px){section{grid-template-columns:1fr}main{margin:0;padding:20px}h1{font-size:32px}}</style>
<main><small>QSITE 2026 · COMPUTATIONAL TRACK</small><h1>HorizonRoute</h1><p>Follow each logical qubit through a verified route. Green is a program gate. Orange is a SWAP.</p>
<select id="c"></select><select id="m"><option>HorizonRoute</option><option>Official baseline</option><option>100 random placements</option></select><p id="stats"></p>
<section><div class="card"><svg id="graph" viewBox="0 0 400 470"></svg><button id="play">Play</button><button id="back">Back</button><button id="next">Next</button><input type="range" id="step" min="0" value="0"><p id="status"></p></div><div class="card"><strong>Ordered operations</strong><p>Physical qubit indices appear below. The route preserves every input operation and its operand order.</p><pre id="ops"></pre></div></section><p><small>All shown routes pass the original challenge scorer. The local viewer uses the saved test output. It needs no server.</small></p></main>
<script>const data=__DATA__;const $=id=>document.getElementById(id);let timer=null;for(const k of Object.keys(data.public))$('c').add(new Option(k,k));function item(){return data.public[$('c').value].results.find(x=>x.method===$('m').value)}function load(){$('step').max=item().route.length;$('step').value=0;draw()}function draw(){let r=item(),step=+$('step').value,map={};for(let [l,p]of Object.entries(r.placement))map[p]=l;for(let op of r.route.slice(0,step))if(op[0]==='SWAP')[map[op[1]],map[op[2]]]=[map[op[2]],map[op[1]]];let op=r.route[step],out='';let xy=p=>{let z=data.hardware_positions[p];return[65+90*z[0],410-85*z[1]]};for(let[a,b]of data.hardware_edges){let[x,y]=xy(a),[u,v]=xy(b),active=op&&op.length===3&&op.slice(1).includes(a)&&op.slice(1).includes(b);out+=`<line x1="${x}" y1="${y}" x2="${u}" y2="${v}" stroke="${active?(op[0]==='SWAP'?'#d88135':'#008879'):'#bacbd0'}" stroke-width="${active?7:3}"/>`}for(let p of Object.keys(data.hardware_positions)){let[x,y]=xy(p),l=map[p];out+=`<circle cx="${x}" cy="${y}" r="23" fill="${l===undefined?'#e4ecef':l==='0'?'#008879':'#d1eae6'}"/><text x="${x}" y="${y+5}" text-anchor="middle" font-size="14" fill="${l==='0'?'white':'#153c43'}">${l===undefined?'·':'L'+l}</text><text x="${x+24}" y="${y+25}" font-size="10" fill="#637d85">${p}</text>`}$('graph').innerHTML=out;$('stats').innerHTML=`<strong>${r.score} score</strong> · ${r.swap_count} SWAPs · ${r.depth} layers · Official check: PASS`;$('status').textContent=step===r.route.length?'Route complete.':`Before operation ${step+1} of ${r.route.length}: ${op.join(' ')}`;$('ops').innerHTML=r.route.map((v,i)=>`<span class="${i===step?'active':''}">${String(i+1).padStart(2)}  ${v.join(' ')}</span>`).join('\n');let active=$('ops').querySelector('.active');if(active)$('ops').scrollTop=Math.max(0,active.offsetTop-$('ops').offsetTop-100)}$('c').onchange=$('m').onchange=load;$('step').oninput=draw;$('next').onclick=()=>{$('step').value=Math.min(+$('step').value+1,+$('step').max);draw()};$('back').onclick=()=>{$('step').value=Math.max(0,+$('step').value-1);draw()};$('play').onclick=()=>{if(timer){clearInterval(timer);timer=null;$('play').textContent='Play'}else{$('play').textContent='Pause';timer=setInterval(()=>{if(+$('step').value===+$('step').max){clearInterval(timer);timer=null;$('play').textContent='Play'}else $('next').click()},800)}};load();</script></html>'''
(ROOT/'demo.html').write_text(TEMPLATE.replace('__DATA__',json.dumps(DATA)))
print('Saved score_comparison.png, ghz_route.png, and demo.html')
