const NS='http://www.w3.org/2000/svg';
const COLORS={teal:'#087f87',teal2:'#43b6bc',orange:'#e55f2a',orange2:'#f3a27e',ink:'#102533',muted:'#6f828c',grid:'#e8eef0',line:'#9aabb2'};
let DATA;
const $=s=>document.querySelector(s); const $$=s=>[...document.querySelectorAll(s)];

fetch('assets/study-data.json').then(r=>r.json()).then(d=>{DATA=d;init();}).catch(err=>{console.error(err);document.body.classList.add('data-error');});

function init(){
  setupReveal(); drawHero(); drawWave(); renderComponents(); drawDecile(); setupMetricTabs(); renderModelBars('auroc'); renderForest(); setupRisk();
}

function setupReveal(){
  const io=new IntersectionObserver(es=>es.forEach(e=>{if(e.isIntersecting)e.target.classList.add('visible')}),{threshold:.08});
  $$('.reveal').forEach(x=>io.observe(x));
}

function svgEl(tag,attrs={}){const e=document.createElementNS(NS,tag);for(const[k,v]of Object.entries(attrs))e.setAttribute(k,v);return e;}
function sx(v,x0,x1){return x0+v*(x1-x0)}
function sy(v,y0,y1,min,max){return y1-(v-min)/(max-min)*(y1-y0)}
function path(x,y,x0,x1,y0,y1,min,max){return y.map((v,i)=>`${i?'L':'M'}${sx(x[i],x0,x1).toFixed(1)} ${sy(v,y0,y1,min,max).toFixed(1)}`).join(' ')}
function bandPath(x,mean,sd,x0,x1,y0,y1,min,max){
  const hi=mean.map((v,i)=>v+sd[i]),lo=mean.map((v,i)=>v-sd[i]);
  let p=hi.map((v,i)=>`${i?'L':'M'}${sx(x[i],x0,x1).toFixed(1)} ${sy(v,y0,y1,min,max).toFixed(1)}`).join(' ');
  for(let i=lo.length-1;i>=0;i--)p+=` L${sx(x[i],x0,x1).toFixed(1)} ${sy(lo[i],y0,y1,min,max).toFixed(1)}`;
  return p+' Z';
}
function addText(svg,text,x,y,cls,anchor='middle'){const t=svgEl('text',{x,y,class:cls,'text-anchor':anchor});t.textContent=text;svg.append(t);return t}

function drawHero(){
  const svg=$('#heroWave'),ph=DATA.phase,w=DATA.waveforms['No CTRCD'].mean;svg.innerHTML='';
  [0,.25,.5,.75,1].forEach(v=>svg.append(svgEl('line',{x1:sx(v,25,695),x2:sx(v,25,695),y1:22,y2:330,class:'hero-gridline'})));
  [-10,0,10].forEach(v=>svg.append(svgEl('line',{x1:25,x2:695,y1:sy(v,22,330,-14,12),y2:sy(v,22,330,-14,12),class:'hero-gridline'})));
  svg.append(svgEl('path',{d:path(ph,w,25,695,22,330,-14,12),class:'hero-wave'}));
  const i0=ph.findIndex(x=>x>=.45),i1=ph.findIndex(x=>x>=.69);
  svg.append(svgEl('path',{d:path(ph.slice(i0,i1),w.slice(i0,i1),sx(ph[i0],25,695),sx(ph[i1-1],25,695),22,330,-14,12),class:'hero-hot'}));
}

function axes(svg,{w=920,h=430,min=-15,max=12,x0=68,x1=900,y0=24,y1=382}={}){
  [0,.2,.4,.6,.8,1].forEach(v=>{let x=sx(v,x0,x1);svg.append(svgEl('line',{x1:x,x2:x,y1:y0,y2:y1,class:'grid'}));addText(svg,v.toFixed(1),x,h-22,'tick')});
  [-10,-5,0,5,10].forEach(v=>{let y=sy(v,y0,y1,min,max);svg.append(svgEl('line',{x1:x0,x2:x1,y1:y,y2:y,class:v===0?'zero':'grid'}));addText(svg,String(v),x0-10,y+4,'tick','end')});
  addText(svg,'Normalized cardiac cycle',(x0+x1)/2,h-2,'axis-label');
  const yt=addText(svg,'Myocardial velocity (cm/s)',14,(y0+y1)/2,'axis-label');yt.setAttribute('transform',`rotate(-90 14 ${(y0+y1)/2})`);
}

let groupMode='both';
function drawWave(){
  const svg=$('#waveChart'),ph=DATA.phase;svg.innerHTML='';axes(svg,{});
  const region=svgEl('rect',{x:sx(.49,68,900),y:24,width:sx(.69,68,900)-sx(.49,68,900),height:358,class:'e-region'});svg.insertBefore(region,svg.firstChild);
  [['No CTRCD','teal'],['CTRCD','orange']].forEach(([g,c])=>{if(groupMode!=='both'&&groupMode!==g)return;const z=DATA.waveforms[g];svg.append(svgEl('path',{d:bandPath(ph,z.mean,z.sd,68,900,24,382,-15,12),class:`band-${c}`}));svg.append(svgEl('path',{d:path(ph,z.mean,68,900,24,382,-15,12),class:`wave-${c}`}));});
  addText(svg,'S',sx(.18,68,900),sy(9.7,24,382,-15,12),'axis-label');
  addText(svg,"E′ region",sx(.60,68,900),sy(-12.5,24,382,-15,12),'axis-label');
}
$('#groupToggle')?.addEventListener('click',e=>{if(e.target.tagName!=='BUTTON')return;$$('#groupToggle button').forEach(b=>b.classList.remove('active'));e.target.classList.add('active');groupMode=e.target.dataset.group;drawWave();});

function renderComponents(){
  const box=$('#componentButtons');box.innerHTML='';DATA.components.forEach((c,i)=>{const b=document.createElement('button');b.textContent=`${i+1}. ${c.label}`;if(i===0)b.classList.add('active');b.onclick=()=>{[...box.children].forEach(x=>x.classList.remove('active'));b.classList.add('active');showComponent(i)};box.append(b)});showComponent(0);
}
function showComponent(i){
  const c=DATA.components[i],d=$('#componentDetail');
  d.innerHTML=`<h4>${c.label}</h4><p>${c.interpretation}. MTI orients this component so that higher standardized values indicate greater mechanical-transition abnormality.</p><div class="compare"><span>CTRCD median<b>${c.event_median}</b></span><span>No CTRCD median<b>${c.nonevent_median}</b></span></div>`;
}

function drawDecile(){
  const d=+$('#decileSlider').value,z=DATA.waveforms.deciles[d-1],svg=$('#decileChart');svg.innerHTML='';axes(svg,{w:920,h:290,min:-15,max:12,x0:62,x1:900,y0:18,y1:244});
  svg.append(svgEl('path',{d:path(DATA.phase,z.waveform,62,900,18,244,-15,12),class:'wave-teal'}));
  const opacity=.10+.05*d;const r=svgEl('rect',{x:sx(.46,62,900),y:18,width:sx(.68,62,900)-sx(.46,62,900),height:226,fill:`rgba(229,95,42,${Math.min(opacity,.32)})`});svg.insertBefore(r,svg.firstChild);
  $('#decileValue').textContent=d;$('#decileCaption').textContent=`Decile ${d}: n=${z.n}, median descriptive MTI ${z.median_mti.toFixed(2)}, observed CTRCD rate ${(100*z.event_rate).toFixed(1)}%. Deciles are descriptive, not clinical cutoffs.`;
}
$('#decileSlider')?.addEventListener('input',drawDecile);

function setupMetricTabs(){$('#metricTabs')?.addEventListener('click',e=>{if(e.target.tagName!=='BUTTON')return;$$('#metricTabs button').forEach(b=>b.classList.remove('active'));e.target.classList.add('active');renderModelBars(e.target.dataset.metric);});}
function renderModelBars(metric){
  const el=$('#modelBars');el.innerHTML='';const label=metric==='auroc'?'AUROC':'AUPRC';
  DATA.models.forEach(m=>{const v=m[metric];const row=document.createElement('div');row.className='model-row'+(m.model==='MTI + LVEF'?' hot':'');row.innerHTML=`<span>${m.model}</span><div class="bar-track"><div class="bar" style="width:${Math.max(4,v*100)}%"></div></div><b>${v.toFixed(3)}</b>`;el.append(row)});
  el.setAttribute('aria-label',label+' model comparison');
}

function renderForest(){
  const el=$('#forestPlot');el.innerHTML='';const min=.7,max=4.2;
  DATA.cox.forEach(r=>{const row=document.createElement('div');row.className='forest-row';const pos=x=>100*(x-min)/(max-min);row.innerHTML=`<span>${r.label}</span><div class="forest-axis"><i class="forest-ci" style="left:${pos(r.lo)}%;width:${pos(r.hi)-pos(r.lo)}%"></i><i class="forest-dot" style="left:${pos(r.hr)}%"></i></div><b>${r.hr.toFixed(2)}</b>`;el.append(row)});
}

function setupRisk(){
  const slider=$('#riskSlider');slider.addEventListener('input',()=>updateRisk(+slider.value));updateRisk(20);drawRiskChart();
}
function riskRow(pct){return DATA.risk_capture.find(r=>Math.round(r.fraction*100)===pct)||DATA.risk_capture[0]}
function updateRisk(pct){const r=riskRow(pct);$('#capacityPct').textContent=pct+'%';$('#selectedN').textContent=r.selected_n;$('#capturePct').textContent=(100*r.capture).toFixed(1)+'%';$('#selectedPpv').textContent=(100*r.ppv).toFixed(1)+'%';}
function drawRiskChart(){
  const svg=$('#riskChart');svg.innerHTML='';const x0=62,x1=900,y0=18,y1=235;
  [0,.25,.5,.75,1].forEach(v=>{const x=sx(v,x0,x1);svg.append(svgEl('line',{x1:x,x2:x,y1:y0,y2:y1,class:'grid'}));addText(svg,Math.round(v*100)+'%',x,270,'tick')});
  [0,.25,.5,.75,1].forEach(v=>{const y=sy(v,y0,y1,0,1);svg.append(svgEl('line',{x1:x0,x2:x1,y1:y,y2:y,class:'grid'}));addText(svg,Math.round(v*100)+'%',x0-10,y+4,'tick','end')});
  const x=DATA.risk_capture.map(r=>r.fraction),y=DATA.risk_capture.map(r=>r.capture);svg.append(svgEl('path',{d:path(x,y,x0,x1,y0,y1,0,1),class:'wave-teal'}));svg.append(svgEl('path',{d:`M${x0} ${y1} L${x1} ${y0}`,stroke:'#b8c6cc','stroke-dasharray':'5 5',fill:'none'}));
  addText(svg,'Fraction selected',480,279,'axis-label');const yt=addText(svg,'Events captured',13,126,'axis-label');yt.setAttribute('transform','rotate(-90 13 126)');
}
