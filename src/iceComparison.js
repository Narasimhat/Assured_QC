const esc = s => String(s ?? '').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const number = (v, digits=1) => Number.isFinite(v) ? v.toFixed(digits) : 'Unavailable';
export function iceCases(task) {
  if (task.controls.length !== 1 || task.edited.length !== 1) throw Error('ICE comparison requires one control and one sample read per clone.');
  if (!['knockout','deletion','snp','tag_small'].includes(task.workflow)) throw Error('This workflow is not supported by the ICE comparison.');
  const guides = task.spec.guides.map(g=>g.spacer);
  if (!guides.length || guides.length > 3 || guides.some(g=>! /^[ACGT]{17,23}$/.test(g))) throw Error('ICE needs explicit guide sequences; inferred cuts are not supported.');
  const donors = task.spec.donors || [];
  if (donors.length > 2) throw Error('ICE comparison supports at most two separate donor hypotheses.');
  return (donors.length ? donors : [{name:'No donor',sequence:null}]).map(d=>({label:d.name || 'Donor',settings:{guides,donor:d.sequence || null,nuclease:'SpCas9'}}));
}
export async function runIce(task, files, nuclease='SpCas9') {
  if((task.iceNuclease || nuclease) !== 'SpCas9') throw Error('Public ICE supports SpCas9 only.');
  const cases=iceCases(task), comparisons=[];
  for(const c of cases){
    const form=new FormData();form.append('control',files.get(task.controls[0]));form.append('sample',files.get(task.edited[0]));form.append('settings',JSON.stringify(c.settings));
    const abort=new AbortController(),timer=setTimeout(()=>abort.abort(),65000);
    try{
      const response=await fetch('/api/ice',{method:'POST',body:form,signal:abort.signal});
      const data=await response.json();
      if(!response.ok || !data.ok) throw Error(data.error || 'ICE failed.');
      const hash=async file=>Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256',await file.arrayBuffer()))).map(v=>v.toString(16).padStart(2,'0')).join('');
      if(data.control_sha256!==await hash(files.get(task.controls[0])) || data.sample_sha256!==await hash(files.get(task.edited[0])))throw Error('ICE input identity mismatch.');
      comparisons.push({...data,label:c.label});
    }catch(e){comparisons.push({ok:false,label:c.label,error:e.name==='AbortError'?'ICE timed out.':e.message});}
    finally{clearTimeout(timer);}
  }
  return comparisons;
}
export function iceComparisonHtml(results){
  if(!results.some(r=>r.ice))return '';
  const rows=results.flatMap(r=>(r.ice||[]).map(c=>{
    const i=c.result||{};
    return `<tr><td>${esc(r.name)}</td><td>${esc(c.label)}</td><td>${number(r.summary?.indelPct)}</td><td>${number(r.summary?.koScorePct)}</td><td>${number(i.ice)}</td><td>${number(i.ko_score)}</td><td>${number(i.hdr_pct)}</td><td>${number(i.rsq,2)}</td><td>${esc(c.error||i.notes||'')}</td></tr>`;
  })).join('');
  return `<section><h2>ICE comparison · research use</h2><p>Public Synthego ICE 1.2.0, source commit 310495fe536b7d2eafcee80f145c55f7c65d947c. Separate estimates; no averaging or automatic change to Assured QC recommendations. This is not the current hosted ICE service.</p><div style="overflow:auto"><table><thead><tr><th>Sample</th><th>Donor hypothesis</th><th>QC indel %</th><th>QC KO screening %</th><th>ICE editing %</th><th>ICE KO %</th><th>ICE HDR %</th><th>ICE R²</th><th>Warnings</th></tr></thead><tbody>${rows}</tbody></table></div><p>Each donor is fitted separately, not pooled. ICE HDR and Assured QC intended-edit estimates have different definitions. Missing or failed results are unavailable, never zero. KO scores predict disruption, not confirmed protein loss. <a href="https://github.com/synthego-open/ice/blob/master/LICENSE">ICE research licence</a>.</p></section>`;
}
