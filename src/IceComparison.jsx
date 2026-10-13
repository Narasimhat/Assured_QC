import {useState} from 'react';
import {runIce,iceComparisonHtml} from './iceComparison.js';
export default function IceComparison({outcome,files,nuclease,onComplete,onBusy}){
 const [busy,setBusy]=useState(false),[progress,setProgress]=useState('');
 async function run(){
  setBusy(true);onBusy(true);const results=[];
  try{for(let i=0;i<outcome.tasks.length;i++){
   setProgress(`ICE: ${i+1} / ${outcome.tasks.length}`);let ice;
   try{ice=await runIce(outcome.tasks[i],files,nuclease);}catch(e){ice=[{ok:false,label:'Not analysed',error:e.message}];}
   results.push({...outcome.results[i],ice});
  }onComplete(results);setProgress('ICE comparison finished. Review warnings and download the report.');}
  finally{setBusy(false);onBusy(false);}
 }
 return <section className="card"><h3>Compare with ICE (optional)</h3><p>Assured QC runs locally in your browser. Running ICE sends these control/sample AB1 files and guide/donor sequences to this app’s server. Temporary files are removed after the request. Research/nonprofit use only. <a href="/ICE-LICENSE.txt">Licence</a>.</p><button disabled={busy} onClick={run}>{busy?'Running ICE…':'Run ICE comparison on server'}</button><p role="status">{progress}</p><div dangerouslySetInnerHTML={{__html:iceComparisonHtml(outcome.results)}}/>{outcome.results.some(r=>r.ice)&&<button onClick={()=>{const u=URL.createObjectURL(new Blob([JSON.stringify(outcome.results.map(r=>({sample:r.name,ice:r.ice})),null,2)],{type:'application/json'}));const a=document.createElement('a');a.href=u;a.download='ICE_comparison_provenance.json';a.click();setTimeout(()=>URL.revokeObjectURL(u),2000);}}>Download ICE results + provenance</button>}</section>;
}
