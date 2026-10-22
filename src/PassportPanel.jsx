import { useId, useState } from 'react';
import { checkPassport, evidenceRoles, mintPassport, strategies, verifyEvidence, zygosities } from './passport.js';

const initial = { cloneId: '', parentalLine: '', targetLocus: '', assembly: '', strategy: strategies[0],
  spacer: '', target: '', efficiency: '', zygosity: zygosities[0], tool: '', toolVersion: '', parameters: '{}' };
export default function PassportPanel({ metadata, source }) {
  const id = useId();
  const [mode, setMode] = useState('create');
  const [input, setInput] = useState(initial);
  const [files, setFiles] = useState({});
  const [ledger, setLedger] = useState([]);
  const [selected, setSelected] = useState('');
  const [imported, setImported] = useState(null);
  const [entry, setEntry] = useState(0);
  const [candidate, setCandidate] = useState(null);
  const [busy, setBusy] = useState(false);
  const [message, setMessage] = useState('');
  const certificate = selected === 'import' ? imported : ledger.find(c => c.recordSha256 === selected);
  async function run(action) {
    setBusy(true); setMessage('');
    try { await action(); } catch (error) { setMessage(error.message); } finally { setBusy(false); }
  }
  const update = (key, value) => setInput(current => ({ ...current, [key]: value }));
  function download(cert) {
    const url = URL.createObjectURL(new Blob([JSON.stringify(cert, null, 2)], { type: 'application/json' }));
    const link = document.createElement('a'); link.href = url; link.download = 'assured-ipsc-passport.json';
    document.body.appendChild(link); link.click(); link.remove(); setTimeout(() => URL.revokeObjectURL(url), 1000);
  }
  function text(key, label, type = 'text') {
    return <label key={key} style={{ display: 'grid', gap: 4 }} htmlFor={`${id}-${key}`}>{label}<input id={`${id}-${key}`} type={type} value={input[key]} onChange={e => update(key, e.target.value)} /></label>;
  }
  return <section style={{ borderTop: '1px solid #b9ceca', marginTop: 20, paddingTop: 16 }}>
    <h3>iPSC Edit Passport</h3>
    <p>Link parental line → design → analysis → QC evidence. Complete the FAIR metadata fields above for the passport. Records stay in this page’s memory until downloaded; reloading clears the session history.</p>
    <fieldset disabled={busy} style={{ border: 0, padding: 0 }}><legend>Passport workspace</legend>
      <label><input type="radio" name={`${id}-mode`} checked={mode === 'create'} onChange={() => { setMode('create'); setMessage(''); }} /> Creator</label>{' '}
      <label><input type="radio" name={`${id}-mode`} checked={mode === 'audit'} onChange={() => { setMode('audit'); setMessage(''); }} /> Auditor</label>
      {mode === 'create' ? <>
        <p>SpCas9 only. A 20-base spacer alone cannot establish a PAM site or genome-wide specificity. Target locus accepts a gene symbol or named locus (for example AAVS1).</p>
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(210px, 1fr))', gap: 12 }}>
          {text('cloneId', 'Clone identifier')}{text('parentalLine', 'Parental line / RRID / hPSCreg ID')}
          {text('targetLocus', 'Target gene or locus')}{text('assembly', 'Reference assembly and accession')}
          <label>Edit strategy<select value={input.strategy} onChange={e => update('strategy', e.target.value)}>{strategies.map(v => <option key={v}>{v}</option>)}</select></label>
          {text('spacer', 'SpCas9 spacer, 5′→3′ (20 bases)')}{text('target', 'Optional matching spacer + NGG target (23 bases)')}
          {text('efficiency', 'Reported editing fraction %, blank if unknown', 'number')}
          <label>Reported zygosity<select value={input.zygosity} onChange={e => update('zygosity', e.target.value)}>{zygosities.map(v => <option key={v}>{v}</option>)}</select></label>
          {text('tool', 'Analysis tool')}{text('toolVersion', 'Analysis tool / engine version')}
        </div>
        <label style={{ display: 'grid', margin: '12px 0' }}>Analysis parameters (JSON object)<textarea value={input.parameters} onChange={e => update('parameters', e.target.value)} rows={3} /></label>
        <p>Attach exact exported files; their contents are fingerprinted, not interpreted. Design, raw sequencing, and analysis files are required. Optional QC attachments are recorded as “not assessed”, never as passed. Up to 20 files, 50 MB each, 200 MB total.</p>
        {evidenceRoles.map(role => <label key={role} style={{ display: 'block', margin: '8px 0' }}>{role}{['design', 'raw-sequencing', 'analysis'].includes(role) ? ' (required)' : ' (optional)'} <input type="file" multiple onChange={e => setFiles(current => ({ ...current, [role]: Array.from(e.target.files) }))} /></label>)}
        <button type="button" onClick={() => run(async () => {
          const attachments = Object.entries(files).flatMap(([role, values]) => values.map(file => ({ role, file })));
          const cert = await mintPassport(input, metadata, attachments, source);
          setLedger(current => [...current, cert]); setSelected(cert.recordSha256); setEntry(0);
          setMessage(`Passport created for ${cert.record.cloneId}. Self-issued and unanchored. Download it before leaving.`);
        })}>Create iPSC passport</button>
      </> : <>
        <label style={{ display: 'block' }}>Import passport JSON <input type="file" accept=".json,application/json" onChange={e => {
          const file = e.target.files[0]; setImported(null); setSelected(''); setEntry(0); setMessage('');
          if (file) run(async () => {
            if (file.size > 1024 * 1024) throw new Error('Passport limit is 1 MB.');
            const cert = JSON.parse(await file.text()); await checkPassport(cert);
            setImported(cert); setSelected('import'); setMessage('Passport integrity checked. Issuer and biological claims remain unverified.');
          });
        }} /></label>
      </>}
      <label style={{ display: 'block', margin: '12px 0' }}>Session passport / imported record <select value={selected} onChange={e => { setSelected(e.target.value); setEntry(0); setMessage(''); }}>
        <option value="">Select a passport</option>
        {ledger.map((cert, i) => <option key={cert.recordSha256} value={cert.recordSha256}>{cert.record.cloneId} · record {i + 1} · {cert.record.timestampUtc}</option>)}
        {imported && <option value="import">Imported: {imported.record.cloneId}</option>}
      </select></label>
      {certificate && <>
        <p style={{ overflowWrap: 'anywhere' }}>Record SHA-256: <code>{certificate.recordSha256}</code></p>
        <button type="button" onClick={() => download(certificate)}>Download iPSC passport JSON</button>
        <details><summary>Inspect recorded claims and QC evidence</summary><pre style={{ whiteSpace: 'pre-wrap', overflowWrap: 'anywhere' }}>{JSON.stringify(certificate.record, null, 2)}</pre></details>
        {mode === 'audit' && <>
          <label style={{ display: 'block' }}>Evidence entry <select value={entry} onChange={e => { setEntry(Number(e.target.value)); setMessage(''); }}>{certificate.record.evidence.map((ev, i) => <option value={i} key={i}>{ev.role} · {ev.fileName}</option>)}</select></label>
          <label style={{ display: 'block' }}>File to compare <input type="file" onChange={e => { setCandidate(e.target.files[0] || null); setMessage(''); }} /></label>
          <button type="button" disabled={!candidate} onClick={() => run(async () => {
            if (candidate.size > 50 * 1024 * 1024) throw new Error('File limit is 50 MB.');
            const result = await verifyEvidence(await candidate.arrayBuffer(), certificate, entry);
            setMessage(`${result.matches ? 'MATCH: file bytes match this evidence entry.' : 'MISMATCH: file differs from this evidence entry.'} Expected SHA-256: ${result.expectedHash}. Received: ${result.actualHash}. This checks one file, not the whole evidence bundle.`);
          })}>Verify evidence file</button>
        </>}
      </>}
    </fieldset>
    <p role="status" aria-live="polite" style={{ overflowWrap: 'anywhere' }}>{busy ? 'Processing passport…' : message}</p>
    <p>No biological identity, measured editing efficiency, zygosity, pluripotency, karyotype, or off-target claim is established by a matching hash. This is versioned JSON, not an established JSON-LD standard, immutable ledger, or independently signed certificate.</p>
  </section>;
}
