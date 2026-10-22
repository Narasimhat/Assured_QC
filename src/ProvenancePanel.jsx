import PassportPanel from './PassportPanel.jsx';
import { useId, useState } from 'react';
import { fields, screen, mint, verify } from './provenance.js';

export default function ProvenancePanel({ source }) {
  const id = useId();
  const [metadata, setMetadata] = useState({});
  const [file, setFile] = useState(null);
  const [proof, setProof] = useState(null);
  const [candidate, setCandidate] = useState(null);
  const [certificate, setCertificate] = useState(null);
  const [message, setMessage] = useState('');
  const [busy, setBusy] = useState(false);
  const assessment = screen(metadata);
  async function run(action) {
    setBusy(true); setMessage('');
    try { await action(); } catch (error) { setMessage(error.message); } finally { setBusy(false); }
  }
  function download() {
    const url = URL.createObjectURL(new Blob([JSON.stringify(certificate, null, 2)], { type: 'application/json' }));
    const link = document.createElement('a'); link.href = url; link.download = 'assured-provenance.json';
    document.body.appendChild(link); link.click(); link.remove(); setTimeout(() => URL.revokeObjectURL(url), 1000);
  }
  return <details style={{ background: '#fff', color: '#153b38', border: '1px solid #b9ceca', borderRadius: 12, padding: 18, margin: '16px 0' }}>
    <summary style={{ cursor: 'pointer', fontWeight: 700 }}>FAIR metadata & provenance notary</summary>
    <p>Register an exported design, QC report, or dataset. All files stay in your browser. Records are self-issued and unanchored.</p>
    <p>Completeness: <strong>{assessment.score}%</strong> · 12 equal-weight checks for field presence and HTTPS URL syntax. This is not formal FAIR certification, scientific QC validation, or regulatory approval. Batch IDs are local identifiers, not registered persistent identifiers.</p>
    <fieldset disabled={busy} style={{ border: 0, padding: 0 }}>
      <legend>Register a file</legend>
      <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(210px, 1fr))', gap: 12 }}>
        {fields.map(([pillar, key, label]) => <label key={key} htmlFor={`${id}-${key}`} style={{ display: 'grid', gap: 4 }}>
          {pillar} · {label}
          <input id={`${id}-${key}`} value={metadata[key] || ''} onChange={e => { setMetadata(m => ({ ...m, [key]: e.target.value })); setCertificate(null); }} style={{ width: '100%', boxSizing: 'border-box' }} />
        </label>)}
      </div>
      <label style={{ display: 'block', margin: '12px 0' }}>File to fingerprint (up to 50 MB) <input type="file" onChange={e => { setFile(e.target.files[0] || null); setCertificate(null); }} /></label>
      <button type="button" disabled={!file || assessment.score !== 100} onClick={() => run(async () => {
        if (file.size > 50 * 1024 * 1024) throw new Error('Choose a file of 50 MB or smaller.');
        setCertificate(await mint(await file.arrayBuffer(), file.name, metadata, source));
        setMessage('Record created. Anchor status: unanchored.');
      })}>Create provenance record</button>
      {certificate && <div><p style={{ overflowWrap: 'anywhere' }}>Record SHA-256: <code>{certificate.recordSha256}</code></p><button type="button" onClick={download}>Download JSON certificate</button></div>}
    </fieldset>
    <fieldset disabled={busy} style={{ marginTop: 18 }}>
      <legend>Verify a received file</legend>
      <label style={{ display: 'block' }}>JSON certificate <input type="file" accept=".json,application/json" onChange={e => { setProof(e.target.files[0] || null); setMessage(''); }} /></label>
      <label style={{ display: 'block' }}>File to compare <input type="file" onChange={e => { setCandidate(e.target.files[0] || null); setMessage(''); }} /></label>
      <button type="button" disabled={!proof || !candidate} onClick={() => run(async () => {
        if (proof.size > 1024 * 1024 || candidate.size > 50 * 1024 * 1024) throw new Error('Limit: 1 MB certificate, 50 MB dataset.');
        const matches = await verify(await candidate.arrayBuffer(), JSON.parse(await proof.text()));
        setMessage(matches ? 'MATCH: file bytes match the supplied record.' : 'MISMATCH: file differs from the supplied record. Check the file and certificate.');
      })}>Verify fingerprints</button>
    </fieldset>
    <p role="status" aria-live="polite">{busy ? 'Computing fingerprints…' : message}</p>
    <p>A hash proves equality to a supplied record, not authorship or data truth. Anyone can replace both a record and its hash. Keep a trusted original independently. Timestamps use your computer clock. No blockchain, digital signature, or external timestamp service is connected. Nothing is retained after leaving the page.</p>
    <PassportPanel metadata={metadata} source={source} />
  </details>;
}
