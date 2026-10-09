import { useEffect, useMemo, useRef, useState } from "react";
import { WORKFLOWS, headline } from "./batch.js";
import Detail from "./Detail.jsx";
import QuickStart from "./QuickStart.jsx";
import { buildQcReportHtml, buildResultsCsv } from "./report.js";

const FLOW_KEYS = Object.keys(WORKFLOWS);
const FLOW_HINTS = {
  knockout: "Pool or clone edited with one guide. Gives editing %, knockout score, indel profile, clone genotype.",
  deletion: "Two or three guides cut together. Also finds the deletion between the cuts.",
  snp: "ssODN with silent blocking changes. Separates the intended edit from blocking-changes-only conversion.",
  tag_small: "Insert up to about 100 bp in an ssODN (donor up to 300 nt with 15 nt arms).",
  reporter: "Tags and reporters too large for deconvolution. Reads junction PCR products against the expected knock-in sequence; checks out-out band sizes.",
};

function download(name, text, type) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const link = document.createElement("a");
  link.href = url; link.download = name; document.body.appendChild(link); link.click(); link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}

const readFile = async (file) => ({ name: file.name, buffer: await file.arrayBuffer() });

export default function App() {
  const [workflow, setWorkflow] = useState("snp");
  const [designSpec, setDesignSpec] = useState(null);
  const [designName, setDesignName] = useState("");
  const [designError, setDesignError] = useState("");
  const [guides, setGuides] = useState("");
  const [donor, setDonor] = useState("");
  const [gene, setGene] = useState("");
  const [controlFile, setControlFile] = useState(null);
  const [editedFiles, setEditedFiles] = useState([]);
  const [sampleType, setSampleType] = useState("clone");
  const [bandsText, setBandsText] = useState("");
  const [running, setRunning] = useState(false);
  const [progress, setProgress] = useState(null);
  const [outcome, setOutcome] = useState(null);
  const [picked, setPicked] = useState(0);
  const worker = useRef(null);
  const flow = WORKFLOWS[workflow];
  const isJunction = flow.method === "junction";

  useEffect(() => () => { if (worker.current) worker.current.terminate(); }, []);

  const onDesign = async (event) => {
    const file = event.target.files[0];
    setDesignError(""); setDesignSpec(null); setDesignName("");
    if (!file) return;
    try { setDesignSpec(JSON.parse(await file.text())); setDesignName(file.name); } catch (error) { setDesignError(`That file is not valid JSON: ${error.message}`); }
  };

  const canRun = useMemo(() => {
    if (running || !editedFiles.length) return false;
    if (isJunction) return Boolean(designSpec);
    return Boolean(controlFile) && (Boolean(designSpec) || guides.trim().length > 0);
  }, [running, editedFiles, isJunction, designSpec, controlFile, guides]);

  const run = async () => {
    setRunning(true); setOutcome(null); setPicked(0); setProgress({ done: 0, total: editedFiles.length, name: "" });
    try {
      const input = {
        workflow, designSpec, manual: { guides, donor, gene: gene || "sample" }, sampleType, bandsText,
        control: controlFile && !isJunction ? await readFile(controlFile) : null,
        samples: await Promise.all(editedFiles.map(readFile)),
      };
      if (worker.current) worker.current.terminate();
      worker.current = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });
      worker.current.onmessage = (event) => {
        if (event.data.progress) setProgress(event.data.progress);
        if (event.data.result) { setOutcome(event.data.result); setRunning(false); setProgress(null); }
      };
      worker.current.onerror = (event) => { setOutcome({ ok: false, error: `The analysis stopped: ${event.message || "worker error"}` }); setRunning(false); setProgress(null); };
      worker.current.postMessage({ id: 1, input });
    } catch (error) {
      setOutcome({ ok: false, error: error.message }); setRunning(false); setProgress(null);
    }
  };

  const designLabel = outcome && outcome.ok ? `${outcome.design.gene || ""} ${outcome.design.label || ""}`.trim() : "";
  const stamp = new Date().toISOString().slice(0, 10);

  return (
    <div className="page">
      <h1>Assured QC</h1>
      <p className="muted">Sanger-trace QC for CRISPR knockouts, deletions, SNP corrections and knock-ins. Files are read and analysed in your browser; nothing is uploaded.</p>

      <QuickStart />

      <details className="card"><summary><b>Or set everything by hand (one control, one design, one workflow)</b></summary>
      <div>
        <h2>1. What was edited</h2>
        <div className="workflows">
          {FLOW_KEYS.map((key) => (
            <button key={key} type="button" className={`workflow${key === workflow ? " active" : ""}`} onClick={() => { setWorkflow(key); setOutcome(null); }}>
              <b>{WORKFLOWS[key].label}</b><br /><span className="muted">{FLOW_HINTS[key]}</span>
            </button>
          ))}
        </div>
      </div>

      <div className="card">
        <h2>2. Design</h2>
        <label className="field" htmlFor="design">Design file exported from the design app (recommended)</label>
        <input id="design" type="file" accept=".json,application/json" onChange={onDesign} />
        {designName && <p className="muted">Loaded {designName}.</p>}
        {designError && <p className="error">{designError}</p>}
        {!isJunction && !designSpec && (
          <>
            <p className="muted">Or enter the guide sequence(s) yourself (the sequence the control trace is searched for):</p>
            <label className="field" htmlFor="guides">Guide sequence(s), comma separated (up to three)</label>
            <input id="guides" type="text" value={guides} onChange={(e) => setGuides(e.target.value)} placeholder="e.g. AACCAGTTGCAGGCGCCCCA" />
            {flow.needsDonor && (
              <>
                <label className="field" htmlFor="donor">Donor sequence (up to 300 nt, at least 15 nt of homology on both ends)</label>
                <textarea id="donor" value={donor} onChange={(e) => setDonor(e.target.value)} />
              </>
            )}
            <label className="field" htmlFor="gene">Label (optional)</label>
            <input id="gene" type="text" value={gene} onChange={(e) => setGene(e.target.value)} />
          </>
        )}
        {isJunction && !designSpec && <p className="note">This check compares reads with the expected knock-in sequence, which comes from the exported design file.</p>}
      </div>

      <div className="card">
        <h2>3. Traces (.ab1)</h2>
        {!isJunction && (
          <>
            <label className="field" htmlFor="control">Control (unedited cells, same primers)</label>
            <input id="control" type="file" accept=".ab1" onChange={(e) => setControlFile(e.target.files[0] || null)} />
          </>
        )}
        <label className="field" htmlFor="edited">{isJunction ? "Junction or out-out reads" : "Edited samples (pools or clones)"}</label>
        <input id="edited" type="file" accept=".ab1" multiple onChange={(e) => setEditedFiles(Array.from(e.target.files))} />
        <p className="muted">{editedFiles.length} file{editedFiles.length === 1 ? "" : "s"} selected.</p>
        {!isJunction && (
          <>
            <label className="field" htmlFor="stype">Sample type</label>
            <select id="stype" value={sampleType} onChange={(e) => setSampleType(e.target.value)}>
              <option value="clone">Clones (call homozygous, heterozygous, compound)</option>
              <option value="pool">Pools (report fractions only)</option>
            </select>
          </>
        )}
        {isJunction && (
          <>
            <label className="field" htmlFor="bands">Out-out PCR band sizes in bp (optional; one line per sample, e.g. “clone1: 841, 631”)</label>
            <textarea id="bands" value={bandsText} onChange={(e) => setBandsText(e.target.value)} />
          </>
        )}
        <p><button type="button" className="primary" disabled={!canRun} onClick={run}>{running ? "Analysing..." : "Analyse"}</button></p>
        {progress && <p className="muted">{progress.done} of {progress.total} done {progress.name ? `(${progress.name})` : ""}</p>}
      </div>

      {outcome && !outcome.ok && <div className="card"><p className="error">{outcome.error}</p></div>}
      {outcome && outcome.ok && (
        <div className="card">
          <h2>Results</h2>
          {outcome.notes.map((n) => <p key={n} className="note">{n}</p>)}
          <p>
            <button type="button" className="link" onClick={() => download(`assured_qc_${stamp}.html`, buildQcReportHtml({ workflow, design: designLabel, results: outcome.results }), "text/html")}>Download report (HTML)</button>
            <button type="button" className="link" onClick={() => download(`assured_qc_${stamp}.csv`, buildResultsCsv(workflow, outcome.results), "text/csv")}>Download table (CSV)</button>
          </p>
          <table>
            <thead><tr><th>Sample</th><th>Result</th><th>Notes</th></tr></thead>
            <tbody>
              {outcome.results.map((r, i) => (
                <tr key={r.name + i} className={`pick${i === picked ? " active" : ""}`} onClick={() => setPicked(i)}>
                  <td>{r.name}</td>
                  <td className={!r.ok ? "fail" : r.kind === "junction" ? r.verdict.status : (r.warnings.length ? "warn" : "pass")}>{headline(workflow, r)}{r.bands ? ` | bands: ${r.bands.call}` : ""}</td>
                  <td>{r.ok ? (r.warnings && r.warnings.length ? `${r.warnings.length} note${r.warnings.length === 1 ? "" : "s"}` : "") : "not analysed"}</td>
                </tr>
              ))}
            </tbody>
          </table>
          <p className="muted">Click a row for the allele breakdown. These are computed estimates from Sanger traces: check the traces before accepting a clone.</p>
        </div>
      )}
      {outcome && outcome.ok && <Detail workflow={workflow} result={outcome.results[picked]} />}
      </details>
    </div>
  );
}
