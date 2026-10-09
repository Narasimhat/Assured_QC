import { useMemo, useRef, useState } from "react";
import { classifyFiles, buildPairs, readDesign, planTasks, runPool } from "./quick.js";
import { parseAbif } from "./abif.js";
import { headline } from "./batch.js";
import { buildQcReportHtml, buildResultsCsv } from "./report.js";
import { buildPlateReportHtml, resequenceListCsv } from "./plateReport.js";
import { tally, DECISION_ORDER } from "./decision.js";
import { NUCLEASES } from "./nucleases.js";
import Detail from "./Detail.jsx";

const LABEL = { accept: "Accept", hold: "Hold", "re-sequence": "Sequence again", reject: "Reject" };
const download = (name, data, type) => { const url = URL.createObjectURL(new Blob([data], { type })); const a = document.createElement("a"); a.href = url; a.download = name; document.body.appendChild(a); a.click(); a.remove(); setTimeout(() => URL.revokeObjectURL(url), 2000); };
const arrayBufferOf = (file) => file.arrayBuffer();
const interval = (iv) => (iv ? `${iv[0]}-${iv[1]}%` : "");

export default function QuickStart() {
  const files = useRef(new Map()); const buffers = useRef(new Map());
  const [traceNames, setTraceNames] = useState([]); const [sheetText, setSheetText] = useState(""); const [designText, setDesignText] = useState(""); const [designName, setDesignName] = useState("");
  const [guidesText, setGuidesText] = useState(""); const [donorText, setDonorText] = useState(""); const [nuclease, setNuclease] = useState("SpCas9"); const [typeOverride, setTypeOverride] = useState(""); const [goal, setGoal] = useState("homozygous");
  const [pairsEdit, setPairsEdit] = useState({}); const [running, setRunning] = useState(false); const [progress, setProgress] = useState(null);
  const [outcome, setOutcome] = useState(null); const [picked, setPicked] = useState(0); const [figureBusy, setFigureBusy] = useState(false); const [error, setError] = useState("");
  const [dragging, setDragging] = useState(false);

  const addFiles = async (list) => {
    setError(""); const incoming = Array.from(list); const kinds = classifyFiles(incoming.map((f) => f.name));
    incoming.forEach((f) => files.current.set(f.name, f));
    setTraceNames(Array.from(new Set([...traceNames, ...kinds.traces]))); setPairsEdit({}); setOutcome(null);
    const sheet = incoming.find((f) => kinds.sheets.includes(f.name)); if (sheet) setSheetText(await sheet.text());
    const design = incoming.find((f) => kinds.designs.includes(f.name)); if (design) { setDesignText(await design.text()); setDesignName(design.name); }
  };

  const built = useMemo(() => (traceNames.length ? buildPairs({ traceNames, sheetText }) : { pairs: [], problems: [] }), [traceNames, sheetText]);
  const design = useMemo(() => { try { return { ...readDesign(designText), error: "" }; } catch (e) { return { kind: "none", error: e.message }; } }, [designText]);
  const pairs = built.pairs.map((p, i) => ({ ...p, ...(pairsEdit[i] || {}) }));
  const guides = guidesText.split(/[,;\s]+/).filter(Boolean);
  const lowCount = pairs.filter((p) => p.confidence === "low").length;
  const guideSource = design.kind === "spec" ? `the design file (${design.label})` : design.kind === "report" ? `the design report (${design.label}: ${design.guides.length} guide${design.guides.length === 1 ? "" : "s"})` : guides.length ? "the guides typed below" : "none: the cut site will be inferred from where the traces depart from the control";

  const run = async () => {
    setRunning(true); setOutcome(null); setPicked(0); setError("");
    try {
      const needed = new Set(pairs.flatMap((p) => [...p.control, ...p.edited]));
      for (const name of needed) if (!buffers.current.has(name)) buffers.current.set(name, await arrayBufferOf(files.current.get(name)));
      const parsed = new Map(); const getTrace = (name) => { if (!parsed.has(name)) parsed.set(name, parseAbif(buffers.current.get(name).slice(0), name)); return parsed.get(name); };
      const { tasks, problems } = planTasks({ pairs, getTrace, design, guides, donor: donorText.replace(/\s+/g, ""), nuclease, sampleType: typeOverride, goal });
      if (!tasks.length) { setError(problems.join(" ") || "Nothing to analyse."); setRunning(false); return; }
      setProgress({ done: 0, total: tasks.length });
      const size = Math.max(1, Math.min(8, (navigator.hardwareConcurrency || 4) - 1));
      const results = await runPool({ tasks, size, spawn: () => new Worker(new URL("./worker.js", import.meta.url), { type: "module" }), onProgress: setProgress,
        payload: (task, index) => ({ task: { ...task, index, controls: task.controls.map((n) => ({ name: n, buffer: buffers.current.get(n).slice(0) })), edited: task.edited.map((n) => ({ name: n, buffer: buffers.current.get(n).slice(0) })) } }) });
      setOutcome({ results, tasks, problems: [...built.problems, ...problems], workflow: tasks[0].workflow });
    } catch (e) { setError(e.message || String(e)); }
    setRunning(false); setProgress(null);
  };

  const result = outcome ? outcome.results[picked] : null; const counts = outcome ? tally(outcome.results.map((r) => r.decision)) : null;
  const stamp = new Date().toISOString().slice(0, 10);
  const figure = async (format) => {
    if (!outcome || !result || !result.ok) return; setFigureBusy(true);
    try {
      const task = outcome.tasks[picked]; const worker = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });
      const svg = await new Promise((resolve, reject) => { worker.onmessage = (e) => (e.data.figure ? resolve(e.data.figure) : reject(new Error(e.data.figureError))); worker.postMessage({ id: 1, figure: { name: task.name, spec: task.spec, sampleType: task.sampleType, workflow: task.workflow, control: { name: task.controls[0], buffer: buffers.current.get(task.controls[0]).slice(0) }, edited: { name: task.edited[0], buffer: buffers.current.get(task.edited[0]).slice(0) } } }); });
      worker.terminate();
      if (format === "svg") download(`${task.name}_figure.svg`, svg, "image/svg+xml");
      else { const scale = 300 / 96; const url = URL.createObjectURL(new Blob([svg], { type: "image/svg+xml" })); const img = new Image(); await new Promise((ok, bad) => { img.onload = ok; img.onerror = bad; img.src = url; }); const canvas = document.createElement("canvas"); canvas.width = Math.round(img.width * scale); canvas.height = Math.round(img.height * scale); const ctx = canvas.getContext("2d"); ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, canvas.width, canvas.height); ctx.drawImage(img, 0, 0, canvas.width, canvas.height); URL.revokeObjectURL(url); canvas.toBlob((blob) => download(`${task.name}_figure_300dpi.png`, blob, "image/png")); }
    } catch (e) { setError(`Could not make the figure: ${e.message}`); }
    setFigureBusy(false);
  };

  return (
    <div className="card quick">
      <h2>Quick start: drop your files</h2>
      <div className={`drop${dragging ? " over" : ""}`} onDragOver={(e) => { e.preventDefault(); setDragging(true); }} onDragLeave={() => setDragging(false)} onDrop={(e) => { e.preventDefault(); setDragging(false); addFiles(e.dataTransfer.files); }}>
        <p><b>Drop the .ab1 files here</b> (controls and samples together, forward and reverse reads), plus, if you have them, the design (report .html or file .json) and a sample sheet (.csv).</p>
        <input type="file" multiple accept=".ab1,.abi,.csv,.tsv,.json,.html,.htm" onChange={(e) => addFiles(e.target.files)} aria-label="Choose files" />
        <input type="file" multiple webkitdirectory="" directory="" onChange={(e) => addFiles(e.target.files)} aria-label="Choose a folder" />
        <p className="muted">Nothing is uploaded: the files are read and analysed in this browser.</p>
      </div>
      {designName && <p className="muted">Design: {designName}{design.error ? <span className="error"> {design.error}</span> : null}</p>}
      {traceNames.length > 0 && (
        <>
          <p>{traceNames.length} trace{traceNames.length === 1 ? "" : "s"}; {pairs.length} sample{pairs.length === 1 ? "" : "s"} paired by {built.source === "sheet" ? "the sample sheet" : "file names"}. Guides from {guideSource}.</p>
          {built.problems.map((p) => <p key={p} className="note">{p}</p>)}
          {lowCount > 0 && <p className="note">{lowCount} pairing{lowCount === 1 ? " is" : "s are"} a guess (marked low): check the control column before running.</p>}
          <div className="scroll">
            <table><thead><tr><th>Sample</th><th>Reads</th><th>Control</th><th>Type</th><th>Pairing</th></tr></thead><tbody>
              {pairs.map((p, i) => (
                <tr key={p.sample + i}><td>{p.sample}</td><td>{p.edited.length}</td>
                  <td><select value={p.control[0]} onChange={(e) => setPairsEdit({ ...pairsEdit, [i]: { control: [e.target.value], confidence: "set by you" } })}>{traceNames.map((n) => <option key={n} value={n}>{n}</option>)}</select></td>
                  <td><select value={p.type} onChange={(e) => setPairsEdit({ ...pairsEdit, [i]: { ...(pairsEdit[i] || {}), type: e.target.value } })}><option value="clone">clone</option><option value="pool">pool</option></select></td>
                  <td className={p.confidence === "low" ? "warn" : ""} title={p.reason}>{p.confidence}</td></tr>
              ))}
            </tbody></table>
          </div>
          <details><summary>Guides, donor, nuclease, goal (optional)</summary>
            <label className="field" htmlFor="q-guides">Guide sequence(s) if there is no design</label>
            <input id="q-guides" type="text" value={guidesText} onChange={(e) => setGuidesText(e.target.value)} placeholder="leave empty to infer the cut from the traces" />
            <label className="field" htmlFor="q-donor">Donor (ssODN, up to 300 nt) for SNP corrections and small tags</label>
            <textarea id="q-donor" value={donorText} onChange={(e) => setDonorText(e.target.value)} />
            <label className="field" htmlFor="q-nuc">Nuclease</label>
            <select id="q-nuc" value={nuclease} onChange={(e) => setNuclease(e.target.value)}>{Object.entries(NUCLEASES).map(([k, v]) => <option key={k} value={k}>{v.label}</option>)}</select>
            <label className="field" htmlFor="q-type">Sample type (overrides the guess for every sample)</label>
            <select id="q-type" value={typeOverride} onChange={(e) => setTypeOverride(e.target.value)}><option value="">as paired</option><option value="clone">clones</option><option value="pool">pools</option></select>
            <label className="field" htmlFor="q-goal">Accept</label>
            <select id="q-goal" value={goal} onChange={(e) => setGoal(e.target.value)}><option value="homozygous">homozygous edits only (heterozygous clones are held)</option><option value="any">heterozygous edits too</option></select>
          </details>
          <p><button type="button" className="primary" disabled={running || !pairs.length} onClick={run}>{running ? "Analysing..." : `Analyse ${pairs.length} sample${pairs.length === 1 ? "" : "s"}`}</button> {progress && <span className="muted">{progress.done} of {progress.total}</span>}</p>
        </>
      )}
      {error && <p className="error">{error}</p>}
      {outcome && (
        <>
          <h2>Decisions</h2>
          <p className="chips">{DECISION_ORDER.map((k) => <span key={k} className={`chip d-${k.replace("-", "")}`}>{LABEL[k]} {counts[k]}</span>)}</p>
          <p>
            <button type="button" className="link" onClick={() => download(`assured_qc_${stamp}.html`, buildQcReportHtml({ workflow: outcome.workflow, results: outcome.results, notes: outcome.problems }), "text/html")}>Report (HTML)</button>
            <button type="button" className="link" onClick={() => download(`assured_qc_${stamp}.csv`, buildResultsCsv(outcome.workflow, outcome.results), "text/csv")}>Table (CSV)</button>
            <button type="button" className="link" onClick={() => download(`assured_qc_plate_${stamp}.html`, buildPlateReportHtml({ title: "Plate", results: outcome.results }), "text/html")}>Plate view (HTML)</button>
            <button type="button" className="link" onClick={() => download(`assured_qc_resequence_${stamp}.csv`, resequenceListCsv(outcome.results), "text/csv")}>Re-sequence list (CSV)</button>
          </p>
          <div className="scroll">
            <table><thead><tr><th>Decision</th><th>Sample</th><th>Result</th><th>Edited, 95% interval</th><th>Confidence</th></tr></thead><tbody>
              {outcome.results.map((r, i) => (
                <tr key={r.name + i} className={`pick${i === picked ? " active" : ""}`} onClick={() => setPicked(i)}>
                  <td><span className={`chip d-${r.decision.decision.replace("-", "")}`}>{LABEL[r.decision.decision]}</span></td><td>{r.name}</td>
                  <td>{r.ok ? headline(outcome.workflow, r) : r.error}</td><td>{r.ok ? `${r.summary.editedPct}% (${interval(r.intervals?.editedPct)})` : ""}</td><td>{r.confidence ? r.confidence.tier : ""}</td></tr>
              ))}
            </tbody></table>
          </div>
          {result && <div className="card"><p><b>{LABEL[result.decision.decision]}.</b> {result.decision.reasons.join(" ")} <span className="muted">{result.decision.next}</span></p>
            <p><button type="button" className="link" disabled={figureBusy || !result.ok} onClick={() => figure("svg")}>Figure (SVG)</button><button type="button" className="link" disabled={figureBusy || !result.ok} onClick={() => figure("png")}>Figure (PNG, 300 dpi)</button></p></div>}
          {result && <Detail workflow={outcome.workflow} result={result} />}
        </>
      )}
    </div>
  );
}
