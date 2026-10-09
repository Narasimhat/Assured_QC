import { headline } from "./batch.js";

const pct = (value) => (value === null || value === undefined ? "n/a" : `${(value * 100).toFixed(1)}%`);

export default function Detail({ workflow, result }) {
  if (!result) return null;
  if (!result.ok) return <div className="card"><h2>{result.name}</h2><p className="error">{result.error}</p></div>;
  if (result.kind === "junction") {
    return (
      <div className="card">
        <h2>{result.name}</h2>
        <p className={result.verdict.status}>{result.verdict.summary}</p>
        <ul>{result.verdict.details.map((d) => <li key={d}>{d}</li>)}</ul>
        <p className="muted">Orientation {result.orientation}; {result.alignedBases} aligned bases at {(result.identity * 100).toFixed(1)}% identity. 5′ junction spanned: {String(result.spans.fivePrimeJunction)}; 3′ junction spanned: {String(result.spans.threePrimeJunction)}; genomic flank bases: {result.spans.genomicFlankBases}.</p>
        {result.events.length > 0 && (
          <table>
            <thead><tr><th>Difference</th><th>Region</th><th>Allele position</th><th>Detail</th><th>Quality</th></tr></thead>
            <tbody>{result.events.map((e, i) => <tr key={i}><td>{e.type}{e.expectedChange ? " (designed change)" : ""}</td><td>{e.region}</td><td>{e.alleleIndex}</td><td>{e.detail}</td><td>{e.quality}</td></tr>)}</tbody>
          </table>
        )}
        {result.markers.some((m) => m.covered) && (
          <table>
            <thead><tr><th>Position</th><th>Change</th><th>Read shows</th></tr></thead>
            <tbody>{result.markers.map((m) => <tr key={m.position}><td>{m.position}</td><td>{m.change} {m.label}</td><td>{m.covered ? (m.present ? `designed base (${m.base})` : `not the designed base (${m.base})`) : "not covered by this read"}</td></tr>)}</tbody>
          </table>
        )}
        {result.bands && <p><b>Out-out PCR:</b> {result.bands.call}. {result.bands.summary}</p>}
      </div>
    );
  }
  const s = result.summary;
  return (
    <div className="card">
      <h2>{result.name}</h2>
      <p>{headline(workflow, result)}</p>
      {result.genotype && <ul>{result.genotype.flags.map((f) => <li key={f} className="muted">{f}</li>)}</ul>}
      <table>
        <thead><tr><th>Wild type</th><th>Edited</th><th>Indels</th><th>Between cuts</th><th>Intended edit</th><th>Blocking-only</th><th>KO score</th><th>Unexplained</th></tr></thead>
        <tbody><tr><td>{s.wtPct}%</td><td>{s.editedPct}%</td><td>{s.indelPct}%</td><td>{s.betweenCutsPct}%</td><td>{s.intendedEditPct}%</td><td>{s.partialConversionPct}%</td><td>{s.koScorePct}%</td><td>{s.unexplainedPct ?? 0}%</td></tr></tbody>
      </table>
      <div className="chart" dangerouslySetInnerHTML={{ __html: result.charts.contributions }} />
      <div className="chart" dangerouslySetInnerHTML={{ __html: result.charts.discordance }} />
      {(workflow === "snp" || workflow === "tag_small") && result.markers.some((m) => m.covered) && (
        <table>
          <thead><tr><th>Position</th><th>Change</th><th>Role</th><th>Donor base in edited</th><th>in control</th><th>Net</th></tr></thead>
          <tbody>{result.markers.map((m) => (
            <tr key={m.position}><td>{m.position}</td><td>{m.change} {m.label}</td><td>{m.role}</td>{m.covered ? <><td>{pct(m.altEdited)}</td><td>{pct(m.altControl)}</td><td>{pct(m.altNet)}</td></> : <td colSpan={3} className="muted">{m.reason}</td>}</tr>
          ))}</tbody>
        </table>
      )}
      <p className="muted">Fit R² {result.quality.r2 ?? "n/a"}; background {pct(result.quality.noiseFraction)}; {result.quality.positionsUsed} positions in the window.</p>
      {result.warnings.length > 0 && <ul>{result.warnings.map((w) => <li key={w} className="warn">{w}</li>)}</ul>}
    </div>
  );
}
