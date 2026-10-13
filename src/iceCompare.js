// Comparison with Synthego ICE, run separately on the same traces (command line: --ice <ICE folder>).
// ICE is a comparator only: it changes no call, rank, score or figure of Assured QC. A comparison is reported as agree, differ, or ICE fit poor
// (ICE's own guidance is that R squared below 0.8 should be read with caution, and a score from a poor fit is not comparable), with the alleles
// that Assured QC reads and ICE does not place. ICE takes one donor per run, so a point edit or small tag is run once per donor and the best donor is shown.

export const ICE_FIT_OK = 0.8;        // R squared at or above which an ICE fit is taken as sound
export const ICE_DIFFERENCE = 15;     // points between the two scores that are called a difference
export const ICE_DONOR_MIN = 20;      // ICE takes a donor of 20 to 300 bases
export const ICE_DONOR_MAX = 300;
const MIN_ALLELE = 0.05;              // alleles below this share are not listed
const KNOCK = new Set(["knockout", "deletion"]);
const EDIT = new Set(["snp", "tag_small"]);
const esc = (v) => String(v ?? "").replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
const num = (v) => (typeof v === "number" && Number.isFinite(v) ? v : null);
const spacerOf = (g) => String((g && (g.spacer || g.sequence)) || g || "").toUpperCase().replace(/[^ACGT]/g, "");
const donorSeq = (d) => String((d && (d.sequence || d.seq)) || d || "").toUpperCase().replace(/[^ACGT]/g, "");

/** The number Assured QC puts against ICE's: the KO score for a knockout or deletion, the share of intended edit for a point edit or tag. */
export function appScore(r) {
  if (!r || !r.ok || !r.summary) return null;
  const s = r.summary, iv = r.intervals || {};
  if (KNOCK.has(r.workflow)) {
    const value = Math.max(0, Math.round(((s.koScorePct ?? 0) * (100 - (s.wtPct ?? 0))) / 100));
    return { value, interval: iv.koScorePct || null, label: "KO score" };
  }
  if (EDIT.has(r.workflow)) return { value: Math.round(s.intendedEditPct ?? 0), interval: iv.intendedEditPct || null, label: "HDR score" };
  return null;
}

/** Alleles Assured QC reads (not wild type), largest first: [{ size, share }]. */
export function appAlleles(r) {
  const bySize = new Map(); // a deletion that can slide along a repeat is listed twice by the library: one allele per size
  ((r && r.contributions) || []).filter((c) => c.kind !== "wt" && Number.isFinite(c.size) && c.size !== 0).forEach((c) => bySize.set(c.size, (bySize.get(c.size) || 0) + (c.fraction ?? 0)));
  return [...bySize].filter(([, share]) => share >= MIN_ALLELE).map(([size, share]) => ({ size, share })).sort((a, b) => b.share - a.share);
}

/** One ICE job per sample (knockout, deletion) or per sample and donor (point edit, small tag). A sample ICE cannot take is listed with a reason. */
export function iceJobsFor(tasks, { nuclease = "SpCas9" } = {}) {
  const jobs = [], skipped = [];
  tasks.forEach((t, index) => {
    const spec = t.spec || {}; const guides = (spec.guides || []).map(spacerOf);
    const why = (reason) => skipped.push({ sample: t.name, index, reason });
    if (!KNOCK.has(t.workflow) && !EDIT.has(t.workflow)) return why(`ICE is not run for the ${t.workflow} workflow`);
    if (String(t.nuclease || nuclease) !== "SpCas9") return why("ICE supports SpCas9 only");
    if (guides.length < 1 || guides.length > 3 || guides.some((g) => g.length < 17 || g.length > 23)) return why("ICE needs one to three guides of 17 to 23 bases without the PAM");
    const base = { sample: t.name, index, workflow: t.workflow, control: t.controls[0], edited: t.edited[0], guides };
    if (KNOCK.has(t.workflow)) { jobs.push({ ...base, donor: "", donorName: "" }); return; }
    const donors = (spec.donors || []).map((d, i) => ({ name: d.name || `donor ${i + 1}`, sequence: donorSeq(d) })).filter((d) => d.sequence);
    const usable = donors.filter((d) => d.sequence.length >= ICE_DONOR_MIN && d.sequence.length <= ICE_DONOR_MAX);
    if (!usable.length) return why(donors.length ? `ICE takes a donor of ${ICE_DONOR_MIN} to ${ICE_DONOR_MAX} bases` : "no donor sequence in the design");
    usable.forEach((d) => jobs.push({ ...base, donor: d.sequence, donorName: d.name }));
  });
  return { jobs, skipped };
}

const STATUS_TEXT = { agree: "Agree", differ: "Differ: check the trace", ice_fit_poor: "ICE fit poor", ice_failed: "ICE failed", not_run: "Not run" };
export const iceStatusText = (status) => STATUS_TEXT[status] || status;
const sizeList = (items) => items.map((a) => `${a.size > 0 ? "+" : ""}${a.size} (${Math.round(a.share * 100)}%)`).join(", ");

/**
 * Compare one analysed sample with its ICE results.
 * @param {object} r Assured QC result of the sample
 * @param {Array} iceRuns one entry per ICE run of the sample: { donorName, ok, summary: { ice, rsq, ko_score, hdr_pct }, contribs: [{ rel, indel }], error }
 * @param {string} skipReason why ICE was not run, if it was not
 */
export function compareWithIce(r, iceRuns, skipReason = "") {
  const app = appScore(r);
  const row = { sample: r.name, workflow: r.workflow, appLabel: app ? app.label : "", app: app ? app.value : null, appInterval: app ? app.interval : null, ice: null, iceR2: null, iceIndel: null, donor: "", difference: null, status: "not_run", appAlleles: [], iceAlleles: [], missing: [], note: skipReason };
  if (!app || !iceRuns || !iceRuns.length) return row;
  const good = iceRuns.filter((x) => x && x.ok && x.summary);
  if (!good.length) { row.status = "ice_failed"; row.note = (iceRuns.find((x) => x && x.error) || {}).error || "ICE returned no result"; return row; }
  const knock = KNOCK.has(r.workflow);
  const valueOf = (x) => num(knock ? x.summary.ko_score : x.summary.hdr_pct) ?? -1;
  const best = good.slice().sort((a, b) => valueOf(b) - valueOf(a))[0];
  row.ice = valueOf(best) < 0 ? null : Math.round(valueOf(best)); row.iceR2 = num(best.summary.rsq); row.iceIndel = num(best.summary.ice); row.donor = knock ? "" : best.donorName || "";
  row.appAlleles = appAlleles(r);
  row.iceAlleles = (best.contribs || []).filter((c) => Number.isFinite(c.indel) && c.indel !== 0 && (c.rel ?? 0) >= MIN_ALLELE).map((c) => ({ size: c.indel, share: c.rel })).sort((a, b) => b.share - a.share);
  if (knock) row.missing = row.appAlleles.filter((a) => a.share >= 0.1 && !row.iceAlleles.some((c) => c.size === a.size));
  if (row.ice === null) { row.status = "ice_failed"; row.note = "ICE gave no score"; return row; }
  row.difference = row.app - row.ice;
  const poor = row.iceR2 !== null && row.iceR2 < ICE_FIT_OK;
  const close = Math.abs(row.difference) < ICE_DIFFERENCE;
  if (poor) {
    row.status = "ice_fit_poor";
    row.note = `ICE explains the trace poorly (R\u00b2 ${row.iceR2.toFixed(2)}, below ${ICE_FIT_OK}), so its score is not comparable${close ? "; the scores are within " + ICE_DIFFERENCE + " points anyway" : ""}.` + (row.missing.length ? ` ICE places ${row.iceAlleles.length ? sizeList(row.iceAlleles) : "no allele"}; Assured QC also reads ${sizeList(row.missing)}, which ICE does not.` : "");
  } else if (close) { row.status = "agree"; row.note = ""; }
  else { row.status = "differ"; row.note = `The scores differ by ${Math.abs(row.difference)} points with a sound ICE fit (R\u00b2 ${row.iceR2 === null ? "n/a" : row.iceR2.toFixed(2)}): look at the trace before relying on either.` + (row.missing.length ? ` Assured QC reads ${sizeList(row.missing)}, which ICE does not place.` : ""); }
  if (!knock && row.donor && good.length > 1) row.note = (row.note ? row.note + " " : "") + `Best of ${good.length} donors run one at a time (${good.map((x) => `${x.donorName}: ${Math.round(valueOf(x))}`).join(", ")}).`;
  return row;
}

/** JSON.parse that also accepts the NaN and Infinity that Python's json.dump writes for non-finite numbers (read as null). */
export const parseIceJson = (text) => JSON.parse(String(text).replace(/(?<=[:\[,]\s*)-?(?:NaN|Infinity)(?=\s*[,}\]])/g, "null"));
const baseName = (p) => String(p).split(/[\\/]/).pop();
/** The record kept for every ICE run (ice_raw.json): what was run, on which traces (SHA-256), and what ICE returned. `shaOf(path)` gives the hash. */
export function iceRawRecords(jobs, results, shaOf) {
  return jobs.map((j, k) => ({ sample: j.sample, donorName: j.donorName || "", guides: j.guides, donor: j.donor || "", control: baseName(j.control), controlSha256: shaOf(j.control), edited: baseName(j.edited), editedSha256: shaOf(j.edited), result: results.get(k) || { ok: false, error: "not run" } }));
}
/** Reuse earlier ICE results (--ice-raw) only for the same sample, guides, donor and traces. Returns { results: Map(job position -> result), problems }. */
export function matchIceRaw(jobs, records, shaOf) {
  const results = new Map(), problems = [];
  jobs.forEach((j, k) => {
    const label = `${j.sample}${j.donorName ? ` (${j.donorName})` : ""}`;
    const rec = (records || []).find((r) => r.sample === j.sample && (r.donorName || "") === (j.donorName || ""));
    if (!rec) return problems.push(`no earlier ICE result for ${label}`);
    if ((rec.guides || []).join(",") !== j.guides.join(",") || (rec.donor || "") !== (j.donor || "")) return problems.push(`${label}: the earlier ICE result used other guides or another donor`);
    if (rec.controlSha256 !== shaOf(j.control) || rec.editedSha256 !== shaOf(j.edited)) return problems.push(`${label}: the earlier ICE result was made on other traces`);
    results.set(k, rec.result);
  });
  return { results, problems };
}

/** Rows for a whole run. `iceByIndex` maps a sample's index (position in `tasks`) to its ICE runs; `skipped` is from iceJobsFor. */
export function compareRun(results, iceByIndex, skipped = []) {
  const skip = new Map(skipped.map((s) => [s.index, s.reason]));
  return results.map((r, i) => (r && r.ok && appScore(r) ? compareWithIce(r, iceByIndex.get(i), skip.get(i) || "") : null)).filter(Boolean);
}

export function iceSummary(rows) {
  const count = (s) => rows.filter((r) => r.status === s).length;
  return { total: rows.length, agree: count("agree"), differ: count("differ"), ice_fit_poor: count("ice_fit_poor"), ice_failed: count("ice_failed"), not_run: count("not_run") };
}

const cell = (v) => { const t = String(v ?? ""); return /[",\n]/.test(t) ? `"${t.replace(/"/g, '""')}"` : t; };
export function iceCsv(rows) {
  const header = ["sample", "workflow", "assured_qc_score", "assured_qc_interval", "ice_score", "ice_indel_pct", "ice_r2", "difference", "status", "donor_used", "assured_qc_alleles", "ice_alleles", "alleles_not_found_by_ice", "note"];
  const lines = rows.map((r) => [r.sample, r.workflow, r.app ?? "", r.appInterval ? `${r.appInterval[0]}-${r.appInterval[1]}` : "", r.ice ?? "", r.iceIndel ?? "", r.iceR2 ?? "", r.difference ?? "", r.status, r.donor, sizeList(r.appAlleles), sizeList(r.iceAlleles), sizeList(r.missing), r.note].map(cell).join(","));
  return [header.join(","), ...lines].join("\n") + "\n";
}

const COLOR = { agree: "#1f6fae", differ: "#c0392b", ice_fit_poor: "#d98a1c" };
/** Assured QC score against the ICE score, one point per clone. Clones that agree are counted, not labelled; the others are named. */
export function iceScatterSvg(rows, { title = "Assured QC against ICE" } = {}) {
  const pts = rows.filter((r) => r.app !== null && r.ice !== null); if (!pts.length) return "";
  const W = 560, H = 580, L = 70, T = 52, S = 400; const x = (v) => L + (v / 100) * S, y = (v) => T + S - (v / 100) * S;
  const knock = pts.every((r) => KNOCK.has(r.workflow)); const name = knock ? "KO score" : "HDR score";
  const out = [`<svg xmlns="http://www.w3.org/2000/svg" width="${W}" height="${H}" viewBox="0 0 ${W} ${H}" font-family="Arial, Helvetica, sans-serif" role="img" aria-label="${esc(title)}"><title>${esc(title)}</title><rect width="${W}" height="${H}" fill="#fff"/>`];
  out.push(`<text x="14" y="26" font-size="15" font-weight="bold" fill="#101828">${esc(title)}</text>`);
  const D = ICE_DIFFERENCE; const band = [[0, 0], [0, D], [100 - D, 100], [100, 100], [100, 100 - D], [D, 0]].map(([a, b]) => `${x(a)},${y(b)}`).join(" ");
  out.push(`<polygon points="${band}" fill="#eef2f7"/>`);
  out.push(`<rect x="${L}" y="${T}" width="${S}" height="${S}" fill="none" stroke="#667085"/>`);
  for (const t of [0, 25, 50, 75, 100]) {
    out.push(`<line x1="${x(t)}" y1="${T + S}" x2="${x(t)}" y2="${T + S + 5}" stroke="#667085"/><text x="${x(t)}" y="${T + S + 19}" font-size="11" text-anchor="middle" fill="#344054">${t}</text>`);
    out.push(`<line x1="${L - 5}" y1="${y(t)}" x2="${L}" y2="${y(t)}" stroke="#667085"/><text x="${L - 9}" y="${y(t) + 4}" font-size="11" text-anchor="end" fill="#344054">${t}</text>`);
  }
  out.push(`<line x1="${x(0)}" y1="${y(0)}" x2="${x(100)}" y2="${y(100)}" stroke="#667085" stroke-dasharray="4 3"/>`);
  out.push(`<text x="${L + S / 2}" y="${T + S + 40}" font-size="12" text-anchor="middle" fill="#101828">Assured QC ${name}</text>`);
  out.push(`<text x="20" y="${T + S / 2}" font-size="12" text-anchor="middle" fill="#101828" transform="rotate(-90 20 ${T + S / 2})">ICE ${name}</text>`);
  pts.forEach((r) => { const c = COLOR[r.status] || "#667085"; const hollow = r.status === "ice_fit_poor";
    out.push(`<circle cx="${x(r.app).toFixed(1)}" cy="${y(r.ice).toFixed(1)}" r="5.5" fill="${hollow ? "#fff" : c}" fill-opacity="${hollow ? 1 : 0.8}" stroke="${c}" stroke-width="${hollow ? 2 : 1}"><title>${esc(`${r.sample}: Assured QC ${r.app}, ICE ${r.ice}${r.iceR2 !== null ? `, R\u00b2 ${r.iceR2.toFixed(2)}` : ""}`)}</title></circle>`); });
  // name the clones that do not agree, stacked so the labels do not overlap
  const named = pts.filter((r) => r.status !== "agree").sort((a, b) => a.ice - b.ice); let lastY = -99;
  named.forEach((r) => { let ty = y(r.ice) + 4; if (ty - lastY < 13) ty = lastY + 13; lastY = ty; const c = COLOR[r.status] || "#667085"; const tx = Math.min(x(r.app) + 10, L + S - 70);
    out.push(`<text x="${tx.toFixed(1)}" y="${ty.toFixed(1)}" font-size="11" fill="${c}">${esc(r.sample)}${r.iceR2 !== null ? ` (R\u00b2 ${r.iceR2.toFixed(2)})` : ""}</text>`); });
  const n = (s) => pts.filter((r) => r.status === s).length; const ly = T + S + 70;
  out.push(`<circle cx="${L + 6}" cy="${ly}" r="5" fill="${COLOR.agree}" fill-opacity="0.8"/><text x="${L + 16}" y="${ly + 4}" font-size="11" fill="#344054">agree, ${n("agree")} (within ${ICE_DIFFERENCE} points)</text>`);
  out.push(`<circle cx="${L + 190}" cy="${ly}" r="5" fill="#fff" stroke="${COLOR.ice_fit_poor}" stroke-width="2"/><text x="${L + 200}" y="${ly + 4}" font-size="11" fill="#344054">ICE fit poor, ${n("ice_fit_poor")} (R\u00b2 below ${ICE_FIT_OK})</text>`);
  out.push(`<circle cx="${L + 380}" cy="${ly}" r="5" fill="${COLOR.differ}" fill-opacity="0.8"/><text x="${L + 390}" y="${ly + 4}" font-size="11" fill="#344054">differ, ${n("differ")}</text>`);
  out.push(`<text x="14" y="${H - 12}" font-size="10" fill="#667085">Dashed line: equal scores; shaded band: within ${ICE_DIFFERENCE} points. ICE is a comparator and changes no call or rank.</text></svg>`);
  return out.join("");
}

/** The report section: the plot and one row per sample. `meta` is the provenance recorded for the run. */
export function iceSectionHtml(rows, meta = {}) {
  if (!rows.length) return "";
  const s = iceSummary(rows); const svg = iceScatterSvg(rows);
  const head = `<h2>Comparison with ICE</h2><p class="muted">ICE ${esc(meta.version || "")} run separately on the same traces${meta.python ? ` (Python ${esc(meta.python)})` : ""}. It changes no call, rank or score: ${s.agree} agree, ${s.ice_fit_poor} with a poor ICE fit, ${s.differ} differ${s.ice_failed ? `, ${s.ice_failed} failed` : ""}${s.not_run ? `, ${s.not_run} not run` : ""}. A score from an ICE fit with R\u00b2 below ${ICE_FIT_OK} is not comparable; the alleles ICE does not place are named. For a point edit or tag ICE is run once per donor and the best donor is shown. ICE is used under its research and non-profit licence.</p>`;
  const body = rows.map((r) => `<tr><td>${esc(r.sample)}</td><td>${r.app ?? ""}${r.appInterval ? ` (${r.appInterval[0]}&ndash;${r.appInterval[1]})` : ""}</td><td>${r.ice ?? ""}</td><td>${r.iceR2 === null ? "" : r.iceR2.toFixed(2)}</td><td>${esc(iceStatusText(r.status))}</td><td>${esc(r.note)}</td></tr>`).join("");
  return `${head}${svg ? `<div style="overflow-x:auto">${svg}</div>` : ""}<table><tr><th>Sample</th><th>Assured QC (95% interval)</th><th>ICE</th><th>ICE fit R&sup2;</th><th>Status</th><th>Note</th></tr>${body}</table>`;
}
