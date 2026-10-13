// Import of results exported from the ICE web tool (EditCo ICE, ice.editco.bio), for the comparison in iceCompare.js.
// Nothing is run here: the export is read, its samples are matched to the samples of the run by name, and the scores and alleles
// ICE reported are put in the form that compareRun takes. The layout read is the one of an ICE 4.2 export: a zip with
// summary_<id>.csv (Label, ICE, KO-Score, KI-Score, R Squared, Guide Sequences, Notes, Indels, ...) and one folder per sample.
// Only the summary is used, so the summary .csv alone is also accepted. Works in the browser and in Node (zip entries are inflated with
// the platform's DecompressionStream).
import { parseCsv } from "./sheet.js";
import { compareRun } from "./iceCompare.js";

const EDIT_WORKFLOWS = new Set(["snp", "tag_small"]);
const u16 = (v, o) => v[o] | (v[o + 1] << 8);
const u32 = (v, o) => (v[o] | (v[o + 1] << 8) | (v[o + 2] << 16) | (v[o + 3] << 24)) >>> 0;
const bytesOf = (b) => (b instanceof Uint8Array ? b : new Uint8Array(b));

async function inflateRaw(data) {
  const stream = new Blob([data]).stream().pipeThrough(new DecompressionStream("deflate-raw"));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

/** Entries of a zip file: Map(name -> Uint8Array), for the names that `accept` takes. Stored and deflated entries; no zip64, no encryption. */
export async function readZip(input, accept = () => true) {
  const v = bytesOf(input); let eocd = -1;
  for (let i = v.length - 22; i >= Math.max(0, v.length - 22 - 65535); i -= 1) if (u32(v, i) === 0x06054b50) { eocd = i; break; }
  if (eocd < 0) throw new Error("This is not a zip file.");
  const count = u16(v, eocd + 10); let p = u32(v, eocd + 16); const out = new Map();
  for (let k = 0; k < count; k += 1) {
    if (p + 46 > v.length || u32(v, p) !== 0x02014b50) throw new Error("The zip file is damaged (its directory cannot be read).");
    const flags = u16(v, p + 8), method = u16(v, p + 10), csize = u32(v, p + 20), usize = u32(v, p + 24), nlen = u16(v, p + 28), elen = u16(v, p + 30), clen = u16(v, p + 32), offset = u32(v, p + 42);
    const name = new TextDecoder().decode(v.subarray(p + 46, p + 46 + nlen)); p += 46 + nlen + elen + clen;
    if (name.endsWith("/") || !accept(name)) continue;
    if (flags & 1) throw new Error("The zip file is encrypted.");
    if (csize === 0xffffffff || usize === 0xffffffff || offset === 0xffffffff) throw new Error("Zip64 files are not supported.");
    if (u32(v, offset) !== 0x04034b50) throw new Error(`The zip file is damaged (${name}).`);
    const start = offset + 30 + u16(v, offset + 26) + u16(v, offset + 28); const data = v.subarray(start, start + csize);
    if (method === 0) out.set(name, data.slice());
    else if (method === 8) out.set(name, await inflateRaw(data));
    else throw new Error(`Zip compression method ${method} is not supported (${name}).`);
  }
  return out;
}

const key = (h) => String(h).toLowerCase().replace(/[^a-z0-9]+/g, "");
const COLUMNS = { label: ["label", "sample", "samplename", "name", "experiment"], ice: ["ice", "indel", "indelpct"], ko: ["koscore", "ko"], ki: ["kiscore", "ki"], rsq: ["rsquared", "rsq", "r2"], guides: ["guidesequences", "guidesequence", "guides", "guide"], notes: ["notes", "note"], indels: ["indels"] };
const find = (header, names) => header.findIndex((h) => names.includes(key(h)));
const number = (v) => { const t = String(v ?? "").trim(); if (!t || /^(none|nan|na|n\/a|null|-)$/i.test(t)) return null; const n = Number(t); return Number.isFinite(n) ? n : null; };
const dna = (s) => String(s || "").toUpperCase().replace(/U/g, "T").replace(/[^ACGT]/g, "");

/** True if the text starts like an ICE summary table (a Label column and a KO-Score, KI-Score or ICE column). */
export function looksLikeIceSummary(text) {
  const first = String(text || "").replace(/^\uFEFF/, "").split(/\r?\n/, 1)[0] || ""; const header = first.split(",").map(key);
  return header.includes("label") && (header.includes("koscore") || header.includes("kiscore")) && header.includes("rsquared");
}

/** The Indels cell of the summary, a Python dict such as {'-77': 49.0, '-79': 51.0} (size: percent), as [{ size, pct }]; sizes that are not numbers are left out. */
export function parseIndels(cellText) {
  const out = []; const re = /['"]?([+-]?\d+)['"]?\s*:\s*(-?\d+(?:\.\d+)?)/g; let m;
  while ((m = re.exec(String(cellText || "")))) out.push({ size: Number(m[1]), pct: Number(m[2]) });
  return out;
}

/** Read an ICE summary table. Returns { rows: [{ label, ice, ko, ki, rsq, guides, notes, indels }], problems }. */
export function parseIceSummary(text) {
  const table = parseCsv(text); const problems = [];
  if (table.length < 2) return { rows: [], problems: ["The ICE summary has no samples."] };
  const header = table[0]; const at = Object.fromEntries(Object.entries(COLUMNS).map(([k, names]) => [k, find(header, names)]));
  if (at.label < 0) return { rows: [], problems: ["The ICE summary has no Label column."] };
  if (at.ko < 0 && at.ki < 0 && at.ice < 0) return { rows: [], problems: ["The ICE summary has no ICE, KO-Score or KI-Score column."] };
  const cellOf = (r, k) => (at[k] >= 0 ? r[at[k]] : "");
  const rows = table.slice(1).map((r) => ({ label: String(cellOf(r, "label")).trim(), ice: number(cellOf(r, "ice")), ko: number(cellOf(r, "ko")), ki: number(cellOf(r, "ki")), rsq: number(cellOf(r, "rsq")), guides: String(cellOf(r, "guides")).split(/[,;\s]+/).map(dna).filter(Boolean), notes: String(cellOf(r, "notes")).trim(), indels: parseIndels(cellOf(r, "indels")) })).filter((r) => r.label);
  if (!rows.length) problems.push("The ICE summary has no samples.");
  return { rows, problems };
}

const isSummaryEntry = (name) => /(^|\/)summary[^/]*\.csv$/i.test(name);
const decode = (b) => new TextDecoder("utf-8").decode(bytesOf(b));

/**
 * Read what the user gives for the ICE results: the zip of the ICE export, its summary .csv, or several of them.
 * @param {Array<{name: string, bytes?: ArrayBuffer|Uint8Array, text?: string}>} inputs
 * @returns {Promise<{rows: Array, sources: string[], problems: string[]}>}
 */
export async function readIceExport(inputs) {
  const rows = [], sources = [], problems = [];
  for (const input of inputs) {
    try {
      let text = input.text;
      if (text === undefined) {
        const bytes = bytesOf(input.bytes);
        if (/\.zip$/i.test(input.name) || (bytes[0] === 0x50 && bytes[1] === 0x4b)) {
          const entries = await readZip(bytes, isSummaryEntry);
          const names = [...entries.keys()].sort((a, b) => a.split("/").length - b.split("/").length);
          if (!names.length) { problems.push(`${input.name}: no summary .csv in the zip (the ICE export has one at its top level).`); continue; }
          text = decode(entries.get(names[0]));
        } else text = decode(bytes);
      }
      const parsed = parseIceSummary(text);
      if (parsed.problems.length) { problems.push(...parsed.problems.map((p) => `${input.name}: ${p}`)); continue; }
      rows.push(...parsed.rows); sources.push(input.name);
    } catch (e) { problems.push(`${input.name}: ${e.message || e}`); }
  }
  return { rows, sources, problems };
}

const WORDS = (s) => String(s ?? "").replace(/\.(ab1|abi|abif|scf)$/i, "").split(/[^A-Za-z0-9]+/).filter(Boolean);
const at = (hay, needle) => { for (let i = 0; i + needle.length <= hay.length; i += 1) if (needle.every((t, j) => hay[i + j] === t)) return i; return -1; };

/**
 * Match the labels of the ICE rows to sample names. A label belongs to a sample when the sample's name (split at anything that is not a letter or
 * digit) occurs in it as whole words, so "cl4" is found in "cl4_KO_run2" and not in "cl40"; the longest name wins. A label that fits two samples
 * equally well is not assigned. What remains of the label after the sample's name is the donor name for a point edit or tag.
 * Returns { bySample: Map(position -> [{ row, rest }]), unmatched: rows, ambiguous: [{ row, samples }] }.
 */
export function matchIceRows(names, rows) {
  const keys = names.map((n) => WORDS(n).map((w) => w.toLowerCase())); const bySample = new Map(), unmatched = [], ambiguous = [];
  rows.forEach((row) => {
    const words = WORDS(row.label); const lower = words.map((w) => w.toLowerCase());
    const hits = keys.map((k, i) => ({ i, k, pos: k.length ? at(lower, k) : -1 })).filter((h) => h.pos >= 0);
    if (!hits.length) { unmatched.push(row); return; }
    const longest = Math.max(...hits.map((h) => h.k.length)); const top = hits.filter((h) => h.k.length === longest);
    if (top.length > 1) { ambiguous.push({ row, samples: top.map((h) => names[h.i]) }); return; }
    const h = top[0]; const rest = [...words.slice(0, h.pos), ...words.slice(h.pos + h.k.length)].join(" ");
    if (!bySample.has(h.i)) bySample.set(h.i, []); bySample.get(h.i).push({ row, rest });
  });
  return { bySample, unmatched, ambiguous };
}

const spacerOf = (g) => dna((g && (g.spacer || g.sequence)) || g);
function guideWarning(task, row) {
  const used = [...new Set(row.guides)].sort(); const wanted = [...new Set(((task.spec && task.spec.guides) || []).map(spacerOf).filter(Boolean))].sort();
  if (!used.length || !wanted.length || used.join() === wanted.join()) return "";
  return `ICE was run with ${used.length} guide${used.length === 1 ? "" : "s"} (${used.join(", ")}), not the design's (${wanted.join(", ")}): the scores are not like for like.`;
}
function runOf(row, donorName, warning) {
  if (row.ice === null && row.ko === null && row.ki === null) return { donorName, ok: false, error: row.notes || "ICE gave no score for this sample" };
  const note = [warning, row.notes ? `ICE note: ${row.notes}` : ""].filter(Boolean).join(" ");
  return { donorName, ok: true, summary: { ice: row.ice, ko_score: row.ko, hdr_pct: row.ki, rsq: row.rsq }, contribs: row.indels.map(({ size, pct }) => ({ indel: size, rel: pct / 100 })), warning: note, label: row.label };
}

/**
 * Compare the results of a run with an ICE export.
 * @param {{tasks: Array, results: Array, rows: Array}} given tasks and results as the run made them; rows from parseIceSummary
 * @returns {{rows: Array, unmatched: Array, ambiguous: Array, problems: string[]}} rows as compareRun returns them
 */
export function compareWithExport({ tasks, results, rows }) {
  const names = tasks.map((t) => t.name); const { bySample, unmatched, ambiguous } = matchIceRows(names, rows);
  const byIndex = new Map(), skipped = [];
  tasks.forEach((t, i) => {
    const hits = bySample.get(i);
    if (!hits) { skipped.push({ index: i, reason: "no sample of the ICE export has this sample's name" }); return; }
    byIndex.set(i, hits.map(({ row, rest }) => runOf(row, EDIT_WORKFLOWS.has(t.workflow) ? rest : "", guideWarning(t, row))));
  });
  const problems = [];
  if (unmatched.length) problems.push(`${unmatched.length} label${unmatched.length === 1 ? "" : "s"} of the ICE export ${unmatched.length === 1 ? "matches" : "match"} none of the run's samples: ${unmatched.map((r) => r.label).join(", ")}.`);
  ambiguous.forEach((a) => problems.push(`The ICE label "${a.row.label}" fits more than one sample (${a.samples.join(", ")}) and was not used.`));
  return { rows: compareRun(results, byIndex, skipped), unmatched, ambiguous, problems };
}
