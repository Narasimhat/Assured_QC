import test from "node:test";
import assert from "node:assert/strict";
import { inflateSync } from "node:zlib";
import { convertSvg, svgToPdf, parseSvg } from "../src/svgPdf.js";
import { analyseSample } from "../src/batch.js";
import { cloneAlignedFigure, clonesAlignedFigure, clonesAlignmentFigure, alignedGroups } from "../src/alignedFigure.js";
import { cloneSelectionFigure } from "../src/selectionFigure.js";
import { runSummaryFigure } from "../src/summaryFigure.js";
import { drawProteinProducts } from "../src/proteinFigure.js";
import { withRecommendations } from "../src/recommend.js";
import { loadSpec, loadTruth, loadTrace } from "./helpers.js";

const latin1 = (bytes) => Buffer.from(bytes).toString("latin1");
const wrap = (body, w = 200, h = 100) => `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}" viewBox="0 0 ${w} ${h}" font-family="Arial, Helvetica, sans-serif">${body}</svg>`;
const plain = async (svg) => latin1((await convertSvg(svg, { compress: false })).bytes);

// real figures of the app, drawn from the simulated traces of the test fixtures
const truth = loadTruth();
const analyse = (key, specName) => { const t = truth[key]; return analyseSample({ name: key, spec: loadSpec(specName), control: loadTrace(t.control), edited: loadTrace(t.edited), sampleType: "clone" }); };
const decided = (r, decision, name) => ({ ...r, name, workflow: "snp", decision: { decision, reasons: [], next: "" } });
const ranked = withRecommendations([decided(analyse("apoe_snp_clone_hom", "apoe_r176c_snp"), "accept", "clone_hom"), decided(analyse("apoe_snp_clone_het", "apoe_r176c_snp"), "hold", "clone_het")]);
const group = alignedGroups(ranked)[0];
const protein = drawProteinProducts({ proteinLength: 169, regions: [{ name: "Domain", start: 20, end: 60 }], exonEnds: [34, 84, 137, 169], epitope: null,
  products: [{ kind: "truncated", allele: "-78 bp", share: 0.5, kept: 40, total: 120, nmd: "expected", outcome: "frameshift" }, { kind: "inframe", allele: "-84 bp", share: 1, kept: 40, total: 141, removed: 28, outcome: "in-frame" }] }, { title: "Products" });
const figures = {
  "clone selection": cloneSelectionFigure(ranked), "summary graph": runSummaryFigure(ranked), "one clone against the control": cloneAlignedFigure(ranked[0].view, { name: "clone_hom", decision: "accept", headline: "Homozygous for the intended edit" }),
  "all clones stacked": clonesAlignedFigure(group), "all clones as an alignment": clonesAlignmentFigure(group), "protein products (tiling pattern)": protein,
};

for (const [name, svg] of Object.entries(figures)) {
  test(`PDF of the ${name}: valid, same size, every text drawn, nothing skipped`, async () => {
    assert.ok(svg && svg.startsWith("<svg"), "the figure exists");
    const { bytes, lost, unsupported, width, height } = await convertSvg(svg, { compress: false });
    const text = latin1(bytes);
    assert.ok(text.startsWith("%PDF-1.4") && text.trimEnd().endsWith("%%EOF"));
    assert.deepEqual(unsupported, [], "every element of the app's own figures is supported");
    assert.equal(lost, 0, "every character maps to the PDF fonts");
    const w = Number(/width="([\d.]+)"/.exec(svg)[1]), h = Number(/height="([\d.]+)"/.exec(svg)[1]);
    assert.equal(width, w * 0.75); assert.equal(height, h * 0.75); assert.match(text, new RegExp(`/MediaBox \\[0 0 ${width} ${height}\\]`));
    // every <text> with content becomes one Tj
    const texts = (svg.match(/<text\b[^>]*>(?:(?!<\/text>)[\s\S])*<\/text>/g) || []).filter((t) => t.replace(/<title>[\s\S]*?<\/title>/g, "").replace(/<[^>]+>/g, "").trim());
    assert.equal((text.match(/\) Tj\n/g) || []).length, texts.length);
    // the cross-reference table points at every object
    const xref = Number(/startxref\n(\d+)\n/.exec(text)[1]); assert.ok(text.slice(xref).startsWith("xref\n0 "));
    const count = Number(/^xref\n0 (\d+)\n/.exec(text.slice(xref))[1]); const rows = text.slice(xref).split("\n").slice(2, 2 + count);
    rows.slice(1).forEach((row, i) => { assert.match(row, /^\d{10} 00000 n $/); assert.ok(text.slice(Number(row.slice(0, 10))).startsWith(`${i + 1} 0 obj\n`), `object ${i + 1} is where the table says`); });
  });
}

test("the compressed PDF holds the same page content and is smaller", async () => {
  const svg = figures["all clones stacked"];
  const raw = latin1((await convertSvg(svg, { compress: false })).bytes); const zipped = (await convertSvg(svg)).bytes; const z = latin1(zipped);
  assert.ok(zipped.length < raw.length / 2, `${zipped.length} vs ${raw.length}`); assert.match(z, /\/Filter \/FlateDecode/);
  const m = /4 0 obj\n<< \/Filter \/FlateDecode \/Length (\d+) >>\nstream\n/.exec(z); const start = m.index + m[0].length;
  const content = inflateSync(Buffer.from(zipped.slice(start, start + Number(m[1])))).toString("latin1");
  const plainContent = /4 0 obj\n<< \/Length \d+ >>\nstream\n([\s\S]*?)\nendstream/.exec(raw)[1];
  assert.equal(content, plainContent);
});

test("the output is deterministic", async () => {
  const a = await svgToPdf(figures["clone selection"]); const b = await svgToPdf(figures["clone selection"]);
  assert.deepEqual(Buffer.from(a), Buffer.from(b));
});

test("tooltips (<title>) are not drawn as text", async () => {
  const out = await plain(wrap('<text x="10" y="20" font-size="10">A<title>second base G 40%</title></text>'));
  assert.match(out, /\(A\) Tj/); assert.ok(!/second base/.test(out));
});

test("text: anchors use the Helvetica metrics, bold is wider, entities decode, specials are escaped", async () => {
  const out = await plain(wrap('<text x="100" y="20" font-size="10" text-anchor="middle">Hi</text><text x="100" y="40" font-size="10" text-anchor="end">Hi</text><text x="10" y="60" font-size="10" font-weight="bold">A&amp;B (x) \\</text>'));
  // H 722 + i 222 = 944 per 1000 em: 9.44 at 10 px
  assert.match(out, /1 0 0 -1 95\.28 20 Tm\n\(Hi\) Tj/); assert.match(out, /1 0 0 -1 90\.56 40 Tm\n\(Hi\) Tj/);
  assert.match(out, /\/BaseFont \/Helvetica-Bold /); assert.match(out, /\(A&B \\\(x\\\) \\\\\) Tj/);
});

test("text: minus sign, en dash, ellipsis, non-breaking space map to PDF fonts; look-alikes and unknown characters are counted", async () => {
  const r = await convertSvg(wrap('<text x="1" y="10" font-size="9">\u2212 \u2013 \u2026 a\u00a0b</text><text x="1" y="30" font-size="9">\u2265 \u03bc \u4e2d</text>'), { compress: false });
  const out = latin1(r.bytes);
  assert.match(out, /\(\\201 \\226 \\205 a\\240b\) Tj/, "minus is code 0x81 (glyph minus), en dash 0x96, ellipsis 0x85, nbsp 0xA0");
  assert.match(out, /\/Differences \[129 \/minus\]/); assert.match(out, /\(>= \\265 \?\) Tj/); assert.equal(r.lost, 1, "only the character with no look-alike is lost");
});

test("monospace text uses Courier with 600-unit glyphs", async () => {
  const out = await plain(wrap('<text x="100" y="20" font-size="10" font-family="Courier New, monospace" text-anchor="end">ACGT</text>'));
  assert.match(out, /\/BaseFont \/Courier /); assert.match(out, /1 0 0 -1 76 20 Tm/); // 4 x 600 = 2400 -> 24 px
});

test("shapes: rounded rectangle, circle, dashed stroke, opacity, group transform, evenodd, arcs", async () => {
  const out = await plain(wrap('<g transform="translate(10,20) rotate(30)"><rect x="0" y="0" width="40" height="20" rx="4" fill="#ff0000" stroke="#000" stroke-width="2" stroke-dasharray="3 2" fill-opacity="0.5"/></g>'
    + '<circle cx="50" cy="50" r="5" fill="#00f"/><path d="M0 0 A10 10 0 0 1 20 0 Z" fill="#0f0"/><path d="M0 0h10v10h-10z M2 2h6v6h-6z" fill="#333" fill-rule="evenodd"/><rect width="0" height="5" fill="red"/>'));
  assert.match(out, /0\.866 0\.5 -0\.5 0\.866 10 20 cm/); assert.match(out, /\[3 2\] 0 d/); assert.match(out, /\/G1 gs/); assert.match(out, /\/ca 0\.5 \/CA 1/);
  assert.match(out, /1 0 0 rg/); assert.match(out, /0 0 1 rg/); assert.match(out, /f\*\n/, "even-odd fill"); assert.match(out, / c\n/, "curves for the rounded corners and the arc");
  assert.equal((out.match(/ re\b/g) || []).length, 0, "rectangles are paths, so rounded and square ones draw the same way"); assert.equal((out.match(/1 0 0 rg/g) || []).length, 1, "a rectangle with no width is not drawn");
});

test("tiling pattern fill is written as a PDF pattern with its own matrix", async () => {
  const out = await plain(wrap('<defs><pattern id="h" width="6" height="6" patternUnits="userSpaceOnUse" patternTransform="rotate(45)"><rect width="6" height="6" fill="#fff"/><line x1="0" y1="0" x2="0" y2="6" stroke="#036" stroke-width="1.6"/></pattern></defs><rect x="10" y="10" width="80" height="30" fill="url(#h)" stroke="#036"/>'));
  assert.match(out, /\/PatternType 1/); assert.match(out, /\/Pattern cs \/P1 scn/); assert.match(out, /\/Matrix \[/); assert.match(out, /\/XStep 6 \/YStep 6/);
});

test("unsupported elements are skipped and reported, not drawn wrongly", async () => {
  const r = await convertSvg(wrap('<image href="x.png" x="0" y="0" width="5" height="5"/><text x="1" y="9" font-size="9">kept</text>'), { compress: false });
  assert.deepEqual(r.unsupported, ["image"]); assert.match(latin1(r.bytes), /\(kept\) Tj/);
});

test("parser: entities, comments, CDATA, self-closing tags, single-quoted attributes", () => {
  const root = parseSvg("<?xml version='1.0'?><!-- c --><svg width='5'><g/><text x='1'>a &lt; b &#x41; &#66;<![CDATA[ <raw> ]]></text></svg>");
  const svg = root.kids.find((k) => k.tag === "svg"); assert.equal(svg.attrs.width, "5"); assert.equal(svg.kids.length, 2);
  assert.equal(svg.kids[1].kids.map((k) => k.text).join(""), "a < b A B <raw> ");
});

test("not an SVG is an error, not an empty PDF", async () => { await assert.rejects(() => convertSvg("<html></html>"), /not an SVG/); });
