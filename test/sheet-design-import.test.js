import test from "node:test";
import assert from "node:assert/strict";
import { parseDesignReport } from "../src/designReport.js";
import { parseCsv, parseSampleSheet, resolveSheet, autoPair, looksLikeControl, directionOf } from "../src/sheet.js";

const spacer = "ACGTACGTACGTACGTACGA"; const donor = `${"A".repeat(40)}C${"G".repeat(40)}`; const opposite = `${"C".repeat(40)}G${"T".repeat(40)}`;
const letters = (s) => [...s].map((b) => `<span>${b}</span>`).join("");
function fixture({ two = false } = {}) {
  return `<h1>Design: Synthetic SNP</h1><h2>2. gRNA Sequences</h2><table><tr><th>Name</th><th>Sequence</th></tr><tr><td>gRNA1</td><td><span>${spacer}</span> <span>TGG</span></td></tr></table><h2>3. Validation Primers</h2><table><tr><td>Primer</td><td>${"T".repeat(20)}</td></tr></table><h2>4. ssODN Donor Templates</h2><h3>ssODN1</h3><div><div>Order this strand</div><div><div>WT</div><div>${letters("A".repeat(81))}</div></div><div><div>Donor</div><div>${letters(donor)}</div></div></div><div><div>Reference strand</div><div><div>Donor</div><div>${letters(opposite)}</div></div></div>${two ? `<h3>ssODN2</h3><div><div>Donor</div><div>${"T".repeat(40)}A${"G".repeat(40)}</div></div>` : ""}<h2>5. Matched Historical Records</h2><table><tr><td>Historical gRNA</td><td>${"C".repeat(20)}</td></tr></table>`;
}

test("design report: guide, annotated donor, omitting PAM, wild type, primers and history; the opposite strand is not a second donor", () => {
  const [d] = parseDesignReport(fixture()); assert.equal(d.guides.length, 1); assert.equal(d.guides[0].sequence, spacer);
  assert.equal(d.donors.length, 1); assert.equal(d.donors[0].sequence, donor); assert.equal(d.donors[0].recommended, true);
});
test("design report: distinct donors are kept and designs in one file are separated", () => {
  const ds = parseDesignReport(fixture({ two: true }) + fixture()); assert.equal(ds.length, 2); assert.equal(ds[0].donors.length, 2); assert.equal(ds[1].donors.length, 1);
});
test("design report: a knockout report has guides and no donor; an unrelated page is refused; code in the report is dropped, not read", () => {
  const ko = `${fixture().split("<h2>4.")[0]}<h2>4. Knockout Design</h2><p>No donor required.</p>`; assert.equal(parseDesignReport(ko)[0].donors.length, 0);
  assert.throws(() => parseDesignReport("<p>not a design</p>"), /No supported/);
  const hostile = `${fixture()}<script>const g = "ACGTACGTACGTACGTACGTAA";</script><img src="https://invalid.test/x" onerror="boom()">`; const [d] = parseDesignReport(hostile); assert.equal(d.guides.length, 1);
});

test("csv: quotes, tabs and a byte-order mark", () => {
  assert.deepEqual(parseCsv('\uFEFFa,b\n"x, y",2\n'), [["a", "b"], ["x, y", "2"]]); assert.deepEqual(parseCsv("a\tb\n1\t2"), [["a", "b"], ["1", "2"]]);
});
test("sample sheet: aliases, guides, donor, pools; files are matched without extension or case", () => {
  const sheet = parseSampleSheet("Sample,File,WT,gRNA,ssODN,Type\nclone1,Clone1.ab1,ctrl,\"ACGTACGTACGTACGTACGT, TTTTACGTACGTACGTACGT\",,clone\npool1,pool,ctrl,ACGTACGTACGTACGTACGT,AC GT,Pool\n");
  assert.equal(sheet.errors.length, 0); assert.deepEqual(sheet.rows[0].guides.length, 2); assert.equal(sheet.rows[1].type, "pool"); assert.equal(sheet.rows[1].donor, "ACGT");
  const r = resolveSheet(sheet.rows, ["CLONE1.AB1", "Ctrl.ab1", "Pool.ab1"]); assert.equal(r.problems.length, 0); assert.deepEqual(r.jobs[0].edited, ["CLONE1.AB1"]); assert.deepEqual(r.jobs[0].control, ["Ctrl.ab1"]);
  const missing = resolveSheet(sheet.rows, ["Clone1.ab1"]); assert.ok(missing.problems.some((p) => /control "ctrl" was not found/.test(p)));
});

test("auto pairing: run code groups, a control-looking name, forward/reverse, pools", () => {
  assert.ok(looksLikeControl("GENE_A_WT.744.ab1") && looksLikeControl("XB17-C.674.ab1") && looksLikeControl("GENE-UTF.766.ab1") && !looksLikeControl("GENE_B_2.1.ab1"));
  assert.deepEqual(directionOf("clone3_R.ab1"), { stem: "clone3", direction: "R" });
  const files = ["X_WT.740.ab1", "X_1.740.ab1", "X_2.740.ab1", "Y_ctrl_F.810.ab1", "Y_ctrl_R.810.ab1", "Y_4_F.810.ab1", "Y_4_R.810.ab1", "orphan.999.ab1"];
  const r = autoPair(files); const byName = Object.fromEntries(r.pairs.map((p) => [p.sample, p]));
  assert.deepEqual(byName["X_1.740"].control, ["X_WT.740.ab1"]); assert.equal(byName["X_1.740"].confidence, "high");
  const y = r.pairs.find((p) => p.sample.startsWith("Y_4")); assert.equal(y.edited.length, 2); assert.equal(y.control.length, 2);
  assert.ok(r.unpaired.some((u) => u.name === "orphan.999.ab1"));
});
test("auto pairing: two controls that tie are paired with a low confidence and reported, never silently", () => {
  const r = autoPair(["S_WT.5.ab1", "S_wt.5.ab1", "S_1.5.ab1"]); assert.equal(r.pairs[0].confidence, "low"); assert.ok(r.problems.length >= 1);
});
test("auto pairing: the parental line name is the control when nothing else looks like one", () => {
  const r = autoPair(["LINE-A.812.ab1", "LINE-A-50_KO_Cl2.ab1"]); assert.deepEqual(r.pairs[0].control, ["LINE-A.812.ab1"]); assert.equal(r.pairs[0].confidence, "low");
});
