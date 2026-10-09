// The classic form (App.jsx) takes an HTML design report as well as a JSON design file: the guides and the recommended donor fill the form.
import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { readDesignForForm } from "../src/quick.js";

const spacer = "ACGTACGTACGTACGTACGA"; const donor = `${"A".repeat(40)}C${"G".repeat(40)}`;
const letters = (s) => [...s].map((b) => `<span>${b}</span>`).join("");
const report = (title = "Synthetic SNP", withDonor = true) => `<h1>Design: ${title}</h1><h2>2. gRNA Sequences</h2><table><tr><th>Name</th><th>Sequence</th></tr><tr><td>gRNA1</td><td><span>${spacer}</span> <span>TGG</span></td></tr></table><h2>3. Validation Primers</h2><table><tr><td>Primer</td><td>${"T".repeat(20)}</td></tr></table>${withDonor ? `<h2>4. ssODN Donor Templates</h2><h3>ssODN1</h3><div><div>Order this strand</div><div><div>Donor</div><div>${letters(donor)}</div></div></div>` : "<h2>4. Knockout Design</h2><p>No donor required.</p>"}<h2>5. Matched Historical Records</h2>`;

test("classic form: an HTML design report fills the guides and the recommended donor", () => {
  const r = readDesignForForm(report(), { needsDonor: true });
  assert.equal(r.kind, "report"); assert.equal(r.guides, spacer); assert.equal(r.donor, donor);
  assert.match(r.note, /Read 1 guide and 1 donor from the design report "Design: Synthetic SNP"/);
});

test("classic form: a report without a donor says so when the workflow needs one, and stays quiet when it does not", () => {
  const needs = readDesignForForm(report("KO", false), { needsDonor: true });
  assert.equal(needs.donor, ""); assert.match(needs.note, /no ssODN donor/);
  assert.doesNotMatch(readDesignForForm(report("KO", false), { needsDonor: false }).note, /no ssODN donor/);
});

test("classic form: several designs in one report use the first and say so; the junction check points to the .json file", () => {
  const two = readDesignForForm(report("A") + report("B"), {});
  assert.equal(two.guides, spacer); assert.match(two.note, /2 designs; the first was used/);
  assert.match(readDesignForForm(report(), { isJunction: true }).note, /needs the \.json design file/);
});

test("classic form: a JSON design file is passed through unchanged; other files are refused with a reason", () => {
  const spec = { design: { gene: "X" }, guides: [] };
  const r = readDesignForForm(JSON.stringify(spec));
  assert.equal(r.kind, "spec"); assert.deepEqual(r.spec, spec); assert.equal(r.note, "");
  assert.throws(() => readDesignForForm("<p>not a design report</p>"), /No supported/);
  assert.throws(() => readDesignForForm(""), /empty/);
});

test("classic form: the design file box lists .html and .htm in its accepted types", () => {
  const app = readFileSync(new URL("../src/App.jsx", import.meta.url), "utf8");
  const input = app.match(/<input id="design"[^>]*>/)[0];
  assert.match(input, /accept="[^"]*\.html[^"]*\.htm/); assert.match(input, /\.json/);
});
