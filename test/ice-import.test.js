import test from 'node:test';
import assert from 'node:assert/strict';
import { readIceCsv, readIceMap, suggestIceColumns, compareIce, comparisonCsv } from '../src/iceImport.js';
const map = { schema: 'assured-qc-ice-export/1', rows: [1, 2].map(i => ({ label: `sample_D${i}`, sample: 'sample', read: 1, donorName: `Donor ${i}`, assuredQc: { ok: true, summary: { intendedEditPct: 30, indelPct: 10 }, quality: { r2: 0.96 } } })) };
test('separate donor results match exact export labels and preserve zero/missing values', () => {
  const table = readIceCsv('Label,HDR,ICE,KO,R2,Status\nsample_D1,0,42,10,0.98,OK\nsample_D2,25,,N/A,0.91,"Low fit, review"');
  const out = compareIce(table, readIceMap(JSON.stringify(map)), suggestIceColumns(table.headers));
  assert.equal(out.rows[0].ice.hdr, 0); assert.equal(out.rows[1].ice.hdr, 25);
  assert.equal(out.rows[1].ice.editing, null); assert.equal(out.rows[1].ice.ko, null);
  assert.equal(out.rows[1].source.donorName, 'Donor 2'); assert.equal(out.missing.length, 0);
  assert.match(comparisonCsv(out), /Low fit, review/);
});
test('duplicate and unknown labels never receive a guessed comparison', () => {
  const t = readIceCsv('Label,HDR\nsample_D1,30\nsample_D1,40\nsample,50');
  const out = compareIce(t, map, suggestIceColumns(t.headers));
  assert.ok(out.rows.every(r => r.source === null)); assert.equal(out.missing.length, 2);
});
test('units are explicit, invalid scores flagged and missing values not zero', () => {
  const t = readIceCsv('Label,HDR,ICE,KO,R2\nsample_D1,0.3,50%,1,98\nsample_D2,2,error,,101');
  const out = compareIce(t, map, suggestIceColumns(t.headers), 'fraction', 'percent');
  assert.deepEqual(out.rows[0].ice, { hdr: 30, editing: 50, ko: 100, r2: .98, status: '' });
  assert.equal(out.rows[1].issues.length, 3); assert.equal(out.rows[1].ice.hdr, null);
  assert.throws(() => compareIce(t, map, { label: 0, hdr: 0 }), /different column/);
});
test('invalid maps and ragged input fail; CSV output neutralizes spreadsheet formulas', () => {
  assert.throws(() => readIceMap('{"rows":[]}'), /original ICE/);
  assert.throws(() => readIceMap(JSON.stringify({ ...map, rows: [map.rows[0], map.rows[0]] })), /duplicate/);
  assert.throws(() => readIceCsv('Label,HDR\na,1,2'), /Inconsistent/);
  const t = readIceCsv('Label,HDR,Status\nsample_D1,30,=1+1');
  assert.match(comparisonCsv(compareIce(t, map, suggestIceColumns(t.headers))), /'\=1\+1/);
});
