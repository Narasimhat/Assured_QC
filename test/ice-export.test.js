import test from 'node:test';
import assert from 'node:assert/strict';
import JSZip from 'jszip';
import { planIceExport, buildIceBundle, ICE_HEADERS } from '../src/iceExport.js';
const task = () => ({ name: '=sample <A>', workflow: 'snp', sampleType: 'clone', nuclease: 'SpCas9', controls: ['wt.ab1'], edited: ['read.ab1'], spec: { guides: [{ spacer: 'ACGTACGTACGTACGTACGT' }, { spacer: 'TGCATGCATGCATGCATGCA' }], donors: [{ name: 'D1', sequence: 'A'.repeat(80) }, { name: 'D2', sequence: 'C'.repeat(80) }] } });
const bytes = name => new TextEncoder().encode('ABIF' + name).buffer;
test('two donors give two models per read with all guides and matched controls', () => {
  const t = task(); t.controls.push('wt_rev.ab1'); t.edited.push('rev.ab1');
  const p = planIceExport([t]); assert.equal(p.samples, 1); assert.equal(p.rows.length, 4);
  assert.equal(new Set(p.rows.map(r => r.label)).size, 4);
  assert.equal(p.rows[1].control, 'wt_rev.ab1'); assert.equal(p.rows[1].edited, 'rev.ab1');
  assert.ok(p.rows.every(r => r.alternativeModel && r.guides.includes(',')));
  assert.deepEqual(p.rows.map(r => r.donor), ['A'.repeat(80), 'A'.repeat(80), 'C'.repeat(80), 'C'.repeat(80)]);
});
test('deduplicates identical donors and explicitly excludes incompatible samples', () => {
  const t = task(); t.spec.donors.push(t.spec.donors[0]); assert.equal(planIceExport([t]).rows.length, 2);
  const wrong = { ...t, nuclease: 'Cas12a' }; const inferred = task(); inferred.spec.guides = [{ spacer: '' }];
  assert.equal(planIceExport([wrong, inferred]).excluded.length, 2);
  const ko = { ...task(), workflow: 'knockout', spec: { guides: task().spec.guides, donors: [] } }; assert.equal(planIceExport([ko]).rows[0].donor, '');
});
test('package contains valid nested zips, official workbook headers and exact trace bytes', async () => {
  const out = await buildIceBundle({ tasks: [task()], results: [{ ok: false, warnings: ['poor quality'] }], getBytes: bytes, rowsPerBatch: 1 });
  assert.equal(out.batches, 2); const outer = await JSZip.loadAsync(out.bytes);
  for (const batch of ['batch_001', 'batch_002']) {
    const ab1 = await JSZip.loadAsync(await outer.file(`${batch}/ab1_files.zip`).async('uint8array'));
    assert.equal(Object.keys(ab1.files).length, 2);
    assert.deepEqual(await ab1.file('trace_0001.ab1').async('uint8array'), new Uint8Array(bytes('wt.ab1')));
    const book = await JSZip.loadAsync(await outer.file(`${batch}/ice_definitions.xlsx`).async('uint8array'));
    const sheet = await book.file('xl/worksheets/sheet1.xml').async('string'); ICE_HEADERS.forEach(h => assert.ok(sheet.includes(h)));
    assert.ok(!sheet.includes('<f>')); assert.ok(sheet.includes('ACGTACGTACGTACGTACGT,TGCATGCATGCATGCATGCA'));
  }
  const manifest = JSON.parse(await outer.file('sample_map.json').async('string')); assert.equal(manifest.rows[0].assuredQc.ok, false); assert.equal(manifest.files[0].sha256.length, 64);
});
test('fails closed on missing bytes and oversize pairs', async () => {
  await assert.rejects(buildIceBundle({ tasks: [task()], getBytes: () => undefined }), /invalid AB1/);
  await assert.rejects(buildIceBundle({ tasks: [task()], getBytes: bytes, maxBytes: 10 }), /too large/);
});
