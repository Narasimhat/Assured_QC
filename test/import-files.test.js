import test from 'node:test';
import assert from 'node:assert/strict';
import JSZip from 'jszip';
import { expandImportFiles } from '../src/importFiles.js';
test('reads nested AB1 files from sequencing ZIP and reports unsupported entries', async () => {
 const z = new JSZip(); z.file('run/sample.ab1', 'ABIFsample'); z.file('run/sample.scf', 'ignored');
 const file = new File([await z.generateAsync({type:'uint8array'})], 'run.zip');
 const out = await expandImportFiles([file]);
 assert.equal(out.files.length, 1); assert.equal(out.files[0].name, 'sample.ab1');
 assert.equal(await out.files[0].text(), 'ABIFsample'); assert.equal(out.skipped.length, 1);
});
test('preserves selected files and handles empty selections', async () => {
 const f = new File(['ABIF'], 'a.ab1');
 assert.equal((await expandImportFiles([f])).files[0], f);
 assert.deepEqual(await expandImportFiles([]), {files:[],skipped:[]});
});
