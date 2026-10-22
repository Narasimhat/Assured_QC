import test from 'node:test';
import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';
import { fields, canonical, sha256 } from '../src/provenance.js';
import { validateGuide, validatePassportInput, mintPassport, verifyEvidence, checkPassport } from '../src/passport.js';
if (!globalThis.crypto) globalThis.crypto = webcrypto;
const metadata = Object.fromEntries(fields.map(([, key]) => [key, key === 'accessUrl' ? 'https://example.org/access' : 'Documented v1']));
const input = { cloneId: 'C04', parentalLine: 'BIHi005-A', targetLocus: 'AAVS1', assembly: 'GRCh38', strategy: 'HDR knock-in',
 spacer: 'GGGGCCACTAGGGACAGGAT', target: '', efficiency: '0', zygosity: 'Not established', tool: 'Assured QC', toolVersion: '0.5.2', parameters: '{"window":30}' };
const file = (name, text) => { const bytes = new TextEncoder().encode(text); return { name, size: bytes.length, arrayBuffer: async () => bytes.buffer }; };
const attachments = ['design', 'raw-sequencing', 'analysis'].map(role => ({ role, file: file(role + '.json', role) }));
test('SpCas9 guide normalization, invalid characters, length and supplied PAM context', () => {
 assert.equal(validateGuide(' '+input.spacer.toLowerCase()+' ').spacer5to3, input.spacer);
 assert.match(validateGuide(input.spacer).pamCheck, /not assessed/);
 assert.match(validateGuide(input.spacer, input.spacer+'AGG').pamCheck, /supplied target/);
 for (const bad of ['ACGT', input.spacer.replace('G','N'), input.spacer+'U']) assert.throws(() => validateGuide(bad));
 assert.throws(() => validateGuide(input.spacer, input.spacer+'AGA'));
 assert.throws(() => validateGuide(input.spacer, 'A'.repeat(20)+'AGG'));
});
test('required lineage, analysis provenance, strategy and numerical boundaries', () => {
 assert.equal(validatePassportInput(input).editingAnalysis.reportedFractionPercent, 0);
 assert.equal(validatePassportInput({...input,efficiency:''}).editingAnalysis.reportedFractionPercent, null);
 for (const patch of [{cloneId:''},{parentalLine:''},{assembly:''},{efficiency:'NaN'},{efficiency:101},{efficiency:-1},{parameters:'[]'},{parameters:'null'},{parameters:'{'},{toolVersion:''},{strategy:'other'},{zygosity:'confirmed'}]) assert.throws(() => validatePassportInput({...input,...patch}));
});
test('multi-file passport survives JSON roundtrip, checks each evidence entry and leaves QC unassessed', async () => {
 const cert = JSON.parse(JSON.stringify(await mintPassport(input,metadata,attachments,'test')));
 assert.equal((await checkPassport(cert)).cloneId,'C04');
 assert.equal((await verifyEvidence(await attachments[0].file.arrayBuffer(),cert,0)).matches,true);
 assert.equal((await verifyEvidence(new TextEncoder().encode('Design'),cert,0)).matches,false);
 assert.match(cert.record.qc.karyotype.status,/not supplied/);
 const withQc=await mintPassport(input,metadata,[...attachments,{role:'karyotype',file:file('qc.pdf','PDF')}],'test');
 assert.equal(withQc.record.qc.karyotype.status,'evidence attached; not assessed');
});
test('missing evidence, invalid roles, empty files and file/count limits rejected', async () => {
 for (const values of [[],attachments.slice(1),[...attachments,{role:'invalid',file:file('a','a')}],[...attachments,{role:'identity',file:file('empty','')}],Array(21).fill(attachments[0]),[...attachments,{role:'identity',file:{size:51*1024*1024}}]]) await assert.rejects(mintPassport(input,metadata,values,'test'));
 await assert.rejects(mintPassport(input,{},attachments,'test'));
});
test('passport metadata tampering, malformed manifests and invalid selections rejected', async () => {
 const cert=await mintPassport(input,metadata,attachments,'test');
 cert.record.parentalLine='different';
 await assert.rejects(checkPassport(cert),/integrity/);
 cert.recordSha256=await sha256(new TextEncoder().encode(canonical(cert.record)));
 assert.equal((await checkPassport(cert)).parentalLine,'different'); // Self-issued is not authentication.
 await assert.rejects(verifyEvidence(new Uint8Array(),cert,99));
 cert.record.evidence=[];
 cert.recordSha256=await sha256(new TextEncoder().encode(canonical(cert.record)));
 await assert.rejects(checkPassport(cert),/manifest/);
 for (const bad of [null,{},[]]) await assert.rejects(checkPassport(bad));
});
