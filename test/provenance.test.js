import test from 'node:test';
import assert from 'node:assert/strict';
import { webcrypto } from 'node:crypto';
import { fields, screen, mint, verify, canonical, sha256 } from '../src/provenance.js';
if (!globalThis.crypto) globalThis.crypto = webcrypto;
const bytes = new TextEncoder().encode('QC,97.2');
const metadata = Object.fromEntries(fields.map(([, key]) => [key, key === 'accessUrl' ? 'https://example.org/access' : 'Documented v1']));
test('metadata screening rejects incomplete and non-HTTPS records', async () => {
  assert.equal(screen({}).score, 0);
  assert.equal(screen(metadata).score, 100);
  await assert.rejects(mint(bytes, 'qc.csv', { ...metadata, accessUrl: 'http://example.org' }, 'test'));
});
test('certificate roundtrip and one-character file change', async () => {
  const cert = await mint(bytes, 'qc.csv', metadata, 'test');
  assert.equal(await verify(bytes, JSON.parse(JSON.stringify(cert))), true);
  assert.equal(await verify(new TextEncoder().encode('QC,97.3'), cert), false);
});
test('altered metadata and malformed certificates rejected', async () => {
  const cert = await mint(bytes, 'qc.csv', metadata, 'test');
  cert.record.metadata.creator = 'Other lab';
  await assert.rejects(verify(bytes, cert), /integrity/);
  for (const bad of [null, {}, [], { record: {} }]) await assert.rejects(verify(bytes, bad));
});
test('canonical ordering and explicit self-issued trust boundary', async () => {
  assert.equal(canonical({ b: 2, a: 'µ' }), canonical({ a: 'µ', b: 2 }));
  const cert = await mint(bytes, 'qc.csv', metadata, 'test');
  cert.record.metadata.creator = 'Other lab';
  cert.recordSha256 = await sha256(new TextEncoder().encode(canonical(cert.record)));
  assert.equal(await verify(bytes, cert), true); // Integrity is not issuer authentication.
});
