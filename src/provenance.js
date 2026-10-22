export const fields = [
  ['F', 'identifier', 'Dataset / batch identifier'], ['F', 'title', 'Title'], ['F', 'keywords', 'Search keywords'],
  ['A', 'accessUrl', 'HTTPS repository or access-policy URL'], ['A', 'accessPolicy', 'Access conditions'], ['A', 'contact', 'Data steward contact'],
  ['I', 'schema', 'Schema name and version'], ['I', 'vocabulary', 'Controlled vocabulary'], ['I', 'units', 'Units or not-applicable rationale'],
  ['R', 'license', 'License or documented usage terms'], ['R', 'protocol', 'SOP / protocol version'], ['R', 'creator', 'Responsible organization'],
];
export const version = 'assured-provenance/1';
export function canonical(value) {
  if (value === null || typeof value !== 'object') return JSON.stringify(value);
  if (Array.isArray(value)) return '[' + value.map(canonical).join(',') + ']';
  return '{' + Object.keys(value).sort().map(k => JSON.stringify(k) + ':' + canonical(value[k])).join(',') + '}';
}
export async function sha256(bytes) {
  if (!globalThis.crypto?.subtle) throw new Error('SHA-256 requires HTTPS or localhost.');
  return Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256', bytes)), b => b.toString(16).padStart(2, '0')).join('');
}
const hashRecord = record => sha256(new TextEncoder().encode(canonical(record)));
export function screen(metadata) {
  const checks = fields.map(([pillar, key, label]) => {
    let passed = typeof metadata[key] === 'string' && metadata[key].trim().length > 0;
    if (key === 'accessUrl' && passed) {
      try { const url = new URL(metadata[key]); passed = url.protocol === 'https:' && !!url.hostname; } catch { passed = false; }
    }
    return { pillar, key, label, passed };
  });
  return { profile: version, score: Math.round(100 * checks.filter(c => c.passed).length / checks.length), checks };
}
export async function mint(bytes, name, metadata, source) {
  const clean = Object.fromEntries(fields.map(([, key]) => [key, String(metadata[key] || '').trim()]));
  const assessment = screen(clean);
  if (assessment.score !== 100) throw new Error('Complete all metadata checks before creating a record.');
  const record = { version, source, timestampUtc: new Date().toISOString(), fileName: name, fileSize: bytes.byteLength,
    fileSha256: await sha256(bytes), metadata: clean, assessment, anchorStatus: 'unanchored' };
  return { record, recordSha256: await hashRecord(record) };
}
export async function verify(bytes, certificate) {
  const record = certificate?.record;
  if (!record || record.version !== version || !/^[a-f0-9]{64}$/.test(certificate.recordSha256 || '') || !/^[a-f0-9]{64}$/.test(record.fileSha256 || '')) throw new Error('Invalid or unsupported certificate.');
  if (await hashRecord(record) !== certificate.recordSha256) throw new Error('Certificate integrity check failed: record changed.');
  return await sha256(bytes) === record.fileSha256;
}
