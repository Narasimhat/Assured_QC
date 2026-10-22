import { canonical, sha256, screen } from './provenance.js';

export const passportVersion = 'assured-ipsc-passport/1';
export const evidenceRoles = ['design', 'raw-sequencing', 'analysis', 'identity', 'pluripotency', 'karyotype', 'mycoplasma', 'off-target'];
const requiredRoles = ['design', 'raw-sequencing', 'analysis'];
export const strategies = ['HDR knock-in', 'Indel knockout', 'Base editing', 'Prime editing'];
export const zygosities = ['Not established', 'Homozygous (reported)', 'Heterozygous (reported)', 'Compound heterozygous (reported)', 'Bulk pool'];
export function validateGuide(spacer, target = '') {
  const sequence = String(spacer || '').trim().toUpperCase();
  if (!/^[ACGT]{20}$/.test(sequence)) throw new Error('SpCas9 spacer must contain exactly 20 A/C/G/T bases, without PAM.');
  const context = String(target || '').trim().toUpperCase();
  if (context && (!/^[ACGT]{21}GG$/.test(context) || context.slice(0, 20) !== sequence)) {
    throw new Error('Optional target must be the matching 20-base spacer followed by an NGG PAM (23 bases, target-strand 5′→3′).');
  }
  return { spacer5to3: sequence, nuclease: 'SpCas9', target5to3: context || null,
    pamCheck: context ? 'matches supplied target; genomic mapping not checked' : 'not assessed; no target context supplied' };
}
function required(input, key) {
  if (typeof input[key] !== 'string' || !input[key].trim()) throw new Error(`Missing ${key}.`);
  return input[key].trim();
}
export function validatePassportInput(input) {
  const cloneId = required(input, 'cloneId');
  const parentalLine = required(input, 'parentalLine');
  const targetLocus = required(input, 'targetLocus');
  const assembly = required(input, 'assembly');
  if (!strategies.includes(input.strategy) || !zygosities.includes(input.zygosity)) throw new Error('Choose a supported strategy and zygosity.');
  const guide = validateGuide(input.spacer, input.target);
  const efficiency = input.efficiency === '' || input.efficiency == null ? null : Number(input.efficiency);
  if (efficiency !== null && (!Number.isFinite(efficiency) || efficiency < 0 || efficiency > 100)) throw new Error('Reported editing fraction must be 0–100, or blank if not established.');
  const tool = required(input, 'tool');
  const toolVersion = required(input, 'toolVersion');
  let parameters;
  try { parameters = JSON.parse(input.parameters); } catch { throw new Error('Analysis parameters must be a JSON object.'); }
  if (!parameters || Array.isArray(parameters) || typeof parameters !== 'object') throw new Error('Analysis parameters must be a JSON object.');
  return { cloneId, parentalLine, design: { targetLocus, assembly, strategy: input.strategy, guide },
    editingAnalysis: { reportedFractionPercent: efficiency, zygosity: input.zygosity,
      tool, toolVersion, parameters, assessment: 'user-reported; not independently validated' } };
}
const recordHash = record => sha256(new TextEncoder().encode(canonical(record)));
export async function mintPassport(input, metadata, attachments, source) {
  const core = validatePassportInput(input);
  const clean = Object.fromEntries(Object.keys(metadata).filter(k => typeof metadata[k] === 'string').map(k => [k, metadata[k].trim()]));
  const assessment = screen(clean);
  if (assessment.score !== 100) throw new Error('Complete the twelve FAIR metadata fields above first.');
  if (!Array.isArray(attachments) || attachments.length > 20) throw new Error('Use at most 20 evidence files.');
  for (const role of requiredRoles) if (!attachments.some(a => a.role === role)) throw new Error(`Attach ${role} evidence.`);
  let total = 0;
  const evidence = [];
  for (const attachment of attachments) {
    const { role, file } = attachment;
    if (!evidenceRoles.includes(role) || !file || file.size > 50 * 1024 * 1024 || file.size === 0) throw new Error('Each evidence file must be non-empty, at most 50 MB, with a supported role.');
    total += file.size;
    if (total > 200 * 1024 * 1024) throw new Error('Evidence total exceeds 200 MB.');
    evidence.push({ role, fileName: file.name, fileSize: file.size, sha256: await sha256(await file.arrayBuffer()) });
  }
  const qc = Object.fromEntries(evidenceRoles.filter(r => !requiredRoles.includes(r)).map(role => [role, {
    status: evidence.some(e => e.role === role) ? 'evidence attached; not assessed' : 'not supplied; not assessed',
  }]));
  const record = { schemaVersion: passportVersion, source, timestampUtc: new Date().toISOString(), ...core,
    metadata: clean, metadataAssessment: assessment, evidence, qc, anchorStatus: 'unanchored',
    limitations: 'Self-issued record. File equality does not establish cell identity, biological quality, or issuer authenticity.' };
  return { record, recordSha256: await recordHash(record) };
}
export async function checkPassport(certificate) {
  const record = certificate?.record;
  if (!record || record.schemaVersion !== passportVersion || typeof certificate.recordSha256 !== 'string' || !/^[a-f0-9]{64}$/.test(certificate.recordSha256)) throw new Error('Invalid or unsupported iPSC passport.');
  if (await recordHash(record) !== certificate.recordSha256) throw new Error('Passport integrity check failed: record changed.');
  if (typeof record.cloneId !== 'string' || !record.cloneId.trim() || typeof record.parentalLine !== 'string' || typeof record.timestampUtc !== 'string') throw new Error('Invalid passport identity fields.');
  if (!Array.isArray(record.evidence) || record.evidence.length < 3 || record.evidence.length > 20 || record.evidence.some(e => !e || !evidenceRoles.includes(e.role) || typeof e.sha256 !== 'string' || !/^[a-f0-9]{64}$/.test(e.sha256) || typeof e.fileName !== 'string' || !Number.isSafeInteger(e.fileSize) || e.fileSize < 1)) throw new Error('Invalid evidence manifest.');
  for (const role of requiredRoles) if (!record.evidence.some(e => e.role === role)) throw new Error(`Missing ${role} evidence manifest.`);
  return record;
}
export async function verifyEvidence(bytes, certificate, index) {
  const record = await checkPassport(certificate);
  if (!Number.isInteger(index) || !record.evidence[index]) throw new Error('Select a valid evidence entry.');
  const expected = record.evidence[index];
  const actualHash = await sha256(bytes);
  return { matches: actualHash === expected.sha256 && bytes.byteLength === expected.fileSize,
    expectedHash: expected.sha256, actualHash };
}
