"""Stateless research-only wrapper around public ICE 1.2.0."""
import hashlib
import json
import math
import re
import sys
import tempfile
from pathlib import Path
from flask import Flask, request, jsonify
from werkzeug.exceptions import HTTPException

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'vendor'))
COMMIT = '310495fe536b7d2eafcee80f145c55f7c65d947c'
app = Flask(__name__)
app.config['MAX_CONTENT_LENGTH'] = 4_000_000

@app.after_request
def private_response(response):
    response.headers['Cache-Control'] = 'no-store'
    return response

@app.errorhandler(413)
def too_large(error):
    return jsonify(error='ICE request exceeds 4 MB. Submit one control/sample pair.'), 413

def validate_settings(settings):
    if not isinstance(settings, dict):
        raise ValueError('Settings must be an object.')
    if settings.get('nuclease', 'SpCas9') != 'SpCas9':
        raise ValueError('Public ICE integration supports SpCas9 only.')
    guides = settings.get('guides')
    if not isinstance(guides, list) or not 1 <= len(guides) <= 3 or any(not isinstance(g, str) or not re.fullmatch('[ACGT]{17,23}', g) for g in guides):
        raise ValueError('Supply one to three DNA guide sequences, 17–23 bases, without PAM.')
    donor = settings.get('donor') or None
    if donor is not None and (not isinstance(donor, str) or not re.fullmatch('[ACGT]{20,300}', donor)):
        raise ValueError('ICE takes one donor per comparison, 20–300 DNA bases.')
    return guides, donor

def analyze(control, sample, settings):
    guides, donor = validate_settings(settings)
    for data in [control, sample]:
        if not data.startswith(b'ABIF') or len(data) > 1_950_000:
            raise ValueError('Use AB1 files below 1.95 MB each.')
    # Import the core directly: avoid the upstream CLI's pandas dependency.
    from ice.classes.sanger_analysis import SangerAnalysis
    with tempfile.TemporaryDirectory(prefix='assured-ice-') as folder:
        folder = Path(folder)
        (folder/'control.ab1').write_bytes(control)
        (folder/'sample.ab1').write_bytes(sample)
        engine = SangerAnalysis(verbose=False)
        engine.debug = False
        engine.allprops = False
        engine.initialize_with(control_path=str(folder/'control.ab1'), edited_path=str(folder/'sample.ab1'),
                               gRNA_sequences=','.join(guides), indel_max_size=20,
                               base_outputname=str(folder/'result'), donor=donor, allprops=False)
        engine.analyze_sample()
        raw = engine.results.to_json(engine.guide_targets, engine.warnings)
    if not all(isinstance(raw.get(k), (float, int)) and math.isfinite(raw[k]) for k in ('ice', 'rsq', 'ko_score')):
        raise ValueError('ICE did not produce a usable fit.')
    for key, value in list(raw.items()):
        if isinstance(value, float) and not math.isfinite(value):
            raw[key] = None
            raw['notes'] += f' Non-finite {key} diagnostic is unavailable.'
    if raw['rsq'] < 0.8:
        raw['notes'] += ' Low ICE fit (R-squared <0.8); review before comparing scores.'
    return {'ok': True, 'engine': 'Synthego public ICE', 'version': '1.2.0', 'commit': COMMIT, 'compatibility_patch': 'Biopython Clustal output formatting only',
            'control_sha256': hashlib.sha256(control).hexdigest(), 'sample_sha256': hashlib.sha256(sample).hexdigest(),
            'settings': {'guides': guides, 'donor': donor, 'nuclease': 'SpCas9'}, 'result': raw}

@app.route('/api/ice', methods=['GET', 'POST'])
def endpoint():
    if request.method == 'GET':
        return jsonify(engine='Synthego public ICE', version='1.2.0', commit=COMMIT, use='research/nonprofit')
    try:
        settings = json.loads(request.form.get('settings', '{}'))
        validate_settings(settings)
        control, sample = request.files.get('control'), request.files.get('sample')
        if control is None or sample is None:
            raise ValueError('A control and a sample AB1 are required.')
        result = analyze(control.read(1_950_001), sample.read(1_950_001), settings)
        # Reject non-finite values rather than silently turning failures into zeros.
        return app.response_class(json.dumps(result, allow_nan=False), mimetype='application/json')
    except HTTPException:
        raise
    except (ValueError, TypeError, KeyError) as error:
        return jsonify(ok=False, error=str(error)[:300]), 422
    except Exception:
        return jsonify(ok=False, error='ICE could not fit this trace pair. Check read quality and guide coverage.'), 422

if __name__ == '__main__':
    app.run(host='127.0.0.1', port=8786)
