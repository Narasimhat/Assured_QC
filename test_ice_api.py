import io,json,unittest
from pathlib import Path
from api.ice import app,validate_settings

class ApiTests(unittest.TestCase):
    def test_public_upstream_native_example(self):
        fixtures=Path(__file__).parent/'test'/'ice-fixtures'
        response=app.test_client().post('/api/ice',data={
            'settings':json.dumps({'guides':['AACCAGTTGCAGGCGCCCCA']}),
            'control':(io.BytesIO((fixtures/'good_example_control.ab1').read_bytes()),'control.ab1'),
            'sample':(io.BytesIO((fixtures/'good_example_edited.ab1').read_bytes()),'sample.ab1')})
        self.assertEqual(response.status_code,200,response.json)
        self.assertEqual(response.json['result']['ice'],77)
        self.assertEqual(response.json['result']['ko_score'],76)
        self.assertEqual(response.json['result']['rsq'],0.98)
        self.assertEqual(len(response.json['control_sha256']),64)

    def test_rejects_unsupported_inputs(self):
        for s in [{},{'guides':['../file']},{'guides':['A'*20],'nuclease':'Cas12a'},{'guides':['A'*20],'donor':['A'*30]}]:
            with self.assertRaises(ValueError): validate_settings(s)
    def test_invalid_upload_and_missing_pair(self):
        client=app.test_client()
        self.assertEqual(client.post('/api/ice',data={'settings':'[]'}).status_code,422)
        r=client.post('/api/ice',data={'settings':json.dumps({'guides':['A'*20]}),'control':(io.BytesIO(b'not abif'),'wt.ab1'),'sample':(io.BytesIO(b'not abif'),'sample.ab1')})
        self.assertEqual(r.status_code,422)
        self.assertIn('AB1',r.json['error'])
        self.assertEqual(r.headers['Cache-Control'],'no-store')
    def test_size_limit(self):
        self.assertEqual(app.test_client().post('/api/ice',data=b'x'*4_000_001,content_type='application/octet-stream').status_code,413)

if __name__=='__main__':unittest.main()
