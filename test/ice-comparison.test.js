import test from 'node:test';
import assert from 'node:assert/strict';
import {iceCases,iceComparisonHtml} from '../src/iceComparison.js';
const task=()=>({controls:['wt'],edited:['clone'],workflow:'snp',spec:{guides:[{spacer:'ACGTACGTACGTACGTACGT'}],donors:[{name:'one',sequence:'A'.repeat(40)},{name:'two',sequence:'C'.repeat(40)}]}});
test('ICE donor hypotheses are separate and inferred or combined reads are rejected',()=>{
 const t=task();assert.equal(iceCases(t).length,2);assert.equal(iceCases(t)[0].settings.donor,'A'.repeat(40));
 t.edited.push('reverse');assert.throws(()=>iceCases(t),/one control/);t.edited.pop();t.spec.guides[0].spacer='';assert.throws(()=>iceCases(t),/explicit guide/);
});
test('ICE failures stay unavailable and imported strings are escaped',()=>{
 const html=iceComparisonHtml([{name:'<script>',summary:{indelPct:100,koScorePct:100},ice:[{ok:false,label:'<img>',error:'failed'}]}]);
 assert.match(html,/Unavailable/);assert.doesNotMatch(html,/<script>|<img>/);assert.match(html,/&lt;script&gt;/);
});
