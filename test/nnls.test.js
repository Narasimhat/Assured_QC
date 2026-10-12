import test from "node:test";
import assert from "node:assert/strict";
import { nnlsGram } from "../src/nnls.js";

function gram(columns, b) {
  const n = columns.length;
  const AtA = Array.from({ length: n }, (_, i) => Float64Array.from({ length: n }, (_, k) => columns[i].reduce((s, v, r) => s + v * columns[k][r], 0)));
  const Atb = Float64Array.from(columns, (c) => c.reduce((s, v, r) => s + v * b[r], 0));
  return { AtA, Atb };
}

test("recovers exact non-negative weights", () => {
  const columns = [[1, 0, 0, 1], [0, 1, 0, 1], [0, 0, 1, 0], [1, 1, 1, 0]];
  const truth = [0.5, 0.25, 0, 0.25];
  const b = [0, 1, 2, 3].map((r) => columns.reduce((s, c, i) => s + truth[i] * c[r], 0));
  const { AtA, Atb } = gram(columns, b);
  const w = nnlsGram(AtA, Atb);
  truth.forEach((value, i) => assert.ok(Math.abs(w[i] - value) < 1e-6, `w[${i}]=${w[i]} vs ${value}`));
});

test("never returns a negative weight and matches the known solution of a constrained problem", () => {
  // minimise |w1*[1,0] + w2*[0,1] - [-1, 2]|: unconstrained w=(-1,2); with w>=0 the answer is (0,2)
  const columns = [[1, 0], [0, 1]];
  const { AtA, Atb } = gram(columns, [-1, 2]);
  const w = nnlsGram(AtA, Atb);
  assert.ok(Math.abs(w[0]) < 1e-9 && Math.abs(w[1] - 2) < 1e-9);
});

test("agrees with the brute-force optimum on random small problems", () => {
  let seed = 12345;
  const rnd = () => { seed = (seed * 1664525 + 1013904223) % 4294967296; return seed / 4294967296; };
  for (let trial = 0; trial < 40; trial += 1) {
    const rows = 6; const n = 3;
    const columns = Array.from({ length: n }, () => Array.from({ length: rows }, () => rnd()));
    const b = Array.from({ length: rows }, () => rnd() * 2 - 0.3);
    const { AtA, Atb } = gram(columns, b);
    const w = nnlsGram(AtA, Atb);
    const cost = (x) => b.reduce((s, v, r) => s + (v - columns.reduce((t, c, i) => t + x[i] * c[r], 0)) ** 2, 0);
    let best = Infinity;
    for (let a = 0; a <= 200; a += 1) for (let c = 0; c <= 200; c += 1) for (let d = 0; d <= 200; d += 1) {
      const x = [a * 0.01, c * 0.01, d * 0.01];
      if (a % 4 || c % 4 || d % 4) continue;
      best = Math.min(best, cost(x));
    }
    assert.ok(cost(w) <= best + 1e-9, `trial ${trial}: ${cost(w)} > grid ${best}`);
    assert.ok([...w].every((v) => v >= 0));
  }
});
