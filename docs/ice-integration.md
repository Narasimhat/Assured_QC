# Optional ICE comparison

For research and nonprofit use only. Copyright and full licence are retained in
`vendor/ice/LICENSE` and linked in the app. Public Synthego ICE version 1.2.0 is
pinned at source commit 310495fe536b7d2eafcee80f145c55f7c65d947c. This is not a claim
of parity with the current hosted ICE service.

After a Quick start analysis, choose **Run ICE comparison on server**. This is
opt-in: it sends each original control/sample AB1 pair and guides/donor to the
same-origin Python endpoint, sequentially. Assured QC itself remains local.
Temporary files are removed when the request finishes; the application does not
save uploads. Responses are non-cacheable. Hosting infrastructure still applies
its normal request logging. One pair per request, 1.95 MB per file, 4 MB total.

SpCas9 with 1–3 explicit guides and single-read pairs is supported. Up to two
donors are evaluated as **separate hypotheses**, never merged or maximized.
Unsupported workflows and failures remain visible as unavailable. ICE does not
change Assured QC calls, ranks or coding predictions. KO screening scores are not
probabilities of protein loss. Low fit (R² <0.8) is explicitly flagged.

The HTML report includes the comparison table. Download ICE results + provenance
to retain input SHA-256 values, source revision, settings, raw scores and warnings.

The vendored core is unchanged except `PairAlignment.align_list_to_clustal` uses
`format(alignment, "clustal")` in place of the removed Biopython `.format` method.
This affects diagnostic output only. UPSTREAM.json records original source hashes.
The wrapper represents non-finite diagnostic values as null with an explicit note;
missing/non-finite primary scores make a fit fail. Python dependencies are pinned.

Local development: install requirements in an isolated Python environment; run
`python api/ice.py`, then `npm run dev`. Vite proxies `/api/ice` to localhost:8786.
On Vercel the Python function and Vite frontend share the same origin.
