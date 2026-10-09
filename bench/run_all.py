"""Run every available tool on one job file and write results/<tool>.json in the common shape.
Usage: python bench/run_all.py jobs.json outdir [--tools assured_qc,ice,peaksplit,traceedit] [--config bench/config.json]
"""
import argparse, json, os, sys, time
sys.path.insert(0, os.path.dirname(__file__))
from adapters import ADAPTERS

ap = argparse.ArgumentParser(); ap.add_argument("jobs"); ap.add_argument("outdir"); ap.add_argument("--tools", default=",".join(ADAPTERS)); ap.add_argument("--config", default=os.path.join(os.path.dirname(__file__), "config.json"))
a = ap.parse_args(); cfg = json.load(open(a.config)); jobs = json.load(open(a.jobs)); os.makedirs(a.outdir, exist_ok=True)
for tool in a.tools.split(","):
    t0 = time.time()
    try: res = ADAPTERS[tool](jobs, cfg, a.outdir)
    except Exception as ex: print(f"{tool}: not run ({ex})"); continue
    json.dump(res, open(os.path.join(a.outdir, f"{tool}.json"), "w"))
    print(f"{tool}: {sum(r['ok'] for r in res)}/{len(res)} analysed in {time.time() - t0:.0f}s")
