#!/usr/bin/env python3
"""Export the node/edge data behind the reused full SHARP architecture figure."""
from pathlib import Path
import argparse, importlib.util, sys
import pandas as pd

EDGES=[
 ("scheduler","observing","solid"),("observing","calibration","solid"),("calibration","nightly_products","solid"),("scheduler_state","scheduler","dotted"),
 ("nightly_products","known_query","solid"),("nightly_products","gaia_mask","solid"),("known_query","known_match","solid"),("known_match","known_store","solid"),("known_store","known_ades","solid"),
 ("gaia_mask","tracklets_links","solid"),("tracklets_links","orbit_subtract","solid"),("orbit_subtract","review_package","solid"),("review_package","human_review","solid"),("human_review","unknown_ades","solid"),("review_state","review_package","dotted"),
 ("known_ades","mpc","solid"),("unknown_ades","mpc","solid"),("human_review","followup_state","dashed"),("followup_state","scheduler","dashed feedback"),
]

def main(a):
    script=Path(a.architecture_script).resolve();sys.path.insert(0,str(script.parent));spec=importlib.util.spec_from_file_location("sharp_architecture",script);mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod);nodes=mod.build_nodes();root=Path(a.root).resolve()
    nd=pd.DataFrame([{"node_id":n.key,"label":n.label.replace("\n"," | "),"x":n.x,"y":n.y,"width":n.width,"height":n.height,"category":n.category,"shape":n.shape,"source_code":"paper_analysis_20260803/scripts/make_fig01_architecture.py"} for n in nodes.values()]);ed=pd.DataFrame(EDGES,columns=["source_node","target_node","edge_style"]);nd.to_csv(root/"figure_data/workflow_nodes.csv",index=False);ed.to_csv(root/"figure_data/workflow_edges.csv",index=False)
    pd.concat([nd.assign(row_type="node"),ed.assign(row_type="edge")],ignore_index=True,sort=False).to_csv(root/"figure_data/workflow_nodes_edges.csv",index=False)
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--root",required=True);p.add_argument("--architecture-script",required=True);main(p.parse_args())
