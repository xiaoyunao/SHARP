#!/usr/bin/env python3
"""Render the selected L1 cutouts with frozen positions for the support audit."""
import gzip
import json

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from exposures import INPUT, ROOT, TABLES


def main():
    with gzip.open(INPUT/"followup_l1_pixel_sections.jsonl.gz","rt") as f:
        inputs=[json.loads(line) for line in f]
    stamps=[{"path":row["path"], **s} for row in inputs for s in row.get("stamps",[])]
    assert len(stamps)==len(inputs)==9
    # Upper row is the three-epoch control. The next rows show the comparisons
    # for three other fields, with the earlier exposure above the later one.
    stamps=[stamps[i] for i in (0,1,2,3,5,7,4,6,8)]
    fig,axes=plt.subplots(3,3,figsize=(7.7,8.5),layout="constrained")
    rows=[]
    colors=("#00e0ed","#ffad3c")
    for ax,s in zip(axes.flat,stamps):
        data=np.asarray(s["pixels"])
        assert data.shape==(41,41)
        center=np.nanmedian(data)
        scale=1.4826*np.nanmedian(np.abs(data-center))
        assert np.isfinite(scale) and scale>0
        stretch=np.arcsinh((data-center)/(3*scale))
        ax.imshow(stretch,origin="lower",cmap="gray",vmin=-.65,vmax=3.7,interpolation="nearest")
        ident,date,seq=s["label"].rsplit(" ",2)
        ax.set_title(f"{ident} · {date[:-4]}-{date[-4:-2]}-{date[-2:]} · {seq}",fontsize=9,pad=5)
        for j,m in enumerate(s["markers"]):
            x,y=m["pixel"]
            x,y=x-s["x0"],y-s["y0"]
            if j==0:
                ax.add_patch(plt.Circle((x,y),3.5,fill=False,color=colors[j],linewidth=1.2))
            else:
                ax.plot(x,y,marker="+",color=colors[j],markersize=9,mew=1.4)
            rows.append({"path":s["path"],"sample":s["label"],"header_date_obs":s["date_obs"],
                         "marker":m["label"],"ra_deg":m["ra"],"dec_deg":m["dec"],"stamp_x":x,"stamp_y":y})
        ax.set_xticks([])
        ax.set_yticks([])
    handles=[Line2D([],[],marker="o",mfc="none",mec=colors[0],ls="none",label="Historical matched sky position"),
             Line2D([],[],marker="+",color=colors[1],ls="none",label="Frozen prediction at each exposure")]
    if hasattr(fig, "set_layout_engine"):
        fig.legend(handles=handles,loc="outside lower center",ncol=2,fontsize=8,frameon=False)
    out=ROOT/"audit_figures"
    out.mkdir(exist_ok=True)
    for fmt in ("pdf","png"):
        fig.savefig(out/f"distant_l1_local_evidence.{fmt}",dpi=300)
    pd.DataFrame(rows).to_csv(TABLES/"distant_l1_stamp_manifest.csv",index=False)


if __name__ == "__main__":
    main()
