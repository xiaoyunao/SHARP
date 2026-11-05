"""Package the already computed nominal plots, with a date-browsable local index."""
import argparse
import hashlib
import json
from pathlib import Path
from PIL import Image
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader


def publish(source, bundle, output):
    index=json.loads((bundle/'night_index.json').read_text())
    nights=index['nights'] if isinstance(index,dict) else index
    dates=[r['date'] for r in nights]
    if len(dates)!=365 or len(set(dates))!=365:raise ValueError('Expected 365 unique nights')
    paths=[source/'pointings_aitoff.png']+sorted(source.glob('monthly_pointings_*.png'))
    if not paths[0].exists():
        paths=[p for p in source.glob('*.png') if p.name.startswith('annual')]+sorted(source.glob('monthly_pointings_*.png'))
    if len(paths)!=13:raise ValueError(f'Expected annual + 12 monthly plots, found {paths}')
    paths.extend(source/f'night_{d}.png' for d in dates)
    if any(not p.exists() for p in paths):raise ValueError('Missing plot: cannot fabricate zero nights')
    output.mkdir(parents=True,exist_ok=True);images=output/'images';images.mkdir(exist_ok=True)
    hashes={}
    for p in paths:
        with Image.open(p) as im:
            im=im.convert('RGB');im.thumbnail((1600,1200));target=images/(p.stem+'.jpg')
            im.save(target,quality=83,optimize=True)
        hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    pdf=canvas.Canvas(str(output/'nominal_annual_atlas.pdf'),pagesize=(900,650))
    pdf.setTitle('SHARP nominal simulated annual schedule — baseline a5cb945')
    for p in paths[:13]:
        pdf.drawImage(ImageReader(str(images/(p.stem+'.jpg'))),10,25,width=880,height=610,preserveAspectRatio=True,anchor='c')
        pdf.setFont('Helvetica',8);pdf.drawString(20,10,'Nominal simulated plan, not completed observations. Baseline a5cb945. Nightly plots: index.html.')
        pdf.showPage()
    pdf.save()
    links='\n'.join(f'<option value="{d}">{d}</option>' for d in dates)
    thumbs='\n'.join(f'<a href="images/{p.stem}.jpg"><img width="360" src="images/{p.stem}.jpg" alt="{p.stem}"></a>' for p in paths[:13])
    html='''<!doctype html><meta charset="utf-8"><title>SHARP nominal annual plan</title>
<style>body{font:16px system-ui;max-width:1400px;margin:30px auto;padding:12px}img{max-width:100%;height:auto}button,select{font:inherit;padding:8px}details img{margin:6px}a{color:#165ca3}</style>
<h1>SHARP：365夜名义模拟计划</h1><p>2026-11-01—2027-10-31；来源 a5cb945。不是已完成观测或天气预报。
13页PDF为全年/月度足迹；下方按北京时间夜次浏览全部365张高度与曝光图。</p>
<p><a href="nominal_annual_atlas.pdf">全年/月度PDF</a> · <a href="../baseline_plot_data/manifest.json">数据及输入哈希</a></p>
<details><summary>全年与12个月足迹</summary>THUMBS</details>
<p><button id="prev">前一夜</button> <select id="night">OPTIONS</select> <button id="next">后一夜</button></p>
<a id="full"><img id="figure" alt="所选夜次高度和曝光条带"></a>
<script>const s=document.getElementById('night');function show(){const u='images/night_'+s.value+'.jpg';document.getElementById('figure').src=u;document.getElementById('full').href=u;location.hash=s.value;}s.onchange=show;document.getElementById('prev').onclick=()=>{s.selectedIndex=Math.max(0,s.selectedIndex-1);show()};document.getElementById('next').onclick=()=>{s.selectedIndex=Math.min(s.length-1,s.selectedIndex+1);show()};if([...s.options].some(x=>x.value===location.hash.slice(1)))s.value=location.hash.slice(1);show();</script>'''.replace('OPTIONS',links).replace('THUMBS',thumbs)
    (output/'index.html').write_text(html)
    (output/'figure_manifest.json').write_text(json.dumps({'scope':'nominal simulation baseline a5cb945; not observed',
        'nights':dates,'input_png_sha256':hashes,'raster_delivery':'1600px JPEG quality83, raw PNGs retained in ignored products',
        'plot_manifest_sha256':hashlib.sha256((bundle/'manifest.json').read_bytes()).hexdigest()},indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('bundle',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();publish(a.source,a.bundle,a.output)
