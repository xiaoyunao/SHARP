"""Package the already computed nominal plots, with a date-browsable local index."""
import argparse
import hashlib
import json
import os
from html import escape
from pathlib import Path
from PIL import Image
from pypdf import PdfWriter


def publish(source, bundle, output):
    manifest_path = bundle / 'manifest.json'
    manifest = json.loads(manifest_path.read_text())
    source_commit = manifest['source']['commit']
    count = manifest['counts']['exposures']
    label = {'4605b372': 'B', 'a5cb9454': 'A'}.get(source_commit[:8], source_commit[:8])
    provenance = json.loads((source / 'plot_provenance.json').read_text())
    bundle_hash = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    if provenance['source_manifest_sha256'] != bundle_hash or provenance['source_commit'] != source_commit:
        raise ValueError('Figure provenance disagrees with bundle')
    index=json.loads((bundle/'night_index.json').read_text())
    nights=index['nights'] if isinstance(index,dict) else index
    dates=[r['date'] for r in nights]
    if len(dates)!=365 or len(set(dates))!=365:raise ValueError('Expected 365 unique nights')
    paths=[source/'pointings_aitoff.png']+sorted(source.glob('monthly_pointings_*.png'))
    if not paths[0].exists():
        paths=[p for p in source.glob('*.png') if p.name.startswith('annual')]+sorted(source.glob('monthly_pointings_*.png'))
    if len(paths)!=13:raise ValueError(f'Expected annual + 12 monthly plots, found {paths}')
    paths.extend(source/f'night_{d}.png' for d in dates)
    if [r['date'] for r in provenance['nights_plotted']] != dates:
        raise ValueError('Figure provenance does not cover the complete night index')
    if any(not p.exists() or not p.with_suffix('.pdf').exists() for p in paths):
        raise ValueError('Missing PNG/PDF plot')
    if sum(r['exposure_rows'] for r in provenance['nights_plotted']) != count:
        raise ValueError('Exposure count does not match bundle')
    output.mkdir(parents=True,exist_ok=True);images=output/'images';images.mkdir(exist_ok=True)
    hashes={}
    for p in paths:
        with Image.open(p) as im:
            im=im.convert('RGB');im.thumbnail((1600,1200));target=images/(p.stem+'.jpg')
            im.save(target,quality=83,optimize=True)
        hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    title = f'SHARP {label} nominal plan · {count:,} exposures'
    for filename, pages in [('nominal_annual_atlas.pdf', paths[:13]), ('nominal_nightly_atlas.pdf', paths[13:])]:
        writer = PdfWriter()
        for p in pages:
            writer.append(str(p.with_suffix('.pdf')), import_outline=False)
        writer.add_metadata({'/Title': title, '/Subject': f'Plan source {source_commit}'})
        writer.write(str(output / filename))
        writer.close()
    options = []
    for night in nights:
        day = night['date']
        has_trojan = any(block['mode'] == 'trojan' for block in night['blocks'])
        reason = 'Trojan 已排程' if has_trojan else '无 Trojan 块 · ' + ('; '.join(night['skips']) or '源计划无 Trojan 块')
        options.append(f'<option value="{day}" data-reason="{escape(reason, quote=True)}">{day}</option>')
    links = '\n'.join(options)
    bundle_link = escape(os.path.relpath(manifest_path, output), quote=True)
    thumbs='\n'.join(f'<a href="images/{p.stem}.jpg"><img width="360" src="images/{p.stem}.jpg" alt="{p.stem}"></a>' for p in paths[:13])
    html='''<!doctype html><meta charset="utf-8"><title>SHARP nominal annual plan</title>
<style>body{font:16px system-ui;max-width:1400px;margin:30px auto;padding:12px}img{max-width:100%;height:auto}button,select{font:inherit;padding:8px}details img{margin:6px}a{color:#165ca3}</style>
<h1>SHARP · PLAN_LABEL 名义年</h1><p>2026-11-01 至 2027-10-31，共 FRAME_COUNT 条名义曝光。排程源提交 SOURCE_SHA。
13页PDF为全年/月度足迹；下方按北京时间夜次浏览全部365张高度与曝光图。</p>
<p><a href="nominal_annual_atlas.pdf">全年/月度PDF</a> · <a href="nominal_nightly_atlas.pdf">365夜 PDF</a> · <a href="BUNDLE_LINK">数据及输入哈希</a></p>
<details><summary>全年与12个月足迹</summary>THUMBS</details>
<p><button id="prev">前一夜</button> <select id="night">OPTIONS</select> <button id="next">后一夜</button></p>
<p id="reason"></p><a id="full"><img id="figure" alt="所选夜次高度和曝光条带"></a>
<script>const s=document.getElementById('night');function show(){const u='images/night_'+s.value+'.jpg';document.getElementById('figure').src=u;document.getElementById('full').href=u;document.getElementById("reason").textContent=s.selectedOptions[0].dataset.reason;location.hash=s.value;}s.onchange=show;document.getElementById('prev').onclick=()=>{s.selectedIndex=Math.max(0,s.selectedIndex-1);show()};document.getElementById('next').onclick=()=>{s.selectedIndex=Math.min(s.length-1,s.selectedIndex+1);show()};if([...s.options].some(x=>x.value===location.hash.slice(1)))s.value=location.hash.slice(1);show();</script>'''.replace('OPTIONS',links).replace('THUMBS',thumbs).replace('PLAN_LABEL',label).replace('FRAME_COUNT',f'{count:,}').replace('SOURCE_SHA',source_commit).replace('BUNDLE_LINK',bundle_link)
    (output/'index.html').write_text(html)
    (output/'figure_manifest.json').write_text(json.dumps({'scope':'nominal simulation', 'strategy':label, 'plan_source_commit':source_commit,
        'exposures':count, 'publication_source':'Git commit containing this manifest',
        'nights':dates,'input_png_sha256':hashes,'raster_delivery':'1600px JPEG quality83, raw PNGs retained in ignored products',
        'plot_manifest_sha256':bundle_hash, 'plot_provenance':provenance},indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('source',type=Path);p.add_argument('bundle',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();publish(a.source,a.bundle,a.output)
