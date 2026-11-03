"""Download public inputs into ignored cache, with URL and content hashes."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import zipfile
from .catalog import BASE,download

SOURCES={
 'li2025.csv':'https://raw.githubusercontent.com/ZhuofuLi/ZTF-Jupiter-Trojan-Rotation-Periods/main/Summary_Table.csv',
 'li2025_license.txt':'https://raw.githubusercontent.com/ZhuofuLi/ZTF-Jupiter-Trojan-Rotation-Periods/main/LICENSE',
 'lcdb.zip':'https://minplanobs.org/MPInfo/datazips/LCLIST_PUB_CURRENT.zip',
 'models2026.html':'https://arxiv.org/html/2610.06082v1',
 'zenodo.json':'https://zenodo.org/api/records/22812438',
 'spin_results':'https://zenodo.org/api/records/22812438/files/spin_results/content',
 'spin_results_partial':'https://zenodo.org/api/records/22812438/files/spin_results_partial/content',
 'damit_asteroids.csv':'https://damit.cuni.cz/projects/damit/exports/table/asteroids',
 'damit_asteroid_models.csv':'https://damit.cuni.cz/projects/damit/exports/table/asteroid_models',
 'iers.all':'https://maia.usno.navy.mil/ser7/finals2000A.all',
 'tssys_README':'https://archive.konkoly.hu/pub/tssys/dr1/README',
 'tssys_release.info':'https://archive.konkoly.hu/pub/tssys/dr1/release.info',
 'tssys_release.merge':'https://archive.konkoly.hu/pub/tssys/dr1/release.merge',
 'tess_edr2.html':'https://arxiv.org/html/2503.13332v1',
}


def main():
    p=argparse.ArgumentParser();p.add_argument('--cache',default=str(BASE/'cache'));a=p.parse_args();cache=Path(a.cache)
    manifest=[]
    for name,url in SOURCES.items():
        try:
            path=download(url,cache/name)
            manifest.append({'file':name,'url':url,'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                             'bytes':path.stat().st_size,'cache_mtime_utc':datetime.fromtimestamp(path.stat().st_mtime,timezone.utc).isoformat(),'status':'downloaded'})
        except Exception as e:
            manifest.append({'file':name,'url':url,'status':'unavailable','error':str(e)})
        print(name,manifest[-1]['status'],flush=True)
    sbparams={'sb-class':'TJN','fields':'spkid,pdes,full_name,epoch,a,e,i,om,w,ma,H','full-prec':'true'}
    path=download('https://ssd-api.jpl.nasa.gov/sbdb_query.api',cache/'trojans.json',sbparams)
    manifest.append({'file':'trojans.json','url':'https://ssd-api.jpl.nasa.gov/sbdb_query.api','parameters':sbparams,
                     'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'status':'downloaded'})
    with zipfile.ZipFile(cache/'lcdb.zip') as z:(cache/'lcdb_summary.txt').write_bytes(z.read('lc_summary_pub.txt'))
    (cache/'sources_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if any(r['status']=='unavailable' for r in manifest):raise SystemExit('Some sources unavailable; inspect manifest; never classify their nonmatches as unknown periods')

if __name__=='__main__':main()
