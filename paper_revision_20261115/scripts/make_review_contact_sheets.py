#!/usr/bin/env python3
"""Server read-only GIF previews for the explicit annual completion-evidence queue."""
import argparse,csv,math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import pandas as pd
from PIL import Image,ImageDraw,ImageOps
ROOT=Path('/pipeline/xiaoyunao')

def preview(row):
 night=str(row.night);folder=ROOT/'heliolincrr/review_packages'/night/'gifs'
 candidates=list(folder.glob(f'{row.trk_sub}_link*.gif'))
 p=candidates[0] if candidates else ROOT/'heliolincrr/plots'/night/f'unknown_link_{int(row.linkage_id):04d}_{night}.gif'
 card=Image.new('RGB',(480,170),'white');draw=ImageDraw.Draw(card)
 draw.text((5,4),f'{night}  {row.trk_sub}  link {int(row.linkage_id)}',fill='black')
 state='existing_gif'
 try:
  with Image.open(p) as gif:
   # Three evenly distributed temporal views. No image measurements are altered.
   for j,fr in enumerate([0,gif.n_frames//2,gif.n_frames-1]):
    gif.seek(fr);im=gif.convert('RGB');im.thumbnail((154,135));card.paste(im,(j*160+(154-im.width)//2,28))
 except Exception:
  draw.text((5,60),'No preserved GIF available; use catalog/member index.',fill='black');state='gif_missing'
 return card,dict(night=night,trk_sub=row.trk_sub,linkage_id=int(row.linkage_id),source_gif=str(p),preview_state=state)

def main():
 a=argparse.ArgumentParser();a.add_argument('--queue',type=Path,required=True);a.add_argument('--output',type=Path,required=True);x=a.parse_args();x.output.mkdir(parents=True,exist_ok=True)
 d=pd.read_csv(x.queue,dtype={'night':str});rows=list(d.sort_values(['night','linkage_id']).itertuples());idx=[];sheet=None
 with ThreadPoolExecutor(max_workers=4) as pool:
  for i,(card,entry) in enumerate(pool.map(preview,rows)):
   page=i//48+1
   if i%48==0:
    if sheet is not None:sheet.save(x.output/f'review_{page-1:03d}.jpg',quality=85)
    sheet=Image.new('RGB',(1440,16*170+35),'#e1e1e1');ImageDraw.Draw(sheet).text((8,8),f'SHARP annual review-completion queue — page {page}; three existing GIF views per link',fill='black')
   sheet.paste(card,((i%3)*480,35+((i%48)//3)*170));entry['sheet']=f'review_{page:03d}.jpg';entry['card_1based']=i%48+1;idx.append(entry)
 if sheet is not None:sheet.save(x.output/f'review_{len(rows)//48+bool(len(rows)%48):03d}.jpg',quality=85)
 pd.DataFrame(idx).to_csv(x.output/'review_preview_index.csv',index=False)
 print('previewed',len(idx),'missing',sum(r['preview_state']=='gif_missing' for r in idx),flush=True)
if __name__=='__main__':main()
