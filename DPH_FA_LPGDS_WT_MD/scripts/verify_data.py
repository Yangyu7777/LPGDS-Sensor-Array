"""Verify original FA_R2 files against the public checksum manifest."""
from pathlib import Path
import csv,hashlib,sys
root=Path(__file__).resolve().parents[1];errors=[];count=0
with (root/'DATA_SHA256.csv').open() as f:
 for row in csv.DictReader(f):
  p=root/row['relative_path']
  if not p.exists():errors.append('Missing (restore first): '+row['relative_path']);continue
  h=hashlib.sha256()
  with p.open('rb') as src:
   for b in iter(lambda:src.read(1024*1024),b''):h.update(b)
  if h.hexdigest()!=row['sha256'] or p.stat().st_size!=int(row['bytes']):errors.append('Mismatch: '+row['relative_path'])
  count+=1
if errors:print('\n'.join(errors));sys.exit(1)
print('All',count,'original FA_R2 files passed SHA256 verification.')
