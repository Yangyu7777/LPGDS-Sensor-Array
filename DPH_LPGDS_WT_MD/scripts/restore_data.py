"""Restore original large files from gzip parts and verify exact SHA256 (standard library only)."""
from pathlib import Path
import json,hashlib,gzip,tempfile,shutil,sys
ROOT=Path(__file__).resolve().parents[1]
def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 manifest=ROOT/'LARGE_FILES.json'
 if not manifest.exists():raise SystemExit('LARGE_FILES.json is absent; use the complete repository download.')
 verify_only='--verify-only' in sys.argv
 for item in json.loads(manifest.read_text())['files']:
  target=ROOT/item['original_path']
  if not verify_only and target.exists():
   if digest(target)!=item['original_sha256']:raise SystemExit('Existing file has wrong SHA256; preserved without overwrite: '+str(target))
   print('Already restored:',item['original_path']);continue
  with tempfile.TemporaryFile() as packed:
   for part in item['parts']:
    f=ROOT/part['path']
    if digest(f)!=part['sha256']:raise SystemExit('Part checksum mismatch: '+str(f))
    with f.open('rb') as chunk:shutil.copyfileobj(chunk,packed,1024*1024)
   packed.seek(0);h=hashlib.sha256();size=0;tmp=None;out=None
   if not verify_only:target.parent.mkdir(parents=True,exist_ok=True);tmp=target.with_name(target.name+'.restoring');out=tmp.open('wb')
   try:
    with gzip.GzipFile(fileobj=packed,mode='rb') as gz:
     for b in iter(lambda:gz.read(1024*1024),b''):
      h.update(b);size+=len(b)
      if out:out.write(b)
   finally:
    if out:out.close()
   if h.hexdigest()!=item['original_sha256'] or size!=item['original_bytes']:raise SystemExit('Original checksum/size mismatch; no original overwritten: '+str(target))
   if tmp:tmp.replace(target)
   print('Verified original:',item['original_path'],size,'bytes')
 print('All large files verified'+(' (no files written).' if verify_only else ' and restored.'))
if __name__=='__main__':main()
