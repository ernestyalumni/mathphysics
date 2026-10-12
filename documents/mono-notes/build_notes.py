#!/usr/bin/env python3
"""Build both editions in isolation, publishing PDFs and a warning report."""
from pathlib import Path
import argparse, hashlib, json, os, re, shutil, subprocess, tempfile
HERE=Path(__file__).resolve().parent

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=HERE)
    parser.add_argument('--repo',type=Path,default=HERE.parents[1])
    args=parser.parse_args()
    output=args.output.resolve()
    report={'editions':{}, 'inputs':{}}
    for p in sorted(HERE.rglob('*')):
        if p.is_file() and p.suffix in {'.tex','.py','.sh','.json'} and p.name!='build-report.json':
            report['inputs'][str(p.relative_to(HERE))]=digest(p)
    with tempfile.TemporaryDirectory(prefix='mathphysics-notes-') as tmp:
        root=Path(tmp)
        stage=root/'documents'/'mono-notes'
        shutil.copytree(HERE,stage,ignore=shutil.ignore_patterns('*.pdf','*.aux','*.log','*.toc','*.out','__pycache__','build-report.json'))
        (root/'LaTeX_and_pdfs').symlink_to(args.repo.resolve()/'LaTeX_and_pdfs',target_is_directory=True)
        for name in ['master','screen']:
            for attempt in range(3):
                result=subprocess.run(['pdflatex','-no-shell-escape','-interaction=nonstopmode','-halt-on-error',name+'.tex'],cwd=stage,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
                if result.returncode:
                    raise SystemExit(result.stdout.decode(errors='replace')[-4000:])
            log=(stage/(name+'.log')).read_text(errors='replace')
            warnings=re.findall(r'LaTeX Warning:.*?(?=\n\n)',log,re.S)
            overflow=re.findall(r'Overfull \\([hv])box \(([0-9.]+)pt .*?(?=\n\n)',log,re.S)
            report['editions'][name]={'pdf_sha256':digest(stage/(name+'.pdf')),
              'pages':int(re.search(r'Output written on .*?\((\d+) pages',log)[1]),
              'warnings':[' '.join(w.split()) for w in warnings],
              'overfull_boxes':[{'axis':axis,'points':float(points)} for axis,points in overflow]}
        output.mkdir(parents=True,exist_ok=True)
        for name in report['editions']:
            shutil.copy2(stage/(name+'.pdf'),output/(name+'.pdf'))
        (output/'build-report.json').write_text(json.dumps(report,indent=2)+'\n')
    for name, entry in report['editions'].items():
        print(f"{name}: {entry['pages']} pages; {len(entry['warnings'])} LaTeX warnings; {len(entry['overfull_boxes'])} overfull boxes")
    print('Published to',output)
if __name__=='__main__': main()
