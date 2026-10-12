#!/usr/bin/env python3
"""Import archival mathematics without modifying the historical manuscripts."""
from pathlib import Path
import hashlib, json, re
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

def closing(text, opening):
    depth = 1
    for pos in range(opening + 1, len(text)):
        if text[pos] == '{' and (pos == 0 or text[pos-1] != '\\'): depth += 1
        if text[pos] == '}' and (pos == 0 or text[pos-1] != '\\'):
            depth -= 1
            if depth == 0: return pos
    raise ValueError('Unbalanced part heading')

def main(root=ROOT):
    dest = HERE / 'legacy' / 'mathematics'
    dest.mkdir(parents=True, exist_ok=True)
    records, bibliographies = [], []
    for prefix, filename in [('agat', 'AGAT_dump.tex'), ('dgdt', 'DGDT_dump.tex')]:
        source = root / 'LaTeX_and_pdfs' / filename
        raw = source.read_text()
        body = raw.split(r'\begin{document}', 1)[1].split(r'\end{document}', 1)[0]
        # Comments are not mathematics; strip them to distinguish inactive layouts.
        body = '\n'.join(line if (r'\url{' in line or r'\href{' in line) and not line.lstrip().startswith('%') else re.sub(r'(?<!\\)%[^\n]*', '', line) for line in body.splitlines())
        body = re.sub(r'\\(?:begin\{multicols\*?\}\{\d+\}|end\{multicols\*?\})', '', body)
        body = re.sub(r'\\(?:label|ref|eqref|pageref|autoref)\{([^}]*)\}',
                      lambda m: m[0].replace('{'+m[1]+'}', '{'+prefix+':'+m[1]+'}'), body)
        body = re.sub(r'\\cite(\[[^\]]*\])?\{([^}]*)\}',
                      lambda m: r'\cite'+(m[1] or '')+'{'+','.join(prefix+':'+k.strip() for k in m[2].split(','))+'}', body)
        body = re.sub(r'\\bibitem(\[[^\]]*\])?\{([^}]*)\}',
                      lambda m: r'\bibitem'+(m[1] or '')+'{'+prefix+':'+m[2]+'}', body)
        body = re.sub(r'\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}', lambda m: re.sub(r'\n(?:[ \t]*\n)+', '\n', m[0]), body, flags=re.S)
        body = body.replace(r'\left. &', r'& \left.')
        body = re.sub(r'\\\[.*?\\\]', lambda m: re.sub(r'\n(?:[ \t]*\n)+', '\n', m[0]), body, flags=re.S)
        body = re.sub(r'\\begin\{(equation\*?|align\*?|gather\*?)\}.*?\\end\{\1\}', lambda m: re.sub(r'\n(?:[ \t]*\n)+', '\n', m[0]), body, flags=re.S)
        body = body.replace('agat:SubSubSect:RecordTypesProduct', 'agat:SubSubSec:RecordTypesProduct')
        body = body.replace('agat:Eq::RightMPSoperation', 'agat:Eq:RightMPSoperation')
        # Repeated worked passages keep distinct targets; refs use the nearest preceding definition.
        counts, targets = {}, {}
        def distinct(match):
            command, key = match[1], match[2]
            if command == 'label':
                counts[key] = counts.get(key, 0) + 1
                targets[key] = key if counts[key] == 1 else key + ':' + str(counts[key])
            return '\\' + command + '{' + targets.get(key, key) + '}'
        body = re.sub(r'\\(label|ref|eqref|pageref|autoref)\{([^}]+)\}', distinct, body)
        bibstart = body.index(r'\begin{thebibliography}')
        bib = body[bibstart:].split(r'\end{thebibliography}', 1)[0]
        bib = re.sub(r'^\\begin\{thebibliography\}\{[^}]*\}', '', bib)
        bibliographies.append(bib)
        body = body[:bibstart]
        parts = list(re.finditer(r'\\part\{', body))
        for i, match in enumerate(parts):
            end = closing(body, match.end()-1)
            title = body[match.end():end]
            chunk = body[end+1:parts[i+1].start() if i+1<len(parts) else len(body)]
            slug = f'{prefix}-{i+1:02d}'
            path = dest / (slug+'.tex')
            path.write_text('% Imported from LaTeX_and_pdfs/'+filename+'; see provenance.json.\n'+r'\part{'+title+'}\n'+chunk)
            records.append({'id':slug,'title':title,'source':'LaTeX_and_pdfs/'+filename,
                            'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                            'output':'legacy/mathematics/'+slug+'.tex'})
    (dest/'bibliography.tex').write_text('\n'.join(bibliographies)+'\n')
    (dest/'provenance.json').write_text(json.dumps({'transformations':['Drop historical front matter and comments','Use master layout in place of embedded multicols','Prefix labels, references, citations and bibliography keys by source','Split at original part headings; preserve mathematical bodies', 'Remove blank paragraphs inside math and TikZ environments', 'Move alignment tab outside left delimiter in set-theory display', 'Repair RecordTypesProduct and RightMPSoperation reference typos', 'Disambiguate repeated equation labels and their subsequent references'], 'parts':records},indent=2)+'\n')
    return records

if __name__ == '__main__':
    import argparse
    parser=argparse.ArgumentParser(); parser.add_argument('--repo',type=Path,default=ROOT)
    main(parser.parse_args().repo)
