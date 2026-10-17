#!/usr/bin/env python3
"""Import Ernest's notes and solutions for John M. Lee, Introduction to Smooth Manifolds.

Source: LeeJM/LeeJM_IntroSmoothManifolds_sol.tex, a standalone amsart document that
\\input's one file per chapter and carries several chapters inline. The manuscript is not
modified. Output: legacy/lee/lee-NN.tex, one file per original \\section, in document order,
plus legacy/lee/provenance.json. master.tex inputs them after the Manifolds part.

Kept separate from import_legacy.py on purpose: re-running that script regenerates the
agat/dgdt files and would discard later edits to them.
"""
from pathlib import Path
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC_DIR = ROOT / 'LeeJM'
MAIN = SRC_DIR / 'LeeJM_IntroSmoothManifolds_sol.tex'
PART_TITLE = r'Smooth Manifolds: notes and solutions for Lee, \emph{Introduction to Smooth Manifolds}'
APPENDIX_LETTERS = {'Review of Topology': 'A', 'Review of Linear Algebra': 'B', 'Review of Calculus': 'C'}


def strip_comments(text):
    # Comments are not mathematics (the chapter files open with long banner comments).
    out = []
    for line in text.splitlines():
        if line.lstrip().startswith('%'):
            continue
        out.append(re.sub(r'(?<!\\)%[^\n]*', '', line) if not (r'\url{' in line or r'\href{' in line) else line)
    return '\n'.join(out)


def expand_inputs(body, used):
    def repl(m):
        name = m[1].strip()
        path = SRC_DIR / (name if name.endswith('.tex') else name + '.tex')
        used.append(path)
        return expand_inputs(strip_comments(path.read_text()), used)
    return re.sub(r'\\input\{([^}]*)\}', repl, body)


def tidy_math(body):
    # Blank lines inside display math are errors in the master's amsmath setup.
    squash = lambda m: re.sub(r'\n(?:[ \t]*\n)+', '\n', m[0])
    body = re.sub(r'\\\[.*?\\\]', squash, body, flags=re.S)
    body = re.sub(r'\\begin\{(equation\*?|align\*?|gather\*?|gathered|aligned|multline\*?)\}.*?\\end\{\1\}',
                  squash, body, flags=re.S)
    return re.sub(r'\\begin\{tikzpicture\}.*?\\end\{tikzpicture\}', squash, body, flags=re.S)


def main():
    raw = MAIN.read_text()
    body = raw.split(r'\begin{document}', 1)[1].split(r'\end{document}', 1)[0]
    body = strip_comments(body)
    # Drop the title page, contact block, acknowledgments and \tableofcontents.
    body = body[body.index(r'\section{'):]
    used = []
    body = expand_inputs(body, used)
    body = re.sub(r'\\(?:begin|end)\{appendix\}', '', body)
    for title, letter in APPENDIX_LETTERS.items():
        body = re.sub(r'\\section\{\s*' + re.escape(title) + r'\s*\}',
                      r'\\section{Appendix ' + letter + ': ' + title + '}', body, count=1)
    body = body.replace(r'\righthalfcup', r'\lrcorner')        # MnSymbol interior-product hook
    body = body.replace(r'\staveXXIX', '')                      # decorative glyph (staves package)
    body = re.sub(r'\\(label|ref|eqref|pageref|autoref)\{([^}]*)\}', r'\\\1{lee:\2}', body)
    body = tidy_math(body)

    dest = HERE / 'legacy' / 'lee'
    dest.mkdir(parents=True, exist_ok=True)
    starts = [m.start() for m in re.finditer(r'\\section\{', body)]
    records = []
    for i, s in enumerate(starts):
        chunk = body[s:starts[i + 1] if i + 1 < len(starts) else len(body)].strip() + '\n'
        title = re.match(r'\\section\{([^}]*)\}', chunk)[1].strip()
        slug = f'lee-{i + 1:02d}'
        head = '% Imported from LeeJM/ by import_lee.py; see legacy/lee/provenance.json.\n'
        if i == 0:
            head += r'\part{' + PART_TITLE + '}\n'
        (dest / f'{slug}.tex').write_text(head + chunk)
        records.append({'id': slug, 'title': title, 'output': f'legacy/lee/{slug}.tex'})
    sources = [MAIN] + used
    (dest / 'provenance.json').write_text(json.dumps({
        'transformations': [
            'Drop the title page, contact block, acknowledgments, table of contents and comments',
            'Expand \\input chapter files in document order',
            'Split at \\section; the first file opens the part',
            'Replace the appendix environment by sections titled Appendix A/B/C',
            'MnSymbol \\righthalfcup (interior product) -> amssymb \\lrcorner; drop decorative \\staveXXIX',
            'Prefix labels and references with lee:',
            'Remove blank paragraphs inside display math and TikZ',
        ],
        'sources': [{'path': str(p.relative_to(ROOT)), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
                    for p in sources],
        'parts': records,
    }, indent=2) + '\n')
    for r in records:
        print(r['id'], r['title'])
    return records


if __name__ == '__main__':
    main()
