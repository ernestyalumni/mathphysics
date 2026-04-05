"""
parse_equations.py
==================
Parse LaTeX output from Cadabra2/Srednicki export scripts to extract equations.
Usage: python3 parse_equations.py
"""

import os
import re
import json
from pathlib import Path
from typing import List, Dict
from dataclasses import dataclass

@dataclass
class Equation:
    label: str
    latex: str
    description: str = ""
    section: str = ""
    subsection: str = ""

@dataclass  
class Chapter:
    id: str
    number: str
    title: str
    status: str
    description: str
    equations: List[Equation]

# Srednicki directory
SR_DIR = Path(__file__).parent

def extract_equations_from_latex(latex_path: Path) -> List[Equation]:
    """Extract equations and labels from a LaTeX chapter file."""
    equations = []
    
    with open(latex_path, 'r') as f:
        content = f.read()
    
    current_section = ""
    current_subsection = ""
    
    # Extract section/subsection markers
    for match in re.finditer(r'\\section\{([^}]+)\}', content):
        current_section = match.group(1)
    
    for match in re.finditer(r'\\subsection\{([^}]+)\}', content):
        current_subsection = match.group(1)
    
    # Find all equations with tags
    for match in re.finditer(r'\\tag\{([^}]+)\}', content):
        eq_num = match.group(1)
        eq_num_str = str(eq_num)
        
        # Get context (200 chars before and after)
        pos = match.start()
        start = max(0, pos - 400)
        end = min(len(content), pos + 400)
        context = content[start:end]
        
        # Extract latex from boxed equation or surrounding text
        # Look for \boxed{...}
        boxed_match = re.search(r'\\boxed\{([^}]+(?:\{[^}]*\}[^}]*)*)\}', context)
        if boxed_match:
            latex = boxed_match.group(1)
        else:
            # Fallback: take the first meaningful math expression
            # Look for \begin{equation} ... \end{equation}
            eqn_start = context.find('\\begin{equation}')
            eqn_end = context.find('\\end{equation}')
            if eqn_start >= 0 and eqn_end > eqn_start:
                latex = context[eqn_start:eqn_end]
            else:
                # Take the tagged line itself
                latex = context[max(0, pos-100):pos+100]
        
        # Clean up latex
        latex = latex.strip()
        
        equations.append(Equation(
            label=eq_num_str,
            latex=latex,
            description="",
            section=current_section,
            subsection=current_subsection
        ))
    
    return equations

def parse_all_chapters() -> List[Chapter]:
    """Parse all Srednicki chapter LaTeX files."""
    chapters = []
    
    chapter_configs = {
        "ch34": {"number": "34", "title": "Left/Right Weyl Spinors", "status": "done"},
        "ch35": {"number": "35", "title": "Majorana and Dirac Spinors", "status": "done"},
        "ch36": {"number": "36", "title": "Weyl Lagrangian & Symmetries", "status": "done"},
        "ch37": {"number": "37", "title": "Canonical Quantization of Weyl Fields", "status": "in-progress"},
        "ch38": {"number": "38", "title": "LSZ for Spinors & Feynman Rules", "status": "in-progress"},
        "ch60": {"number": "60", "title": "MHV Amplitudes (Parke-Taylor)", "status": "done"},
        "bcfw": {"number": "BCFW", "title": "BCFW Recursion Relations", "status": "not-started"},
        "adscft": {"number": "AdS/CFT", "title": "AdS/CFT Correspondence", "status": "not-started"},
    }
    
    for chapter_id, config in chapter_configs.items():
        latex_file = SR_DIR / f"{chapter_id}.tex"
        if not latex_file.exists():
            continue
            
        equations = extract_equations_from_latex(latex_file)
        
        chapter = Chapter(
            id=chapter_id,
            number=config["number"],
            title=config["title"],
            status=config["status"],
            description=chapter_id,
            equations=equations
        )
        chapters.append(chapter)
    
    return chapters

if __name__ == "__main__":
    chapters = parse_all_chapters()
    
    print(f"Parsed {len(chapters)} chapters")
    for ch in chapters:
        print(f"  {ch.id}: {len(ch.equations)} equations")
    
    # Save to JSON for API
    output_path = SR_DIR / "parsed_results.json"
    with open(output_path, 'w') as f:
        json.dump(chapters, f, indent=2)
    print(f"\nSaved to: {output_path}")
