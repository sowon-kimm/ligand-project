#!/usr/bin/env python3
"""Create capped-start poster copies with all original paths as a gray backdrop.

The original payload and energy/expandedEntries/energeticSpan/minimumPath
functions are embedded verbatim. No step-energy filtering is introduced.
"""
from pathlib import Path
import re,json,hashlib
ROOT=Path(__file__).resolve().parents[1]
PARTS=Path(__file__).parent

def build(facet):
    source=ROOT/f'html_files/capped_start_path_{facet}_corrected.html'
    original=source.read_text()
    scripts=re.findall(r'<script\b[^>]*>(.*?)</script>',original,re.S)
    app=scripts[-1]
    marker='const payload = '
    start=app.index(marker)+len(marker)
    payload,length=json.JSONDecoder().raw_decode(app[start:])
    functions=app[app.index('    function energy('):app.index('    function addLevel(')].strip()
    model='const payload = '+app[start:start+length]+';\n'
    model+='const stateByName = new Map(payload.states.map(state => [state.state,state]));\n'+functions
    page=(PARTS/'poster_capped_start/page.html.template').read_text()
    page=page.replace('Accessible Growth Pathways on Mn₃O₄ (001)',f'Capped-start Growth Pathways on Mn₃O₄ ({facet})')
    page=page.replace('grid-template-columns: 1fr 1fr; gap: 50px;','grid-template-columns: 1fr 1fr 1fr; gap: 28px;')
    page=page.replace('grid-template-columns: 125px 1fr 82px; gap: 14px;','grid-template-columns: 95px 1fr 72px; gap: 10px;')
    oa='<label><span>ΔμOA</span><input aria-label="Oleic acid chemical potential" data-role="mu-oa" type="range" min="-4" max="0" step="0.05" value="0"><output data-role="mu-oa-value">0.00 eV</output></label>'
    page=page.replace('</div>\n<div class="exports">',oa+'\n</div>\n<div class="exports">')
    for key,value in {'__PLOTLY_LIBRARY__':scripts[1],'__MODEL__':model,
                      '__VISUALIZATION__':(PARTS/'poster_capped_start/presentation.js').read_text(),
                      '__CONTROLS__':''}.items():page=page.replace(key,value)
    output=ROOT/f'html_files/capped_start_path_{facet}_corrected_poster.html'
    output.write_text(page)
    output.with_suffix('.provenance.json').write_text(json.dumps({
        'source':source.name,'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'paths':len(payload['paths']),'states':len(payload['states']),
        'calculation':'Payload and original energy, expandedEntries, energeticSpan, minimumPath copied verbatim.',
        'background':'Union of all existing expanded paths at their original entry coordinates; identical segments drawn once.',
        'filter':'No additional filter.'},indent=2)+'\n')
    print(output)

if __name__=='__main__':
    for facet in ['001','011']:build(facet)
