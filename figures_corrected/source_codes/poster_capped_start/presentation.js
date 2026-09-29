// Presentation only. All energy/span/path-selection functions above are unchanged.
const nodeKey = e => `${e.x}|${e.state}|${e.kind}`;
const edgeKey = (a, b) => `${nodeKey(a)}>${nodeKey(b)}`;
const allNodes = new Map(), allEdges = new Map();
// Keep each path's original expanded-entry coordinates. Deduplicate identical
// drawn segments so shared edges do not darken through repeated overplotting.
for (const path of payload.paths) {
    const entries = expandedEntries(path, 0, 0, 0);
    for (const e of entries) allNodes.set(nodeKey(e), {x: e.x, state: e.state, kind: e.kind});
    for (let i = 1; i < entries.length; i++) {
        const a = entries[i - 1], b = entries[i];
        allEdges.set(edgeKey(a, b), {a: nodeKey(a), b: nodeKey(b)});
    }
}
const extentX = Math.max(...[...allNodes.values()].map(e => e.x));
const styles = {desorption: {name: 'Desorption', color: '#334155'},
    OL: {name: 'Olation', color: '#2563eb'}, OX: {name: 'Oxolation', color: '#16a34a'},
    C: {name: 'Condensation', color: '#9333ea'}, SR: {name: 'Redox', color: '#dc2626'}};
const config = {responsive: false, displayModeBar: false, displaylogo: false, scrollZoom: false};

function makeCappedPoster(muMn, muH2O, muOA) {
    const result = minimumPath(muMn, muH2O, muOA);
    const entries = expandedEntries(result.path, muMn, muH2O, muOA);
    const nodes = new Map([...allNodes].map(([key, e]) => [key, {...e,
        energy: energy(e.state, e.kind === 'capped' ? 'carboxyl' : 'pure', muMn, muH2O, muOA)}]));
    const selectedNodes = new Set(entries.map(nodeKey));
    const selectedEdges = new Set(entries.slice(1).map((b, i) => edgeKey(entries[i], b)));
    const segments = edges => {
        const x = [], y = [];
        for (const [a,b] of edges) {
            x.push(a.x + payload.levelHalfWidth, b.x - payload.levelHalfWidth, null);
            y.push(a.energy, b.energy, null);
        }
        return {x, y};
    };
    const levels = list => {
        const x = [], y = [];
        for (const e of list) {
            x.push(e.x - payload.levelHalfWidth, e.x + payload.levelHalfWidth, null);
            y.push(e.energy, e.energy, null);
        }
        return {x,y};
    };
    const backgroundNodes = [...nodes].filter(([k]) => !selectedNodes.has(k)).map(([,e]) => e);
    const backgroundEdges = [...allEdges].filter(([k]) => !selectedEdges.has(k))
        .map(([,e]) => [nodes.get(e.a), nodes.get(e.b)]);
    const data = [{...segments(backgroundEdges), type:'scatter', mode:'lines',
        name:'Other pathways', legendrank:6, line:{color:'#b4bac2',width:0.7}, opacity:0.22, hoverinfo:'skip'},
        {...levels(backgroundNodes), type:'scatter', mode:'lines', showlegend:false,
            line:{color:'#b4bac2',width:0.7}, opacity:0.22, hoverinfo:'skip'},
        {x:backgroundNodes.map(e=>e.x), y:backgroundNodes.map(e=>e.energy),
            text:backgroundNodes.map(e=>`${e.state} (${e.kind})`),type:'scatter',mode:'markers',showlegend:false,
            marker:{color:'#b4bac2',size:3},opacity:0.22,
            hovertemplate:'%{text}<br>Expanded step: %{x}<br>ΔG = %{y:.4f} eV<extra></extra>'}];
    for (const [i, kind] of ['desorption','OL','OX','C','SR'].entries()) {
        const edges = entries.slice(1).flatMap((b,j) =>
            (b.kind === 'pure' ? 'desorption' : b.mechanism) === kind ? [[entries[j],b]] : []);
        const xy = segments(edges);
        data.push({x:xy.x.length?xy.x:[null],y:xy.y.length?xy.y:[null],type:'scatter',mode:'lines',
            name:styles[kind].name,legendrank:i+1,line:{color:styles[kind].color,width:3.5},opacity:1,hoverinfo:'skip'});
    }
    data.push({...levels(entries),type:'scatter',mode:'lines',line:{color:'#20252b',width:2.7},showlegend:false,hoverinfo:'skip'});
    data.push({x:entries.map(e=>e.x),y:entries.map(e=>e.energy),text:entries.map(e=>`${e.state} (${e.kind})`),
        type:'scatter',mode:'markers',marker:{color:'#20252b',size:6},showlegend:false,
        hovertemplate:'%{text}<br>Expanded step: %{x}<br>ΔG = %{y:.4f} eV<extra>Minimum-span pathway</extra>'});
    const pair=result.pair, annotations=[], shapes=[];
    const addLabel = (entry,text,ay,ax=0,bold=false) => annotations.push({x:entry.x,y:entry.energy,
        xref:'x',yref:'y',text:bold?`<b>${text}</b>`:text,showarrow:true,arrowhead:0,
        arrowwidth:0.7,arrowcolor:'#64748b',ax,ay,font:{size:bold?15:13,color:'#20252b'},
        bgcolor:'rgba(255,255,255,0.94)',borderpad:2});
    const controlling = new Set([nodeKey(pair.tdi),nodeKey(pair.tdts)]);
    const capped = entries.filter(e=>e.kind==='capped');
    const labelStride = capped.length > 10 ? 2 : 1;
    capped.forEach((entry,i) => {
        if (controlling.has(nodeKey(entry))) return;
        if (i%labelStride===0 || i===capped.length-1) addLabel(entry,`${entry.state} (cap)`, i%2 ? 29 : -29);
    });
    for (const [kind,entry,symbol,ay] of [['TDI',pair.tdi,'circle',33],['TDTS*',pair.tdts,'triangle-up',-33]]) {
        data.push({x:[entry.x],y:[entry.energy],text:[`${kind}: ${entry.state} (${entry.kind})`],
            type:'scatter',mode:'markers',name:kind,showlegend:false,
            marker:{symbol,size:kind==='TDI'?15:13,color:'#111827',line:{color:'white',width:1.1}},
            hovertemplate:'%{text}<br>ΔG = %{y:.4f} eV<extra></extra>'});
        addLabel(entry,`${kind}: ${entry.state}<br>(${entry.kind})`,ay,entry.x===0?35:0,true);
    }
    const bracketX=extentX+0.8, low=pair.tdi.energy,
        high=pair.tdts.energy+pair.appliedReactionEnergy;
    const line=(x0,y0,x1,y1,dash='dot') => shapes.push({type:'line',xref:'x',yref:'y',x0,y0,x1,y1,
        line:{color:'#64748b',width:1,dash}});
    line(pair.tdi.x+0.14,low,bracketX,low);
    line(pair.tdts.x+0.14,high,bracketX,high);
    if(pair.wrapsCycle)line(pair.tdts.x,pair.tdts.energy,pair.tdts.x,high);
    if(Math.abs(high-low)>payload.tieTolerance)annotations.push({xref:'x',yref:'y',axref:'x',ayref:'y',
        x:bracketX,y:high,ax:bracketX,ay:low,text:'',showarrow:true,arrowside:'end+start',
        arrowhead:2,startarrowhead:2,arrowwidth:1.6,arrowcolor:'#111827'});
    else line(bracketX-.12,low,bracketX+.12,low,'solid');
    annotations.push({xref:'x',yref:'y',x:bracketX+.18,y:(low+high)/2,xanchor:'left',
        text:`<b>δE<sub>min</sub> = ${result.deltaE.toFixed(2)} eV</b>`,showarrow:false,
        font:{size:19},bgcolor:'rgba(255,255,255,0.96)',borderpad:4});
    const note=(text,x,y,extra={}) => annotations.push({xref:'paper',yref:'paper',x,y,text,
        showarrow:false,xanchor:'left',yanchor:'top',align:'left',font:{size:14,color:'#475569'},...extra});
    note('Capped(Sᵢ) → Pure(Sᵢ) → Capped(Sⱼ)',0,1.15,{font:{size:16,color:'#475569'}});
    note(`${payload.paths.length.toLocaleString('en-US')} paths · Minimum-span pathway: ${result.path.path_id}`+
        (result.ties>1?` · ${result.ties} co-optimal; one shown`:''),0,1.075);
    note(`ΔμMn(OH)<sub>2</sub> = ${muMn.toFixed(2)} eV<br>ΔμH<sub>2</sub>O = ${muH2O.toFixed(2)} eV<br>ΔμOA = ${muOA.toFixed(2)} eV`,
        .985,.15,{xanchor:'right',bgcolor:'rgba(255,255,255,0.96)',borderpad:7,font:{size:13,color:'#475569'}});
    note('cap = carboxyl-capped · TDTS* = relaxed-state proxy'+
        (pair.wrapsCycle?' · Span upper guide: G(TDTS*) + ΔGᵣ':''),0,-.29,{font:{size:12,color:'#64748b'}});
    const allY=[...nodes.values()].map(e=>e.energy);allY.push(low,high);
    const ymin=Math.min(...allY),ymax=Math.max(...allY),pad=Math.max(1,(ymax-ymin)*.18);
    const layout={width:1360,height:680,paper_bgcolor:'white',plot_bgcolor:'white',
        font:{family:'Arial, Helvetica, sans-serif',size:17,color:'#20252b'},
        title:{text:`Capped-start Growth Pathways on Mn₃O₄ (${payload.facetLabel})`,x:.065,xanchor:'left',y:.98,font:{size:25}},
        margin:{l:95,r:35,t:140,b:145},
        xaxis:{title:{text:'Expanded reaction step',standoff:15},range:[-.65,extentX+Math.max(4.2,extentX*.3)],
            tickmode:'array',tickvals:Array.from({length:Math.floor(extentX/(extentX>18?2:1))+1},(_,i)=>i*(extentX>18?2:1)),
            showgrid:false,zeroline:false,showline:true,linecolor:'#374151',ticks:'outside'},
        yaxis:{title:{text:'Relative free energy, ΔG (eV)',standoff:16},range:[ymin-pad,ymax+pad],
            showgrid:false,zeroline:false,showline:true,linecolor:'#374151',ticks:'outside',nticks:6},
        legend:{orientation:'h',x:.5,xanchor:'center',y:-.19,yanchor:'top',font:{size:15},itemsizing:'constant'},
        hovermode:'closest',annotations,shapes};
    return {data,layout,result,muMn,muH2O,muOA,
        backgroundNodeCount:backgroundNodes.length,backgroundEdgeCount:backgroundEdges.length};
}

const graph=document.getElementById('poster-001');
const inputs=['mu-mn','mu-h2o','mu-oa'].map(role=>document.querySelector(`[data-role="${role}"]`));
let lastRender=Promise.resolve(), pendingFrame=null;
function update(){
    const values=inputs.map(el=>Number(el.value));
    inputs.forEach((el,i)=>document.querySelector(`[data-role="${el.dataset.role}-value"]`).textContent=`${values[i].toFixed(2)} eV`);
    const fig=makeCappedPoster(...values);
    window.posterModel=fig;
    lastRender=lastRender.then(()=>Plotly.react(graph,fig.data,fig.layout,config));
    window.posterReady=lastRender;
    return lastRender;
}
inputs.forEach(input=>input.addEventListener('input',()=>{
    if(pendingFrame!==null)cancelAnimationFrame(pendingFrame);
    pendingFrame=requestAnimationFrame(()=>{pendingFrame=null;update();});
}));
for(const button of document.querySelectorAll('[data-export]'))button.addEventListener('click',async()=>{
    button.disabled=true;
    const status=document.getElementById('export-status');
    try{
        if(pendingFrame!==null){cancelAnimationFrame(pendingFrame);pendingFrame=null;}
        await update();
        const format=button.dataset.export;
        status.textContent='Preparing figure…';
        await Plotly.downloadImage(graph,{format,width:1360,height:680,scale:format==='png'?3:1,
            filename:`capped_start_${payload.facetLabel}_mu_${inputs.map(el=>Number(el.value).toFixed(2)).join('_')}`});
        status.textContent=format==='svg'?'SVG exported.':'PNG exported (4080 × 2040).';
    }catch(error){status.textContent=`Export failed: ${error.message}`;}
    finally{button.disabled=false;}
});
update();
