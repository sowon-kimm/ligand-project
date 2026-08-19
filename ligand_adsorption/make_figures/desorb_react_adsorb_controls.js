(function () {
    "use strict";

    const payload = __PAYLOAD_JSON__;
    const graph = document.getElementById("{plot_id}");
    const stateByName = new Map(payload.states.map((state) => [state.state, state]));

    const controls = document.createElement("div");
    controls.className = "dra-controls";
    controls.innerHTML = `
        <label><span>Δμ Mn(OH)<sub>2</sub></span>
            <input data-role="mu-mn" type="range" min="${payload.muMin}" max="${payload.muMax}" step="0.05" value="0">
            <output data-role="mu-mn-value">0.00 eV</output></label>
        <label><span>Δμ H<sub>2</sub>O</span>
            <input data-role="mu-h2o" type="range" min="${payload.muMin}" max="${payload.muMax}" step="0.05" value="0">
            <output data-role="mu-h2o-value">0.00 eV</output></label>`;
    const style = document.createElement("style");
    style.textContent = `
        .dra-controls { box-sizing:border-box; width:1500px; padding:16px 90px 4px;
            display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:36px;
            background:#fff; color:#111; font:14px Arial,sans-serif; }
        .dra-controls label { display:grid; grid-template-columns:125px minmax(220px,1fr) 72px;
            align-items:center; gap:12px; }
        .dra-controls label > span { font-weight:600; white-space:nowrap; }
        .dra-controls input[type="range"] { width:100%; accent-color:#2563eb; }
        .dra-controls output { font-variant-numeric:tabular-nums; text-align:right; white-space:nowrap; }`;
    document.head.appendChild(style);
    graph.parentElement.insertBefore(controls, graph);

    const muMnInput = controls.querySelector('[data-role="mu-mn"]');
    const muH2OInput = controls.querySelector('[data-role="mu-h2o"]');
    const muMnValue = controls.querySelector('[data-role="mu-mn-value"]');
    const muH2OValue = controls.querySelector('[data-role="mu-h2o-value"]');

    function energy(name, key, muMn, muH2O) {
        const state = stateByName.get(name);
        return state[key] - state.O * muMn + (state.C + 0.5 * state.X) * muH2O;
    }

    function expandedEntries(path, ligand, muMn, muH2O) {
        const entries = [];
        for (const name of path.states) {
            const state = stateByName.get(name);
            entries.push(
                {
                    state: name,
                    kind: "pure",
                    label: `pure(${name})`,
                    x: state.pure_x,
                    energy: energy(name, "pure", muMn, muH2O),
                },
                {
                    state: name,
                    kind: "capped",
                    label: `${ligand}-capped(${name})`,
                    x: state.capping_x,
                    energy: energy(name, ligand, muMn, muH2O),
                },
            );
        }
        return entries;
    }

    function energeticSpan(path, ligand, muMn, muH2O) {
        const entries = expandedEntries(path, ligand, muMn, muH2O);
        const reactionEnergy = entries[entries.length - 1].energy - entries[0].energy;
        const candidates = entries.slice(0, -1);
        let maximum = -Infinity;
        let pairs = [];
        candidates.forEach((tdts, i) => {
            candidates.forEach((tdi, j) => {
                const span = tdts.energy - tdi.energy + (i < j ? reactionEnergy : 0);
                if (span > maximum + payload.tieTolerance) {
                    maximum = span;
                    pairs = [{tdi, tdts, tdiIndex: j, tdtsIndex: i}];
                } else if (Math.abs(span - maximum) <= payload.tieTolerance) {
                    pairs.push({tdi, tdts, tdiIndex: j, tdtsIndex: i});
                }
            });
        });
        return {deltaE: maximum, pairs, reactionEnergy};
    }

    function minimumPath(ligand, muMn, muH2O) {
        let minimum = Infinity;
        let representative = null;
        let representativePair = null;
        let pairTies = 0;
        let ties = 0;
        for (const path of payload.paths) {
            const span = energeticSpan(path, ligand, muMn, muH2O);
            if (span.deltaE < minimum - payload.tieTolerance) {
                minimum = span.deltaE;
                representative = path;
                representativePair = span.pairs[0];
                pairTies = span.pairs.length;
                ties = 1;
            } else if (Math.abs(span.deltaE - minimum) <= payload.tieTolerance) {
                ties += 1;
            }
        }
        return {
            path: representative,
            deltaE: minimum,
            ties,
            pair: representativePair,
            pairTies,
        };
    }

    function addLevel(x, y, xs, ys) {
        xs.push(x - payload.levelHalfWidth, x + payload.levelHalfWidth, null);
        ys.push(y, y, null);
    }

    function pathData(path, ligand, muMn, muH2O) {
        const out = {
            desorption: {x: [], y: []}, adsorption: {x: [], y: []},
            reaction: Object.fromEntries(payload.transitions.map((name) => [name, {x: [], y: []}])),
            cappedLevels: {x: [], y: []}, pureLevels: {x: [], y: []},
            capped: {x: [], y: [], text: [], customdata: []},
            pure: {x: [], y: [], text: [], customdata: []},
        };
        for (let i = 0; i < path.states.length; i += 1) {
            const name = path.states[i];
            const state = stateByName.get(name);
            const pureEnergy = energy(name, "pure", muMn, muH2O);
            const cappedEnergy = energy(name, ligand, muMn, muH2O);
            const adsorptionDelta = cappedEnergy - pureEnergy;

            addLevel(state.pure_x, pureEnergy, out.pureLevels.x, out.pureLevels.y);
            addLevel(state.capping_x, cappedEnergy, out.cappedLevels.x, out.cappedLevels.y);
            out.desorption.x.push(
                state.pure_x + payload.levelHalfWidth,
                state.capping_x - payload.levelHalfWidth,
                null,
            );
            out.desorption.y.push(pureEnergy, cappedEnergy, null);
            out.pure.x.push(state.pure_x);
            out.pure.y.push(pureEnergy);
            out.pure.text.push(name);
            if (i === 0) {
                out.pure.customdata.push([name, "pure initial state", 0]);
            } else {
                const source = path.states[i - 1];
                const sourceCapped = energy(source, ligand, muMn, muH2O);
                const mechanism = path.mechanisms[i - 1];
                out.pure.customdata.push([
                    name,
                    `${mechanism}: after desorption + reaction`,
                    pureEnergy - sourceCapped,
                ]);
            }
            out.capped.x.push(state.capping_x);
            out.capped.y.push(cappedEnergy);
            out.capped.text.push(name);
            out.capped.customdata.push([
                name,
                `${ligand} capped half-step`,
                adsorptionDelta,
            ]);
        }

        for (let i = 0; i < path.states.length - 1; i += 1) {
            const source = path.states[i];
            const target = path.states[i + 1];
            const sourceState = stateByName.get(source);
            const targetState = stateByName.get(target);
            const mechanism = path.mechanisms[i];
            out.reaction[mechanism].x.push(
                sourceState.capping_x + payload.levelHalfWidth,
                targetState.pure_x - payload.levelHalfWidth,
                null,
            );
            out.reaction[mechanism].y.push(
                energy(source, ligand, muMn, muH2O),
                energy(target, "pure", muMn, muH2O),
                null,
            );
        }
        return out;
    }

    function restyleLigand(ligand, result, muMn, muH2O) {
        const data = pathData(result.path, ligand, muMn, muH2O);
        const ix = payload.traceIndices[ligand];
        Plotly.restyle(graph, {x: [data.desorption.x], y: [data.desorption.y]}, [ix.desorption]);
        for (const mechanism of payload.transitions) {
            Plotly.restyle(graph, {x: [data.reaction[mechanism].x], y: [data.reaction[mechanism].y]}, [ix.reaction[mechanism]]);
        }
        Plotly.restyle(graph, {x: [data.adsorption.x], y: [data.adsorption.y]}, [ix.adsorption]);
        Plotly.restyle(graph, {x: [data.cappedLevels.x], y: [data.cappedLevels.y]}, [ix.capped_levels]);
        Plotly.restyle(graph, {
            x: [data.capped.x], y: [data.capped.y],
            customdata: [data.capped.customdata],
        }, [ix.capped_markers]);
        Plotly.restyle(graph, {x: [data.pureLevels.x], y: [data.pureLevels.y]}, [ix.pure_levels]);
        Plotly.restyle(graph, {
            x: [data.pure.x], y: [data.pure.y], text: [data.pure.text],
            customdata: [data.pure.customdata],
        }, [ix.pure_markers]);
        if (ligand === "carboxyl") {
            const pair = result.pair;
            Plotly.restyle(graph, {
                x: [[pair.tdi.x]], y: [[pair.tdi.energy]],
                text: [[`TDI: ${pair.tdi.state}`]],
                customdata: [[[pair.tdi.label, ligand]]],
            }, [ix.tdi_marker]);
            Plotly.restyle(graph, {
                x: [[pair.tdts.x]], y: [[pair.tdts.energy]],
                text: [[`TDTS*: ${pair.tdts.state}`]],
                customdata: [[[pair.tdts.label, ligand]]],
            }, [ix.tdts_marker]);
        }
    }

    function update() {
        const muMn = Number(muMnInput.value);
        const muH2O = Number(muH2OInput.value);
        muMnValue.value = `${muMn.toFixed(2)} eV`;
        muH2OValue.value = `${muH2O.toFixed(2)} eV`;
        const results = {};
        let maximumCoordinate = 0;
        for (const ligand of payload.ligands) {
            results[ligand] = minimumPath(ligand, muMn, muH2O);
            for (const name of results[ligand].path.states) {
                maximumCoordinate = Math.max(maximumCoordinate, stateByName.get(name).capping_x);
            }
            restyleLigand(ligand, results[ligand], muMn, muH2O);
        }
        const tieText = (result) => result.ties > 1 ? ` (+${result.ties - 1} co-optimal)` : "";
        const carboxyl = results.carboxyl;
        const amine = results.amine;
        Plotly.relayout(graph, {
            "title.text": "<b>Pure integer → capped half-step → next pure minimum pathways (001)</b>" +
                `<br><sup>ΔμMn(OH)₂ = ${muMn.toFixed(2)} eV, ΔμH₂O = ${muH2O.toFixed(2)} eV` +
                ` | Carboxyl δE<sub>span</sub> = ${carboxyl.deltaE.toFixed(4)} eV (${carboxyl.path.path_id})${tieText(carboxyl)}` +
                ` | Amine δE<sub>span</sub> = ${amine.deltaE.toFixed(4)} eV (${amine.path.path_id})${tieText(amine)}` +
                `<br>Carboxyl TDI = ${carboxyl.pair.tdi.label}, TDTS* = ${carboxyl.pair.tdts.label}` +
                ` | path: ${carboxyl.path.states.join(" → ")}` +
                `<br>Amine path: ${amine.path.states.join(" → ")}</sup>`,
            "xaxis.range": [-0.35, maximumCoordinate + 0.35],
            "yaxis.autorange": false,
        });
    }

    let pending = null;
    function scheduleUpdate() {
        if (pending !== null) cancelAnimationFrame(pending);
        pending = requestAnimationFrame(() => { pending = null; update(); });
    }
    muMnInput.addEventListener("input", scheduleUpdate);
    muH2OInput.addEventListener("input", scheduleUpdate);
    update();
})();
