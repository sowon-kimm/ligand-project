(function () {
    "use strict";

    const payload = __PAYLOAD_JSON__;
    const graph = document.getElementById("{plot_id}");
    const stateByName = new Map(payload.states.map((state) => [state.state, state]));

    const controls = document.createElement("div");
    controls.className = "chemical-potential-controls";
    controls.innerHTML = `
        <label>
            <span>Δμ Mn(OH)<sub>2</sub></span>
            <input data-role="mu-mn" type="range" min="${payload.muMin}" max="${payload.muMax}" step="0.05" value="0">
            <output data-role="mu-mn-value">0.00 eV</output>
        </label>
        <label>
            <span>Δμ H<sub>2</sub>O</span>
            <input data-role="mu-h2o" type="range" min="${payload.muMin}" max="${payload.muMax}" step="0.05" value="0">
            <output data-role="mu-h2o-value">0.00 eV</output>
        </label>
    `;

    const style = document.createElement("style");
    style.textContent = `
        .chemical-potential-controls {
            box-sizing: border-box;
            width: 1400px;
            padding: 16px 90px 4px;
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            gap: 36px;
            background: #fff;
            color: #111;
            font: 14px Arial, sans-serif;
        }
        .chemical-potential-controls label {
            display: grid;
            grid-template-columns: 125px minmax(220px, 1fr) 72px;
            align-items: center;
            gap: 12px;
        }
        .chemical-potential-controls label > span {
            font-weight: 600;
            white-space: nowrap;
        }
        .chemical-potential-controls input[type="range"] {
            width: 100%;
            accent-color: #2563eb;
        }
        .chemical-potential-controls output {
            font-variant-numeric: tabular-nums;
            text-align: right;
            white-space: nowrap;
        }
    `;
    document.head.appendChild(style);
    graph.parentElement.insertBefore(controls, graph);

    const muMnInput = controls.querySelector('[data-role="mu-mn"]');
    const muH2OInput = controls.querySelector('[data-role="mu-h2o"]');
    const muMnValue = controls.querySelector('[data-role="mu-mn-value"]');
    const muH2OValue = controls.querySelector('[data-role="mu-h2o-value"]');

    function correctedEnergy(state, energyKey, muMn, muH2O) {
        return state[energyKey] - state.O * muMn + (state.C + 0.5 * state.X) * muH2O;
    }

    function linkY(ligand, muMn, muH2O) {
        const values = [];
        for (const state of payload.states) {
            values.push(
                correctedEnergy(state, "pure", muMn, muH2O),
                correctedEnergy(state, ligand, muMn, muH2O),
                null,
            );
        }
        return values;
    }

    function reactionY(ligand, transition, muMn, muH2O) {
        const values = [];
        for (const [fromName, toName, edgeTransition] of payload.edges) {
            if (edgeTransition !== transition) {
                continue;
            }
            values.push(
                correctedEnergy(stateByName.get(fromName), ligand, muMn, muH2O),
                correctedEnergy(stateByName.get(toName), "pure", muMn, muH2O),
                null,
            );
        }
        return values;
    }

    function levelY(energyKey, muMn, muH2O) {
        const values = [];
        for (const state of payload.states) {
            const energy = correctedEnergy(state, energyKey, muMn, muH2O);
            values.push(energy, energy, null);
        }
        return values;
    }

    function markerY(energyKey, muMn, muH2O) {
        return payload.states.map((state) => correctedEnergy(state, energyKey, muMn, muH2O));
    }

    function expandedEntries(path, ligand, muMn, muH2O) {
        const entries = [];
        for (const name of path.states) {
            const state = stateByName.get(name);
            entries.push(
                {
                    state: name,
                    label: `pure(${name})`,
                    x: state.pure_x,
                    energy: correctedEnergy(state, "pure", muMn, muH2O),
                },
                {
                    state: name,
                    label: `${ligand}-capped(${name})`,
                    x: state.capping_x,
                    energy: correctedEnergy(state, ligand, muMn, muH2O),
                },
            );
        }
        return entries;
    }

    function energeticSpan(path, ligand, muMn, muH2O) {
        const entries = expandedEntries(path, ligand, muMn, muH2O);
        const reactionEnergy = entries[entries.length - 1].energy - entries[0].energy;
        const candidates = entries.slice(0, -1);
        let maximumSpan = -Infinity;
        let pairs = [];

        candidates.forEach((tdts, tdtsIndex) => {
            candidates.forEach((tdi, tdiIndex) => {
                const correction = tdtsIndex < tdiIndex ? reactionEnergy : 0;
                const span = tdts.energy - tdi.energy + correction;
                if (span > maximumSpan + payload.tieTolerance) {
                    maximumSpan = span;
                    pairs = [{tdi, tdts}];
                } else if (Math.abs(span - maximumSpan) <= payload.tieTolerance) {
                    pairs.push({tdi, tdts});
                }
            });
        });
        return {deltaE: maximumSpan, pairs};
    }

    function minimumPath(ligand, muMn, muH2O) {
        let minimumSpan = Infinity;
        let representative = null;
        let representativePair = null;
        let cooptimalCount = 0;

        for (const path of payload.allPaths) {
            const span = energeticSpan(path, ligand, muMn, muH2O);
            if (span.deltaE < minimumSpan - payload.tieTolerance) {
                minimumSpan = span.deltaE;
                representative = path;
                representativePair = span.pairs[0];
                cooptimalCount = 1;
            } else if (Math.abs(span.deltaE - minimumSpan) <= payload.tieTolerance) {
                cooptimalCount += 1;
            }
        }
        return {
            path: representative,
            deltaE: minimumSpan,
            cooptimalCount,
            pair: representativePair,
        };
    }

    function sortedStates(names) {
        return Array.from(names).sort((left, right) => {
            const a = stateByName.get(left);
            const b = stateByName.get(right);
            return a.coordinate - b.coordinate || left.localeCompare(right);
        });
    }

    function selectedLinkXY(names, ligand, muMn, muH2O) {
        const x = [];
        const y = [];
        for (const name of names) {
            const state = stateByName.get(name);
            x.push(
                state.pure_x + payload.levelHalfWidth,
                state.capping_x - payload.levelHalfWidth,
                null,
            );
            y.push(
                correctedEnergy(state, "pure", muMn, muH2O),
                correctedEnergy(state, ligand, muMn, muH2O),
                null,
            );
        }
        return {x, y};
    }

    function selectedReactionXY(path, ligand, transition, muMn, muH2O) {
        const x = [];
        const y = [];
        path.mechanisms.forEach((mechanism, index) => {
            if (mechanism !== transition) {
                return;
            }
            const source = stateByName.get(path.states[index]);
            const target = stateByName.get(path.states[index + 1]);
            x.push(
                source.capping_x + payload.levelHalfWidth,
                target.pure_x - payload.levelHalfWidth,
                null,
            );
            y.push(
                correctedEnergy(source, ligand, muMn, muH2O),
                correctedEnergy(target, "pure", muMn, muH2O),
                null,
            );
        });
        return {x, y};
    }

    function selectedLevelXY(names, energyKey, xKey, muMn, muH2O) {
        const x = [];
        const y = [];
        for (const name of names) {
            const state = stateByName.get(name);
            const energy = correctedEnergy(state, energyKey, muMn, muH2O);
            x.push(
                state[xKey] - payload.levelHalfWidth,
                state[xKey] + payload.levelHalfWidth,
                null,
            );
            y.push(energy, energy, null);
        }
        return {x, y};
    }

    function cappedMarkerData(names, ligand, muMn, muH2O) {
        const x = [];
        const y = [];
        const customdata = [];
        for (const name of names) {
            const state = stateByName.get(name);
            x.push(state.capping_x);
            y.push(correctedEnergy(state, ligand, muMn, muH2O));
            customdata.push([
                state.state,
                state.coordinate,
                state[ligand],
                state[ligand] - state.pure,
                state.O,
                state.C,
                state.X,
            ]);
        }
        return {x, y, customdata};
    }

    function pureMarkerData(names, muMn, muH2O) {
        const x = [];
        const y = [];
        const text = [];
        const customdata = [];
        for (const name of names) {
            const state = stateByName.get(name);
            x.push(state.pure_x);
            y.push(correctedEnergy(state, "pure", muMn, muH2O));
            text.push(state.state);
            customdata.push([
                state.state,
                state.coordinate,
                state.pure,
                state.O,
                state.C,
                state.X,
            ]);
        }
        return {x, y, text, customdata};
    }

    function updateMinimumPaths(muMn, muH2O) {
        const results = {};
        const pureStateUnion = new Set();

        for (const ligand of payload.ligands) {
            const result = minimumPath(ligand, muMn, muH2O);
            results[ligand] = result;
            result.path.states.forEach((name) => pureStateUnion.add(name));
            const indices = payload.traceIndices[ligand];
            const link = selectedLinkXY(
                result.path.states, ligand, muMn, muH2O
            );
            Plotly.restyle(
                graph,
                {x: [link.x], y: [link.y]},
                [indices.capping_link],
            );

            for (const transition of payload.transitions) {
                const reaction = selectedReactionXY(
                    result.path, ligand, transition, muMn, muH2O
                );
                Plotly.restyle(
                    graph,
                    {x: [reaction.x], y: [reaction.y]},
                    [indices.reaction[transition]],
                );
            }

            const levels = selectedLevelXY(
                result.path.states, ligand, "capping_x", muMn, muH2O
            );
            Plotly.restyle(
                graph,
                {x: [levels.x], y: [levels.y]},
                [indices.capped_levels],
            );
            const markers = cappedMarkerData(
                result.path.states, ligand, muMn, muH2O
            );
            Plotly.restyle(
                graph,
                {
                    x: [markers.x],
                    y: [markers.y],
                    customdata: [markers.customdata],
                },
                [indices.capped_markers],
            );
        }

        const pureNames = sortedStates(pureStateUnion);
        const pureLevels = selectedLevelXY(
            pureNames, "pure", "pure_x", muMn, muH2O
        );
        const pureMarkers = pureMarkerData(pureNames, muMn, muH2O);
        const sharedIndices = payload.traceIndices[payload.ligands[0]];
        Plotly.restyle(
            graph,
            {x: [pureLevels.x], y: [pureLevels.y]},
            [sharedIndices.pure_levels],
        );
        Plotly.restyle(
            graph,
            {
                x: [pureMarkers.x],
                y: [pureMarkers.y],
                text: [pureMarkers.text],
                customdata: [pureMarkers.customdata],
            },
            [sharedIndices.pure_markers],
        );

        const carboxyl = results.carboxyl;
        const amine = results.amine;
        if (sharedIndices.tdi_marker !== undefined) {
            Plotly.restyle(
                graph,
                {
                    x: [[carboxyl.pair.tdi.x]],
                    y: [[carboxyl.pair.tdi.energy]],
                    text: [[`TDI: ${carboxyl.pair.tdi.state}`]],
                    customdata: [[[carboxyl.pair.tdi.label]]],
                },
                [sharedIndices.tdi_marker],
            );
            Plotly.restyle(
                graph,
                {
                    x: [[carboxyl.pair.tdts.x]],
                    y: [[carboxyl.pair.tdts.energy]],
                    text: [[`TDTS*: ${carboxyl.pair.tdts.state}`]],
                    customdata: [[[carboxyl.pair.tdts.label]]],
                },
                [sharedIndices.tdts_marker],
            );
        }
        const carboxylTie = carboxyl.cooptimalCount > 1
            ? ` (+${carboxyl.cooptimalCount - 1} co-optimal)`
            : "";
        const amineTie = amine.cooptimalCount > 1
            ? ` (+${amine.cooptimalCount - 1} co-optimal)`
            : "";
        const title = "<b>Minimum pure → capped half-step → next pure pathways</b>" +
            `<br><sup>ΔμMn(OH)₂ = ${muMn.toFixed(2)} eV, ` +
            `ΔμH₂O = ${muH2O.toFixed(2)} eV` +
            ` | Carboxyl δE<sub>span</sub> = ${carboxyl.deltaE.toFixed(4)} eV ` +
            `(${carboxyl.path.path_id})${carboxylTie}` +
            ` | Amine δE<sub>span</sub> = ${amine.deltaE.toFixed(4)} eV ` +
            `(${amine.path.path_id})${amineTie}` +
            `<br>Carboxyl TDI = ${carboxyl.pair.tdi.label}, ` +
            `TDTS* = ${carboxyl.pair.tdts.label}` +
            ` | path: ${carboxyl.path.states.join(" → ")}` +
            `<br>Amine path: ${amine.path.states.join(" → ")}</sup>`;
        Plotly.relayout(graph, {
            "title.text": title,
            "yaxis.autorange": false,
        });
    }

    function update() {
        const muMn = Number(muMnInput.value);
        const muH2O = Number(muH2OInput.value);
        muMnValue.value = `${muMn.toFixed(2)} eV`;
        muH2OValue.value = `${muH2O.toFixed(2)} eV`;

        if (payload.layoutMode === "combined-minimum") {
            updateMinimumPaths(muMn, muH2O);
            return;
        }

        for (const ligand of payload.ligands) {
            const indices = payload.traceIndices[ligand];
            Plotly.restyle(graph, {y: [linkY(ligand, muMn, muH2O)]}, [indices.capping_link]);
            for (const transition of payload.transitions) {
                Plotly.restyle(
                    graph,
                    {y: [reactionY(ligand, transition, muMn, muH2O)]},
                    [indices.reaction[transition]],
                );
            }
            Plotly.restyle(
                graph,
                {y: [levelY("pure", muMn, muH2O)]},
                [indices.pure_levels],
            );
            Plotly.restyle(
                graph,
                {y: [levelY(ligand, muMn, muH2O)]},
                [indices.capped_levels],
            );
            Plotly.restyle(
                graph,
                {y: [markerY("pure", muMn, muH2O)]},
                [indices.pure_markers],
            );
            Plotly.restyle(
                graph,
                {y: [markerY(ligand, muMn, muH2O)]},
                [indices.capped_markers],
            );
        }

        const title = payload.layoutMode === "combined"
            ? "<b>Combined carboxyl and amine capped reaction networks</b>" +
                `<br><sup>Solid reaction lines: carboxyl | dashed reaction lines: amine` +
                ` | ΔμMn(OH)₂ = ${muMn.toFixed(2)} eV, ` +
                `ΔμH₂O = ${muH2O.toFixed(2)} eV</sup>`
            : "<b>Chemical-potential dependent ligand-capped reaction networks</b>" +
                `<br><sup>ΔμMn(OH)₂ = ${muMn.toFixed(2)} eV, ` +
                `ΔμH₂O = ${muH2O.toFixed(2)} eV</sup>`;

        Plotly.relayout(graph, {
            "title.text": title,
            "yaxis.autorange": false,
            "yaxis2.autorange": false,
        });
    }

    let pendingFrame = null;
    function scheduleUpdate() {
        if (pendingFrame !== null) {
            cancelAnimationFrame(pendingFrame);
        }
        pendingFrame = requestAnimationFrame(() => {
            pendingFrame = null;
            update();
        });
    }

    muMnInput.addEventListener("input", scheduleUpdate);
    muH2OInput.addEventListener("input", scheduleUpdate);
    update();
})();
