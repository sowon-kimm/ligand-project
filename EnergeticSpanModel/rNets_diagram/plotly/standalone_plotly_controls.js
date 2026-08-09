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
            <input data-role="mu-mn" type="range" min="-4" max="0" step="0.05" value="0">
            <output data-role="mu-mn-value">0.00 eV</output>
        </label>
        <label>
            <span>Δμ H<sub>2</sub>O</span>
            <input data-role="mu-h2o" type="range" min="-4" max="0" step="0.05" value="0">
            <output data-role="mu-h2o-value">0.00 eV</output>
        </label>
    `;

    const style = document.createElement("style");
    style.textContent = `
        .chemical-potential-controls {
            box-sizing: border-box;
            width: 1200px;
            padding: 16px 85px 4px;
            display: grid;
            grid-template-columns: repeat(2, minmax(0, 1fr));
            gap: 36px;
            background: #fff;
            color: #111;
            font: 14px Arial, sans-serif;
        }
        .chemical-potential-controls label {
            display: grid;
            grid-template-columns: 120px minmax(220px, 1fr) 72px;
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

    function correctedEnergy(state, muMn, muH2O) {
        return state.E0 - state.O * muMn + (state.C + 0.5 * state.X) * muH2O;
    }

    function sortStateNames(names) {
        return Array.from(names).sort((left, right) => {
            const coordinateDelta = stateByName.get(left).coordinate - stateByName.get(right).coordinate;
            return coordinateDelta || left.localeCompare(right);
        });
    }

    function updateNetwork(muMn, muH2O) {
        const allowedEdges = payload.edges.filter(([fromName, toName]) => {
            const fromEnergy = correctedEnergy(stateByName.get(fromName), muMn, muH2O);
            const toEnergy = correctedEnergy(stateByName.get(toName), muMn, muH2O);
            return toEnergy - fromEnergy <= payload.barrierThreshold;
        });

        const reachable = new Set(["S000"]);
        let changed = true;
        while (changed) {
            changed = false;
            for (const [fromName, toName] of allowedEdges) {
                if (reachable.has(fromName) && !reachable.has(toName)) {
                    reachable.add(toName);
                    changed = true;
                }
            }
        }

        const reachableEdges = allowedEdges.filter(([fromName]) => reachable.has(fromName));
        payload.transitions.forEach((transition, traceIndex) => {
            const x = [];
            const y = [];
            for (const [fromName, toName, edgeTransition] of reachableEdges) {
                if (edgeTransition !== transition) {
                    continue;
                }
                const fromState = stateByName.get(fromName);
                const toState = stateByName.get(toName);
                x.push(
                    fromState.coordinate + payload.levelHalfWidth,
                    toState.coordinate - payload.levelHalfWidth,
                    null,
                );
                y.push(
                    correctedEnergy(fromState, muMn, muH2O),
                    correctedEnergy(toState, muMn, muH2O),
                    null,
                );
            }
            Plotly.restyle(graph, {x: [x], y: [y]}, [traceIndex]);
        });

        const levelX = [];
        const levelY = [];
        const markerY = [];
        for (const state of payload.states) {
            const energy = correctedEnergy(state, muMn, muH2O);
            levelX.push(
                state.coordinate - payload.levelHalfWidth,
                state.coordinate + payload.levelHalfWidth,
                null,
            );
            levelY.push(energy, energy, null);
            markerY.push(energy);
        }
        Plotly.restyle(graph, {x: [levelX], y: [levelY]}, [4]);
        Plotly.restyle(graph, {y: [markerY]}, [5]);
        Plotly.relayout(graph, {
            "title.text": "Chemical-potential dependent free-energy diagram" +
                `<br><sup>ΔμMn(OH)2 = ${muMn.toFixed(2)} eV, ` +
                `ΔμH2O = ${muH2O.toFixed(2)} eV</sup>`,
            "yaxis.autorange": false,
        });
    }

    function energeticSpan(pathInfo, muMn, muH2O) {
        const names = pathInfo.states;
        const energies = names.map((name) => correctedEnergy(stateByName.get(name), muMn, muH2O));
        const dGr = energies[energies.length - 1] - energies[0];
        const candidateNames = names.slice(0, -1);
        const candidateEnergies = energies.slice(0, -1);
        let maxSpan = -Infinity;
        let pairs = [];

        candidateEnergies.forEach((Gi, i) => {
            candidateEnergies.forEach((Gj, j) => {
                const correction = i < j ? dGr : 0;
                const span = Gi - Gj + correction;
                if (span > maxSpan + payload.tieTolerance) {
                    maxSpan = span;
                    pairs = [[candidateNames[j], candidateNames[i]]];
                } else if (Math.abs(span - maxSpan) <= payload.tieTolerance) {
                    pairs.push([candidateNames[j], candidateNames[i]]);
                }
            });
        });

        return {
            path_id: pathInfo.path_id,
            states: pathInfo.states,
            mechanisms: pathInfo.mechanisms,
            deltaE: maxSpan,
            pairs,
        };
    }

    function updateBestPath(muMn, muH2O) {
        const results = payload.allPaths.map((path) => energeticSpan(path, muMn, muH2O));
        const minimumDeltaE = Math.min(...results.map((result) => result.deltaE));
        const bestPaths = results
            .filter((result) => Math.abs(result.deltaE - minimumDeltaE) <= payload.tieTolerance)
            .sort((left, right) => left.path_id.localeCompare(right.path_id));
        const selectedPaths = payload.showAllCooptimal ? bestPaths : [bestPaths[0]];

        const selectedEdges = Object.fromEntries(payload.transitions.map((name) => [name, new Map()]));
        const selectedStates = new Set();
        const tdiStates = new Set();
        const tdtsStates = new Set();

        for (const result of selectedPaths) {
            result.states.forEach((name) => selectedStates.add(name));
            result.mechanisms.forEach((mechanism, index) => {
                const edge = [result.states[index], result.states[index + 1]];
                selectedEdges[mechanism].set(edge.join("\u0000"), edge);
            });
            result.pairs.forEach(([tdi, tdts]) => {
                tdiStates.add(tdi);
                tdtsStates.add(tdts);
            });
        }

        payload.transitions.forEach((transition, traceIndex) => {
            const edges = Array.from(selectedEdges[transition].values()).sort((left, right) => {
                return left[0].localeCompare(right[0]) || left[1].localeCompare(right[1]);
            });
            const x = [];
            const y = [];
            for (const [fromName, toName] of edges) {
                const fromState = stateByName.get(fromName);
                const toState = stateByName.get(toName);
                x.push(
                    fromState.coordinate + payload.levelHalfWidth,
                    toState.coordinate - payload.levelHalfWidth,
                    null,
                );
                y.push(
                    correctedEnergy(fromState, muMn, muH2O),
                    correctedEnergy(toState, muMn, muH2O),
                    null,
                );
            }
            Plotly.restyle(graph, {x: [x], y: [y]}, [traceIndex]);
        });

        const sortedStates = sortStateNames(selectedStates);
        const levelX = [];
        const levelY = [];
        const markerX = [];
        const markerY = [];
        for (const name of sortedStates) {
            const state = stateByName.get(name);
            const energy = correctedEnergy(state, muMn, muH2O);
            levelX.push(
                state.coordinate - payload.levelHalfWidth,
                state.coordinate + payload.levelHalfWidth,
                null,
            );
            levelY.push(energy, energy, null);
            markerX.push(state.coordinate);
            markerY.push(energy);
        }
        Plotly.restyle(graph, {x: [levelX], y: [levelY]}, [4]);
        Plotly.restyle(graph, {x: [markerX], y: [markerY], text: [sortedStates]}, [5]);

        const sortedTdi = sortStateNames(tdiStates);
        const sortedTdts = sortStateNames(tdtsStates);
        Plotly.restyle(graph, {
            x: [sortedTdi.map((name) => stateByName.get(name).coordinate)],
            y: [sortedTdi.map((name) => correctedEnergy(stateByName.get(name), muMn, muH2O))],
            text: [sortedTdi.map((name) => `TDI: ${name}`)],
        }, [6]);
        Plotly.restyle(graph, {
            x: [sortedTdts.map((name) => stateByName.get(name).coordinate)],
            y: [sortedTdts.map((name) => correctedEnergy(stateByName.get(name), muMn, muH2O))],
            text: [sortedTdts.map((name) => `TDTS: ${name}`)],
        }, [7]);

        const pathText = bestPaths.length === 1
            ? bestPaths[0].path_id
            : `${bestPaths[0].path_id} (+ ${bestPaths.length - 1} co-optimal)`;
        const sequence = bestPaths[0].states.join(" → ");
        Plotly.relayout(graph, {
            "title.text": "<b>Minimum energetic-span pathway</b>" +
                `<br><sup>ΔμMn(OH)₂ = ${muMn.toFixed(2)} eV, ` +
                `ΔμH₂O = ${muH2O.toFixed(2)} eV` +
                ` | δEmin = ${minimumDeltaE.toFixed(4)} eV | ${pathText}` +
                `<br>${sequence}</sup>`,
            "yaxis.autorange": false,
        });
    }

    let pendingFrame = null;
    function update() {
        const muMn = Number(muMnInput.value);
        const muH2O = Number(muH2OInput.value);
        muMnValue.value = `${muMn.toFixed(2)} eV`;
        muH2OValue.value = `${muH2O.toFixed(2)} eV`;
        if (payload.mode === "network") {
            updateNetwork(muMn, muH2O);
        } else {
            updateBestPath(muMn, muH2O);
        }
    }

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
