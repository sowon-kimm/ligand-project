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
                correctedEnergy(stateByName.get(toName), ligand, muMn, muH2O),
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

    function update() {
        const muMn = Number(muMnInput.value);
        const muH2O = Number(muH2OInput.value);
        muMnValue.value = `${muMn.toFixed(2)} eV`;
        muH2OValue.value = `${muH2O.toFixed(2)} eV`;

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
