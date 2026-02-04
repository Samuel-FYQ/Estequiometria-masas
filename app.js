/* ===============================
   app.js (COMPLETO)
   - Ajuste con coeficientes mínimos
   - Masas molares
   - Cálculo de masas (sin detalles en el recuadro)
   - Conversión “modo cuaderno” con fracciones + tachado + colores por pareja
   =============================== */

document.addEventListener("DOMContentLoaded", () => {
  // ===============================
  //  Normalizar subíndices Unicode
  // ===============================
  function normalizeFormula(str) {
    return (str || "")
      .replace(/₀/g, "0").replace(/₁/g, "1").replace(/₂/g, "2")
      .replace(/₃/g, "3").replace(/₄/g, "4").replace(/₅/g, "5")
      .replace(/₆/g, "6").replace(/₇/g, "7").replace(/₈/g, "8")
      .replace(/₉/g, "9");
  }

  // ===============================
  //  Tabla periódica (masas atómicas)
  // ===============================
  const ATOMIC_MASSES = {
    H: 1.008, He: 4.0026, Li: 6.94, Be: 9.0122, B: 10.81, C: 12.011, N: 14.007, O: 15.999,
    F: 18.998, Ne: 20.18, Na: 22.99, Mg: 24.305, Al: 26.982, Si: 28.085, P: 30.974, S: 32.06,
    Cl: 35.45, Ar: 39.948, K: 39.098, Ca: 40.078, Sc: 44.956, Ti: 47.867, V: 50.942, Cr: 51.996,
    Mn: 54.938, Fe: 55.845, Co: 58.933, Ni: 58.693, Cu: 63.546, Zn: 65.38, Ga: 69.723, Ge: 72.63,
    As: 74.922, Se: 78.971, Br: 79.904, Kr: 83.798, Rb: 85.468, Sr: 87.62, Ag: 107.8682,
    Cd: 112.414, Sn: 118.71, Sb: 121.76, I: 126.904, Ba: 137.327, Pt: 195.084, Au: 196.967,
    Hg: 200.592, Pb: 207.2
  };

  // ===============================
  //  Enteros: gcd / lcm
  // ===============================
  function gcdInt(a, b) {
    a = Math.abs(Math.trunc(a));
    b = Math.abs(Math.trunc(b));
    while (b !== 0) {
      const t = b;
      b = a % b;
      a = t;
    }
    return a || 1;
  }
  function gcdArrayInt(arr) {
    return arr.reduce((g, v) => gcdInt(g, v), 0) || 1;
  }
  function lcmInt(a, b) {
    a = Math.abs(Math.trunc(a));
    b = Math.abs(Math.trunc(b));
    if (a === 0 || b === 0) return 0;
    return (a / gcdInt(a, b)) * b;
  }
  function lcmArrayInt(arr) {
    return arr.reduce((l, v) => lcmInt(l, v), 1) || 1;
  }

  // Aproxima x como fracción p/q con q <= maxDen (continued fractions)
  function toFraction(x, maxDen = 2000, eps = 1e-12) {
    if (!isFinite(x)) throw new Error("Número no finito en el ajuste.");
    if (Math.abs(x) < eps) return [0, 1];

    const sign = x < 0 ? -1 : 1;
    x = Math.abs(x);

    let a = Math.floor(x);
    let h1 = 1, k1 = 0;
    let h = a, k = 1;

    let frac = x - a;
    while (k <= maxDen && Math.abs(h / k - x) > eps && frac > eps) {
      frac = 1 / frac;
      a = Math.floor(frac);

      const h2 = h1, k2 = k1;
      h1 = h; k1 = k;
      h = a * h1 + h2;
      k = a * k1 + k2;

      frac = frac - a;
    }

    return [sign * h, k];
  }

  // ===============================
  //  Parser de fórmulas
  // ===============================
  function parseFormula(formula) {
    formula = normalizeFormula(formula);
    const cleaned = formula.replace(/\s+/g, "");
    const stack = [new Map()];
    let i = 0;

    while (i < cleaned.length) {
      const ch = cleaned[i];

      if (ch === "(") {
        stack.push(new Map());
        i++;
        continue;
      }

      if (ch === ")") {
        i++;
        let num = "";
        while (i < cleaned.length && /[0-9]/.test(cleaned[i])) num += cleaned[i++];
        const mult = num ? parseInt(num, 10) : 1;

        const group = stack.pop();
        const top = stack[stack.length - 1];
        for (const [el, c] of group.entries()) {
          top.set(el, (top.get(el) || 0) + c * mult);
        }
        continue;
      }

      if (/[A-Z]/.test(ch)) {
        let sym = ch;
        i++;
        if (i < cleaned.length && /[a-z]/.test(cleaned[i])) sym += cleaned[i++];

        let num = "";
        while (i < cleaned.length && /[0-9]/.test(cleaned[i])) num += cleaned[i++];
        const count = num ? parseInt(num, 10) : 1;

        const top = stack[stack.length - 1];
        top.set(sym, (top.get(sym) || 0) + count);
        continue;
      }

      // Si hay un carácter raro, paramos (hidratados, cargas, etc.)
      break;
    }

    return stack[0];
  }

  function molarMassFromFormula(f) {
    const comp = parseFormula(f);
    let mass = 0;
    for (const [el, c] of comp.entries()) {
      if (!ATOMIC_MASSES[el]) throw new Error("Elemento no soportado: " + el);
      mass += ATOMIC_MASSES[el] * c;
    }
    return mass;
  }

  // ===============================
  //  Matriz elementos × especies
  // ===============================
  function buildMatrix(species) {
    const elements = new Set();
    const parsed = species.map(s => ({ ...s, comp: parseFormula(s.formula) }));
    parsed.forEach(s => s.comp.forEach((_, el) => elements.add(el)));

    const elems = [...elements];
    return elems.map(el =>
      parsed.map(s => (s.comp.get(el) || 0) * (s.side === "reactivo" ? -1 : 1))
    );
  }

  // ===============================
  //  Espacio nulo → coeficientes mínimos
  // ===============================
  function nullspaceVector(A) {
    const m = A.length;
    const n = A[0].length;

    const M = A.map(row => row.slice());
    const pivotCols = [];
    let row = 0;

    // RREF Gauss-Jordan
    for (let col = 0; col < n && row < m; col++) {
      let pivot = -1;
      for (let r = row; r < m; r++) {
        if (Math.abs(M[r][col]) > 1e-12) { pivot = r; break; }
      }
      if (pivot === -1) continue;

      [M[row], M[pivot]] = [M[pivot], M[row]];

      const pv = M[row][col];
      for (let c = col; c < n; c++) M[row][c] /= pv;

      for (let r = 0; r < m; r++) {
        if (r !== row && Math.abs(M[r][col]) > 1e-12) {
          const f = M[r][col];
          for (let c = col; c < n; c++) M[r][c] -= f * M[row][c];
        }
      }

      pivotCols.push(col);
      row++;
    }

    const pivotSet = new Set(pivotCols);
    const freeCols = [];
    for (let c = 0; c < n; c++) if (!pivotSet.has(c)) freeCols.push(c);

    if (freeCols.length === 0) {
      throw new Error("No se puede ajustar la reacción con las especies dadas.");
    }

    // Vector del espacio nulo (float)
    const f = freeCols[0];
    const v = Array(n).fill(0);
    v[f] = 1;

    for (let r = 0; r < pivotCols.length; r++) {
      const pc = pivotCols[r];
      let sum = 0;
      for (let c = pc + 1; c < n; c++) sum += M[r][c] * v[c];
      v[pc] = -sum;
    }

    // Pasar a fracciones p/q y luego a enteros mínimos
    const fracs = v.map(x => toFraction(x, 2000));
    const denoms = fracs.map(([p, q]) => q);
    const L = lcmArrayInt(denoms);

    let ints = fracs.map(([p, q]) => p * (L / q));

    // Reducir por mcd
    const g = gcdArrayInt(ints.map(x => Math.abs(x)));
    ints = ints.map(x => x / g);

    // Convención: que haya positivos
    if (!ints.some(x => x > 0)) ints = ints.map(x => -x);

    if (ints.some(x => x === 0)) {
      throw new Error("Ajuste degenerado: revisa las especies introducidas.");
    }

    return ints;
  }

  // ===============================
  //  Conversión “modo cuaderno” (INLINE)
  // ===============================
  function unitSpan(unitText) {
    const u = (unitText || "").trim().replace(/\s+/g, " ");
    return `<span class="u" data-u="${u}">${u}</span>`;
  }
  function termHTML(valueText, unitText) {
    return `<span class="term"><span class="v">${valueText}</span> ${unitSpan(unitText)}</span>`;
  }
  function fractionHTML(topValue, topUnit, bottomValue, bottomUnit) {
    return `
      <span class="fraction">
        <span class="top">${termHTML(topValue, topUnit)}</span>
        <span class="bar"></span>
        <span class="bottom">${termHTML(bottomValue, bottomUnit)}</span>
      </span>
    `;
  }

  // Colores por pareja cancelada (clases c0..c7)
  function cancelUnitsIn(conversionBoxEl) {
    if (!conversionBoxEl) return;
    const left = conversionBoxEl.querySelector(".conv-left");
    if (!left) return;

    const tops = Array.from(left.querySelectorAll(".top .u, .factor .u"));
    const bottoms = Array.from(left.querySelectorAll(".bottom .u"));

    const mapTop = new Map();
    const mapBottom = new Map();

    function push(map, el) {
      const key = (el.dataset.u || "").trim();
      if (!key) return;
      if (!map.has(key)) map.set(key, []);
      map.get(key).push(el);
    }

    tops.forEach(el => push(mapTop, el));
    bottoms.forEach(el => push(mapBottom, el));

    let colorIndex = 0;
    for (const [unit, arrTop] of mapTop.entries()) {
      const arrBot = mapBottom.get(unit);
      if (!arrBot || !arrBot.length) continue;

      const k = Math.min(arrTop.length, arrBot.length);
      const cls = `c${colorIndex % 8}`;
      colorIndex++;

      for (let i = 0; i < k; i++) {
        arrTop[i].classList.add("cancel", cls);
        arrBot[i].classList.add("cancel", cls);
      }
    }
  }

  // ===============================
  //  DOM
  // ===============================
  const r1Input = document.getElementById("r1");
  const r2Input = document.getElementById("r2");
  const p1Input = document.getElementById("p1");
  const p2Input = document.getElementById("p2");

  const btnBalance = document.getElementById("btn-balance");
  const balanceError = document.getElementById("balance-error");
  const balancedInfo = document.getElementById("balanced-info");
  const balancedEquationText = document.getElementById("balanced-equation-text");
  const molarMassList = document.getElementById("molar-mass-list");

  const knownSpeciesSelect = document.getElementById("known-species");
  const targetSpeciesSelect = document.getElementById("target-species");
  const knownMassInput = document.getElementById("known-mass");
  const btnCalc = document.getElementById("btn-calc");
  const calcError = document.getElementById("calc-error");
  const resultBox = document.getElementById("result-box");
  const resultText = document.getElementById("result-text");

  const btnShowConversion = document.getElementById("btn-show-conversion");
  const conversionBox = document.getElementById("conversion-box");
  const conversionText = document.getElementById("conversion-text");

  const resetReaction = document.getElementById("reset-reaction");
  const resetMass = document.getElementById("reset-mass");

  let currentSpecies = [];

  // ===============================
  //  Ajustar reacción
  // ===============================
  btnBalance.addEventListener("click", () => {
    balanceError.textContent = "";
    balancedInfo.classList.add("hidden");

    knownSpeciesSelect.innerHTML = '<option value="">— Selecciona —</option>';
    targetSpeciesSelect.innerHTML = '<option value="">— Selecciona —</option>';
    currentSpecies = [];

    const r1 = normalizeFormula(r1Input.value.trim());
    const r2 = normalizeFormula(r2Input.value.trim());
    const p1 = normalizeFormula(p1Input.value.trim());
    const p2 = normalizeFormula(p2Input.value.trim());

    const species = [];
    if (r1) species.push({ id: "r1", formula: r1, side: "reactivo" });
    if (r2) species.push({ id: "r2", formula: r2, side: "reactivo" });
    if (p1) species.push({ id: "p1", formula: p1, side: "producto" });
    if (p2) species.push({ id: "p2", formula: p2, side: "producto" });

    if (!species.some(s => s.side === "reactivo") || !species.some(s => s.side === "producto")) {
      balanceError.textContent = "Introduce al menos un reactivo y un producto.";
      return;
    }

    try {
      const A = buildMatrix(species);
      const coeffs = nullspaceVector(A);

      const enriched = species.map((sp, i) => ({
        ...sp,
        coeff: coeffs[i],
        molarMass: molarMassFromFormula(sp.formula)
      }));

      currentSpecies = enriched;

      const left = enriched
        .filter(s => s.side === "reactivo")
        .map(s => `${s.coeff === 1 ? "" : s.coeff + " "}${s.formula}`)
        .join(" + ");

      const right = enriched
        .filter(s => s.side === "producto")
        .map(s => `${s.coeff === 1 ? "" : s.coeff + " "}${s.formula}`)
        .join(" + ");

      balancedEquationText.textContent = left + " → " + right;

      molarMassList.innerHTML = "";
      enriched.forEach(s => {
        const li = document.createElement("li");
        li.textContent = `${s.formula}: ${s.molarMass.toFixed(3)} g/mol`;
        molarMassList.appendChild(li);
      });

      enriched.forEach(s => {
        const opt1 = document.createElement("option");
        opt1.value = s.id;
        opt1.textContent = s.formula;
        knownSpeciesSelect.appendChild(opt1);

        const opt2 = document.createElement("option");
        opt2.value = s.id;
        opt2.textContent = s.formula;
        targetSpeciesSelect.appendChild(opt2);
      });

      balancedInfo.classList.remove("hidden");
    } catch (err) {
      balanceError.textContent = err.message || "Error al ajustar la reacción.";
    }
  });

  // ===============================
  //  Cálculo de masas
  // ===============================
  btnCalc.addEventListener("click", () => {
    calcError.textContent = "";
    resultBox.classList.add("hidden");
    conversionBox.classList.add("hidden");
    btnShowConversion.classList.add("hidden");

    if (!currentSpecies.length) {
      calcError.textContent = "Primero ajusta la reacción.";
      return;
    }

    const knownId = knownSpeciesSelect.value;
    const targetId = targetSpeciesSelect.value;
    const massStr = knownMassInput.value.trim();

    if (!knownId || !targetId) {
      calcError.textContent = "Selecciona ambas sustancias.";
      return;
    }
    if (knownId === targetId) {
      calcError.textContent = "Las sustancias deben ser distintas.";
      return;
    }
    if (!massStr || isNaN(Number(massStr)) || Number(massStr) <= 0) {
      calcError.textContent = "Introduce una masa válida.";
      return;
    }

    const known = currentSpecies.find(s => s.id === knownId);
    const target = currentSpecies.find(s => s.id === targetId);
    if (!known || !target) {
      calcError.textContent = "Error interno: selección inválida.";
      return;
    }

    const mKnown = Number(massStr);
    const nKnown = mKnown / known.molarMass;
    const nTarget = nKnown * (target.coeff / known.coeff);
    const mTarget = nTarget * target.molarMass;

    // ✅ SOLO RESULTADO (sin detalles)
    resultText.textContent = `${mTarget.toFixed(3)} g de ${target.formula}`;

    // Conversión bonita
    conversionText.innerHTML = `
      <span class="conversion-line">
        <span class="conv-left">
          <span class="factor">${termHTML(`${mKnown}`, `g ${known.formula}`)}</span>
          <span class="times">×</span>
          ${fractionHTML(`1`, `mol ${known.formula}`, `${known.molarMass.toFixed(3)}`, `g ${known.formula}`)}
          <span class="times">×</span>
          ${fractionHTML(`${target.coeff}`, `mol ${target.formula}`, `${known.coeff}`, `mol ${known.formula}`)}
          <span class="times">×</span>
          ${fractionHTML(`${target.molarMass.toFixed(3)}`, `g ${target.formula}`, `1`, `mol ${target.formula}`)}
        </span>

        <span class="times eq">=</span>

        <span class="conv-right">
          <span class="factor result">${termHTML(`${mTarget.toFixed(3)}`, `g ${target.formula}`)}</span>
        </span>
      </span>
    `;

    cancelUnitsIn(conversionBox);

    resultBox.classList.remove("hidden");
    btnShowConversion.textContent = "Mostrar conversión";
    btnShowConversion.classList.remove("hidden");
  });

  // ===============================
  //  Mostrar / ocultar conversión
  // ===============================
  btnShowConversion.addEventListener("click", () => {
    const hidden = conversionBox.classList.contains("hidden");
    if (hidden) {
      conversionBox.classList.remove("hidden");
      btnShowConversion.textContent = "Ocultar conversión";
    } else {
      conversionBox.classList.add("hidden");
      btnShowConversion.textContent = "Mostrar conversión";
    }
  });

  // ===============================
  //  Reinicio tarjeta 1
  // ===============================
  resetReaction.addEventListener("click", () => {
    r1Input.value = "";
    r2Input.value = "";
    p1Input.value = "";
    p2Input.value = "";

    balanceError.textContent = "";
    balancedInfo.classList.add("hidden");

    knownSpeciesSelect.innerHTML = '<option value="">— Selecciona —</option>';
    targetSpeciesSelect.innerHTML = '<option value="">— Selecciona —</option>';

    calcError.textContent = "";
    resultBox.classList.add("hidden");
    conversionBox.classList.add("hidden");
    btnShowConversion.classList.add("hidden");
    conversionText.innerHTML = "";

    currentSpecies = [];
  });

  // ===============================
  //  Reinicio tarjeta 2
  // ===============================
  resetMass.addEventListener("click", () => {
    knownSpeciesSelect.value = "";
    targetSpeciesSelect.value = "";
    knownMassInput.value = "";

    calcError.textContent = "";
    resultBox.classList.add("hidden");
    conversionBox.classList.add("hidden");
    btnShowConversion.classList.add("hidden");
    conversionText.innerHTML = "";
  });
});
