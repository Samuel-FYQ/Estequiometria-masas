/* =========================================================
   app.js (COMPLETO)
   - Reactivos/Productos ILIMITADOS (añadir/quitar filas)
   - Ajuste con coeficientes mínimos
   - Masas molares
   - Cálculo: conocido (masa/moles) -> objetivo (masa/moles)
   - Conversión “modo cuaderno” con fracciones + tachado + colores
   - Formato ES: coma decimal y SIN separador de miles
   ========================================================= */

document.addEventListener("DOMContentLoaded", () => {

  // ===============================
  //  Formato español (coma) SIN miles + parser que acepta coma
  // ===============================
  function fmt(num, decimals) {
    return Number(num).toLocaleString("es-ES", {
      useGrouping: false, // sin 1.000
      minimumFractionDigits: decimals,
      maximumFractionDigits: decimals
    });
  }

  function toNumber(val) {
    if (val === null || val === undefined) return NaN;
    if (typeof val === "number") return val;
    return Number(String(val).trim().replace(/\s+/g, "").replace(",", "."));
  }

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

      break; // hidratos/cargas no soportados
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
    if (freeCols.length === 0) throw new Error("No se puede ajustar la reacción con las especies dadas.");

    const f = freeCols[0];
    const v = Array(n).fill(0);
    v[f] = 1;

    for (let r = 0; r < pivotCols.length; r++) {
      const pc = pivotCols[r];
      let sum = 0;
      for (let c = pc + 1; c < n; c++) sum += M[r][c] * v[c];
      v[pc] = -sum;
    }

    const fracs = v.map(x => toFraction(x, 2000));
    const denoms = fracs.map(([, q]) => q);
    const L = lcmArrayInt(denoms);

    let ints = fracs.map(([p, q]) => p * (L / q));

    const g = gcdArrayInt(ints.map(x => Math.abs(x)));
    ints = ints.map(x => x / g);

    if (!ints.some(x => x > 0)) ints = ints.map(x => -x);
    if (ints.some(x => x === 0)) throw new Error("Ajuste degenerado: revisa las especies introducidas.");

    return ints;
  }

  // ===============================
  //  Conversión “modo cuaderno” (fracciones inline)
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

  // Colores por pareja cancelada
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
  //  DOM: listas dinámicas de especies
  // ===============================
  const reactantsList = document.getElementById("reactants-list");
  const productsList = document.getElementById("products-list");
  const addReactantBtn = document.getElementById("add-reactant");
  const addProductBtn = document.getElementById("add-product");

  // Tarjeta 1
  const btnBalance = document.getElementById("btn-balance");
  const balanceError = document.getElementById("balance-error");
  const balancedInfo = document.getElementById("balanced-info");
  const balancedEquationText = document.getElementById("balanced-equation-text");
  const molarMassList = document.getElementById("molar-mass-list");
  const resetReaction = document.getElementById("reset-reaction");

  // Tarjeta 2
  const knownSpeciesSelect = document.getElementById("known-species");
  const targetSpeciesSelect = document.getElementById("target-species");

  const knownType = document.getElementById("known-type");     // "mass" | "moles"
  const targetType = document.getElementById("target-type");   // "mass" | "moles"

  const knownMassBlock = document.getElementById("known-mass-block");
  const knownMolesBlock = document.getElementById("known-moles-block");
  const knownMassInput = document.getElementById("known-mass");
  const knownMolesInput = document.getElementById("known-moles");

  const btnCalc = document.getElementById("btn-calc");
  const calcError = document.getElementById("calc-error");
  const resultBox = document.getElementById("result-box");
  const resultText = document.getElementById("result-text");

  const btnShowConversion = document.getElementById("btn-show-conversion");
  const conversionBox = document.getElementById("conversion-box");
  const conversionText = document.getElementById("conversion-text");

  const resetMass = document.getElementById("reset-mass");

  let currentSpecies = [];

  // ===============================
  //  UI adaptativa (masa/moles)
  // ===============================
  function hide(el) { if (el) el.classList.add("hidden"); }
  function show(el) { if (el) el.classList.remove("hidden"); }

  function updateKnownUI() {
    hide(knownMassBlock);
    hide(knownMolesBlock);
    if (knownType.value === "mass") show(knownMassBlock);
    if (knownType.value === "moles") show(knownMolesBlock);
  }

  if (knownType) knownType.addEventListener("change", updateKnownUI);
  updateKnownUI();

  // ===============================
  //  Añadir / quitar filas
  // ===============================
  const MAX_PER_SIDE = 10;

  function attachRemoveHandler(btn) {
    btn.addEventListener("click", (e) => {
      const row = e.currentTarget.closest(".rowline");
      if (!row) return;
      const parent = row.parentElement;
      if (!parent) return;

      const rows = parent.querySelectorAll(".rowline");
      if (rows.length <= 1) {
        const input = row.querySelector("input");
        if (input) input.value = "";
        return;
      }
      row.remove();
    });
  }

  function createRow(side) {
    const row = document.createElement("div");
    row.className = "rowline";

    const inp = document.createElement("input");
    inp.type = "text";
    inp.className = "species-input";
    inp.dataset.side = side;
    inp.placeholder = side === "reactivo"
      ? "Ej: H2, Fe2(SO4)3..."
      : "Ej: H2O, CO2...";

    const btn = document.createElement("button");
    btn.type = "button";
    btn.className = "icon-btn remove-species";
    btn.title = "Quitar";
    btn.setAttribute("aria-label", "Quitar");
    btn.textContent = "🗑️";
    attachRemoveHandler(btn);

    row.appendChild(inp);
    row.appendChild(btn);
    return row;
  }

  function addRow(side) {
    const parent = side === "reactivo" ? reactantsList : productsList;
    if (!parent) return;
    const count = parent.querySelectorAll(".rowline").length;
    if (count >= MAX_PER_SIDE) return;
    parent.appendChild(createRow(side));
  }

  if (addReactantBtn) addReactantBtn.addEventListener("click", () => addRow("reactivo"));
  if (addProductBtn) addProductBtn.addEventListener("click", () => addRow("producto"));

  // Para botones iniciales del HTML
  document.querySelectorAll(".remove-species").forEach(btn => {
    if (!btn.dataset.bound) {
      attachRemoveHandler(btn);
      btn.dataset.bound = "1";
    }
  });

  // ===============================
  //  Leer especies desde UI (N)
  // ===============================
  function readSpeciesFromUI() {
    const inputs = Array.from(document.querySelectorAll(".species-input"));
    const reactivos = [];
    const productos = [];

    for (const inp of inputs) {
      const side = inp.dataset.side;
      const val = normalizeFormula(inp.value.trim());
      if (!val) continue;
      if (side === "reactivo") reactivos.push(val);
      if (side === "producto") productos.push(val);
    }

    const species = [];
    reactivos.forEach((f, idx) => species.push({ id: `r${idx + 1}`, formula: f, side: "reactivo" }));
    productos.forEach((f, idx) => species.push({ id: `p${idx + 1}`, formula: f, side: "producto" }));

    return { species, reactivosCount: reactivos.length, productosCount: productos.length };
  }

  // ===============================
  //  Ajustar reacción
  // ===============================
  btnBalance.addEventListener("click", () => {
    balanceError.textContent = "";
    balancedInfo.classList.add("hidden");

    knownSpeciesSelect.innerHTML = '<option value="">— Selecciona —</option>';
    targetSpeciesSelect.innerHTML = '<option value="">— Selecciona —</option>';
    currentSpecies = [];

    const { species, reactivosCount, productosCount } = readSpeciesFromUI();

    if (reactivosCount < 1 || productosCount < 1) {
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
        li.textContent = `${s.formula}: ${fmt(s.molarMass, 3)} g/mol`;
        molarMassList.appendChild(li);
      });

      enriched.forEach(s => {
        const o1 = document.createElement("option");
        o1.value = s.id;
        o1.textContent = s.formula;
        knownSpeciesSelect.appendChild(o1);

        const o2 = document.createElement("option");
        o2.value = s.id;
        o2.textContent = s.formula;
        targetSpeciesSelect.appendChild(o2);
      });

      balancedInfo.classList.remove("hidden");
    } catch (err) {
      balanceError.textContent = err.message || "Error al ajustar la reacción.";
    }
  });

  // ===============================
  //  Cálculo: masa/moles -> masa/moles
  // ===============================
  btnCalc.addEventListener("click", () => {
    calcError.textContent = "";
    resultBox.classList.add("hidden");
    conversionBox.classList.add("hidden");
    btnShowConversion.classList.add("hidden");
    conversionText.innerHTML = "";

    if (!currentSpecies.length) {
      calcError.textContent = "Primero ajusta la reacción.";
      return;
    }

    const knownId = knownSpeciesSelect.value;
    const targetId = targetSpeciesSelect.value;

    if (!knownId || !targetId) {
      calcError.textContent = "Selecciona ambas sustancias.";
      return;
    }
    if (knownId === targetId) {
      calcError.textContent = "Las sustancias deben ser distintas.";
      return;
    }

    const known = currentSpecies.find(s => s.id === knownId);
    const target = currentSpecies.find(s => s.id === targetId);
    if (!known || !target) {
      calcError.textContent = "Error interno: selección inválida.";
      return;
    }

    // ---- dato conocido -> nKnown
    let nKnown;
    const leftPieces = [];

    if (knownType.value === "mass") {
      const mKnown = toNumber(knownMassInput.value);
      if (!isFinite(mKnown) || mKnown <= 0) {
        calcError.textContent = "Introduce una masa válida.";
        return;
      }
      nKnown = mKnown / known.molarMass;

      leftPieces.push(`<span class="factor">${termHTML(`${fmt(mKnown, 3)}`, `g ${known.formula}`)}</span>`);
      leftPieces.push(`<span class="times">×</span>`);
      leftPieces.push(
        fractionHTML(`1`, `mol ${known.formula}`, `${fmt(known.molarMass, 3)}`, `g ${known.formula}`)
      );
    } else {
      const n = toNumber(knownMolesInput.value);
      if (!isFinite(n) || n <= 0) {
        calcError.textContent = "Introduce moles válidos.";
        return;
      }
      nKnown = n;

      leftPieces.push(`<span class="factor">${termHTML(`${fmt(nKnown, 6)}`, `mol ${known.formula}`)}</span>`);
    }

    // ---- estequiometría
    const stoich = fractionHTML(
      `${target.coeff}`, `mol ${target.formula}`,
      `${known.coeff}`, `mol ${known.formula}`
    );

    const nTarget = nKnown * (target.coeff / known.coeff);

    // ---- salida
    let resultStr = "";
    let rightHTML = "";
    const tailPieces = [];

    if (targetType.value === "moles") {
      resultStr = `${fmt(nTarget, 6)} mol de ${target.formula}`;
      rightHTML = `<span class="factor result">${termHTML(`${fmt(nTarget, 6)}`, `mol ${target.formula}`)}</span>`;
    } else {
      const mTarget = nTarget * target.molarMass;
      resultStr = `${fmt(mTarget, 3)} g de ${target.formula}`;
      rightHTML = `<span class="factor result">${termHTML(`${fmt(mTarget, 3)}`, `g ${target.formula}`)}</span>`;

      tailPieces.push(`<span class="times">×</span>`);
      tailPieces.push(
        fractionHTML(`${fmt(target.molarMass, 3)}`, `g ${target.formula}`, `1`, `mol ${target.formula}`)
      );
    }

    // Resultado (solo resultado)
    resultText.textContent = resultStr;
    resultBox.classList.remove("hidden");

    // Conversión bonita
    conversionText.innerHTML = `
      <span class="conversion-line">
        <span class="conv-left">
          ${leftPieces.join("")}
          <span class="times">×</span>
          ${stoich}
          ${tailPieces.join("")}
        </span>

        <span class="times eq">=</span>

        <span class="conv-right">
          ${rightHTML}
        </span>
      </span>
    `;

    cancelUnitsIn(conversionBox);

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
    function resetList(listEl, side) {
      if (!listEl) return;
      while (listEl.querySelectorAll(".rowline").length > 2) listEl.lastElementChild.remove();
      while (listEl.querySelectorAll(".rowline").length < 2) listEl.appendChild(createRow(side));
      listEl.querySelectorAll("input.species-input").forEach(i => (i.value = ""));
      listEl.querySelectorAll(".remove-species").forEach(btn => {
        if (!btn.dataset.bound) {
          attachRemoveHandler(btn);
          btn.dataset.bound = "1";
        }
      });
    }

    resetList(reactantsList, "reactivo");
    resetList(productsList, "producto");

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
    calcError.textContent = "";
    resultBox.classList.add("hidden");
    conversionBox.classList.add("hidden");
    btnShowConversion.classList.add("hidden");
    conversionText.innerHTML = "";

    knownSpeciesSelect.value = "";
    targetSpeciesSelect.value = "";

    if (knownMassInput) knownMassInput.value = "";
    if (knownMolesInput) knownMolesInput.value = "";
  });

});
